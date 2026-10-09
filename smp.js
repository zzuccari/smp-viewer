// Parser for Micromeritics 3Flex .SMP files (reverse-engineered, Flex v6.03)
const OMEGA = {"carbon dioxide": 0.2239, "oxygen": 0.0222, "nitrogen": 0.0372, "argon": 0.0, "methane": 0.0114,
               "krypton": 0.0, "hydrogen": -0.219, "ethane": 0.0995, "ethylene": 0.0866, "propane": 0.1523};
// correction of the Tsonopoulos estimate, calibrated against Flex 6.03 reports
const B_SCALE = {"oxygen": 0.9837, "carbon dioxide": 0.9962};
const COLS = ["P (mmHg)", "P (mbar)", "P/Po", "Q (cm3/g STP)", "Q (mmol/g)", "t (min)"];

function utf16(s) {
  const out = new Uint8Array(s.length * 2);
  for (let i = 0; i < s.length; i++) out[2 * i] = s.charCodeAt(i);
  return out;
}

function find(b, pat, from = 0) {
  outer: for (let i = from; i <= b.length - pat.length; i++) {
    for (let j = 0; j < pat.length; j++) if (b[i + j] !== pat[j]) continue outer;
    return i;
  }
  return -1;
}

function rfind(b, pat, end) {
  outer: for (let i = end - pat.length; i >= 0; i--) {
    for (let j = 0; j < pat.length; j++) if (b[i + j] !== pat[j]) continue outer;
    return i;
  }
  return -1;
}

// UTF-16 strings of at least 3 printable characters
function strings(b) {
  const out = [];
  let i = 0;
  while (i < b.length - 1) {
    let j = i, s = "";
    while (j < b.length - 1 && b[j + 1] === 0 && ((b[j] >= 0x20 && b[j] <= 0x7e) || b[j] >= 0xb0)) {
      s += String.fromCharCode(b[j]);
      j += 2;
    }
    if (s.length >= 3) { out.push(s); i = j; } else i++;
  }
  return out;
}

function round(x, n) { return Number(x.toFixed(n)); }

function parseSMP(buf) {
  const d = new Uint8Array(buf);
  if (String.fromCharCode(...d.slice(0, 6)) !== "^vMIC#") throw new Error("Not a Micromeritics SMP file");
  const dv = new DataView(d.buffer, d.byteOffset, d.byteLength);

  // block directory: (id u16, offset u32, length u32) from 0x35
  const blocks = {};
  for (let p = 0x35; ; p += 10) {
    const id = dv.getUint16(p, true);
    if (id === 0) break;
    const o = dv.getUint32(p + 2, true), l = dv.getUint32(p + 6, true);
    blocks[id] = {b: d.subarray(o, o + l), v: new DataView(d.buffer, d.byteOffset + o, l)};
  }
  const f64 = (blk, off) => blk.v.getFloat64(off, true);

  let blk = blocks[301];
  const q = find(blk.b, utf16("Bar Code:")) + 20 + 7 + 5;
  const mass = f64(blk, q), mTube = f64(blk, q + 8), mTubeSample = f64(blk, q + 16);

  // 66-byte records starting at 0xc6, first one is the zero point.
  // Each record stores the total quantity dosed (cm3 STP); the software computes
  // quantity adsorbed = (dosed - gas left in the free space) / mass at report time.
  blk = blocks[303];
  const b = blk.b;
  const raw = [];
  let p = 0xc6;
  for (; b[p] === 0xe0; p += 66)
    raw.push([f64(blk, p + 11), f64(blk, p + 19), f64(blk, p + 27), blk.v.getUint32(p + 35, true)]);
  // the fixed-layout section follows the last record
  const fsAnalysis = f64(blk, p + 314), fsAmbient = f64(blk, p + 322);
  const T = f64(blk, p + 25);
  const g = find(b, utf16(".fpi"));
  const gas = String.fromCharCode(...b.subarray(rfind(b, utf16("\\"), g) + 2, g).filter((_, i) => i % 2 === 0));
  const Tc = f64(blk, g + 28 + 44), Pc = f64(blk, g + 28 + 52);

  // non-ideality of the free-space gas from the second (Tsonopoulos) and
  // third (Orbey-Vera) virial coefficients
  const Tr = T / Tc, w = OMEGA[gas] || 0, R = 62363.6;
  const f0 = 0.1445 - 0.330 / Tr - 0.1385 / Tr ** 2 - 0.0121 / Tr ** 3 - 0.000607 / Tr ** 8;
  const f1 = 0.0637 + 0.331 / Tr ** 2 - 0.423 / Tr ** 3 - 0.008 / Tr ** 8;
  const B = R * Tc / Pc * (f0 + w * f1) * (B_SCALE[gas] || 1);  // cm3/mol
  const C = (R * Tc / Pc) ** 2 * (
    0.01407 + 0.02432 / Tr ** 2.8 - 0.00313 / Tr ** 10.5
    + w * (-0.02676 + 0.0177 / Tr ** 2.8 + 0.040 / Tr ** 3 - 0.003 / Tr ** 6 - 0.00228 / Tr ** 10.5)
  );  // cm6/mol2
  const RT = R * T;
  const k = -B / RT;  // 1/mmHg
  const k2 = k * k - (C - B * B) / RT ** 2;  // 1/mmHg2

  const rows = raw.slice(1).map(([P, Prel, D, t]) => {
    const Q = (D - fsAnalysis * P / 760 * (1 + k * P + k2 * P * P)) / mass;
    return [P, P * 1.333224, Prel, Q, Q / 22.414, t];
  });

  const s301 = strings(blocks[301].b), s303 = strings(b);
  const saved = new Date(dv.getUint32(0x10, true) * 1000).toLocaleString();
  const info =
    `Sample: ${s301[0]}\nMethod: ${s301[5]}\n` +
    `Instrument: ${s303[3]} s/n ${s303[2]}, ${s303[1]}\n` +
    `File saved: ${saved}\n` +
    `Sample mass: ${round(mass, 6)} g (tube ${round(mTube, 6)} g, tube+sample ${round(mTubeSample, 6)} g)\n` +
    `Free space: ambient ${fsAmbient.toFixed(4)} cm3, analysis ${fsAnalysis.toFixed(4)} cm3\n` +
    `Gas: ${gas}, bath temperature ${round(T, 4)} K, B = ${B.toFixed(1)} cm3/mol\n` +
    `Points: ${rows.length}\n\nInstrument log:\n` + (blocks[705] ? strings(blocks[705].b).join("\n") : "");
  return {name: s301[0], rows, info, mass};
}

// rows scaled to a different sample mass (g) and rounded like the Flex report
function formatRows(rows, fileMass, mass) {
  const f = fileMass / mass;
  return rows.map(([P, Pm, Pr, Q, Qm, t]) =>
    [round(P, 6), round(Pm, 6), round(Pr, 9), round(Q * f, 5), round(Qm * f, 6), t]);
}
