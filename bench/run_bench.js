// Benchmark on real stacks: runs the GhostBuster core on FITS files, the way the automated mode does,
// and writes for each image a contact sheet (raw | corrected, brightest stars) and per-star metrics.
//
// Usage:   node bench/run_bench.js <file or folder> [...] [--out bench-results] [--stars 8] [--all]
//   --out    output folder (default: bench-results); results go to <out>/v<version>/
//   --stars  number of stars on each contact sheet (default: 8)
//   --all    also correct the uncertain (orange) stars, which the dialog leaves unselected by default
//
// Requires Node.js >= 14, no other dependency. Reads 2D (mono) or 3D (RGB) FITS, BITPIX 8, 16, 32,
// -32 or -64. Integer data is scaled to [0, 1] like PixInsight does. The FITS headers are only used
// for the pixel scale; nothing from them is written to the results.
"use strict";
const fs = require("fs");
const path = require("path");
const zlib = require("zlib");
const { loadCore, source } = require("../tests/load.js");

// ---------------------------------------------------------------- FITS reader
function readFits(file) {
   const buf = fs.readFileSync(file);
   const kw = [], hdr = {};
   let pos = 0, end = false;
   while (!end && pos + 2880 <= buf.length) {
      for (let i = 0; i < 36; ++i) {
         const card = buf.toString("latin1", pos + i * 80, pos + i * 80 + 80);
         const name = card.slice(0, 8).trim();
         if (name === "END") { end = true; break; }
         if (card[8] !== "=") continue;
         let v = card.slice(10);
         const q = v.indexOf("'");
         if (q >= 0) { const q2 = v.indexOf("'", q + 1); v = v.slice(q + 1, q2 < 0 ? undefined : q2).trim(); }
         else v = v.split("/")[0].trim();
         hdr[name] = v; kw.push({ name, value: v });
      }
      pos += 2880;
   }
   if (!end) throw new Error("no END card in the primary header");
   const bitpix = parseInt(hdr.BITPIX, 10), naxis = parseInt(hdr.NAXIS, 10);
   if (naxis < 2 || naxis > 3) throw new Error("NAXIS " + naxis + " not supported");
   const W = parseInt(hdr.NAXIS1, 10), H = parseInt(hdr.NAXIS2, 10), nc = naxis === 3 ? parseInt(hdr.NAXIS3, 10) : 1;
   if (nc !== 1 && nc !== 3) throw new Error("expected 1 or 3 planes, got " + nc);
   const bzero = parseFloat(hdr.BZERO || "0"), bscale = parseFloat(hdr.BSCALE || "1");
   const bytes = Math.abs(bitpix) / 8, N = W * H;
   if (buf.length < pos + N * nc * bytes) throw new Error("file truncated");
   const range = { 8: 255, 16: 65535, 32: 4294967295 }[bitpix];
   const intOffset = { 8: 0, 16: 32768, 32: 2147483648 }[bitpix];
   const chans = [];
   let fmax = 0;
   for (let c = 0; c < nc; ++c) {
      const a = new Float32Array(N), base = pos + c * N * bytes;
      for (let i = 0; i < N; ++i) {
         const o = base + i * bytes;
         let v;
         switch (bitpix) {
            case 8: v = buf.readUInt8(o); break;
            case 16: v = buf.readInt16BE(o); break;
            case 32: v = buf.readInt32BE(o); break;
            case -32: v = buf.readFloatBE(o); break;
            case -64: v = buf.readDoubleBE(o); break;
            default: throw new Error("BITPIX " + bitpix + " not supported");
         }
         v = v * bscale + bzero;
         // integer data mapped to [0, 1] as PixInsight does: unsigned data (8-bit, or BZERO = 2^(n-1))
         // over its full range, signed data shifted to it first
         if (range) v = (bitpix === 8 || (bzero === intOffset && bscale === 1) ? v : v + intOffset) / range;
         a[i] = v;
         if (v > fmax && isFinite(v)) fmax = v;
      }
      chans.push(a);
   }
   let note = "";
   if (!range && fmax > 1) {
      // floating-point stack stored in ADU: rescaled to [0, 1]
      const k = fmax <= 65535 ? 65535 : fmax;
      for (const a of chans) for (let i = 0; i < N; ++i) a[i] /= k;
      note = "floating-point data above 1, divided by " + k;
   }
   return { W, H, nc, chans, kw, note };
}

// ---------------------------------------------------------------- PNG writer (8-bit RGB)
const CRC = (() => {
   const t = new Uint32Array(256);
   for (let n = 0; n < 256; ++n) { let c = n; for (let k = 0; k < 8; ++k) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; }
   return t;
})();
function crc32(b) { let c = 0xffffffff; for (let i = 0; i < b.length; ++i) c = CRC[(c ^ b[i]) & 255] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; }
function chunk(type, data) {
   const len = Buffer.alloc(4); len.writeUInt32BE(data.length);
   const td = Buffer.concat([Buffer.from(type, "latin1"), data]);
   const crc = Buffer.alloc(4); crc.writeUInt32BE(crc32(td));
   return Buffer.concat([len, td, crc]);
}
// rgb: Uint8Array of w*h*3
function encodePng(rgb, w, h) {
   const raw = Buffer.alloc((w * 3 + 1) * h);
   for (let y = 0; y < h; ++y) { raw[y * (w * 3 + 1)] = 0; Buffer.from(rgb.buffer, rgb.byteOffset + y * w * 3, w * 3).copy(raw, y * (w * 3 + 1) + 1); }
   const ihdr = Buffer.alloc(13);
   ihdr.writeUInt32BE(w, 0); ihdr.writeUInt32BE(h, 4); ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
   return Buffer.concat([Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]), chunk("IHDR", ihdr),
                         chunk("IDAT", zlib.deflateSync(raw, { level: 9 })), chunk("IEND", Buffer.alloc(0))]);
}

// ---------------------------------------------------------------- helpers
function median(arr) {
   const a = Float64Array.from(arr.filter(v => isFinite(v))).sort();
   if (!a.length) return NaN;
   const m = a.length >> 1;
   return a.length % 2 ? a[m] : 0.5 * (a[m - 1] + a[m]);
}
// Linked screen transfer function from the luminance of the whole image (PixInsight AutoStretch formula).
function makeStf(chans, W, H) {
   const sample = [], step = Math.max(1, Math.floor(W * H / 400000));
   for (let i = 0; i < W * H; i += step) { let v = 0; for (const a of chans) v += a[i]; sample.push(v / chans.length); }
   const m = median(sample), mad = 1.4826 * median(sample.map(v => Math.abs(v - m)));
   const lo = Math.max(m - 2.8 * mad, 0), mm = (m - lo) / (1 - lo), t = 0.25;
   const mb = ((t - 1) * mm) / ((2 * t - 1) * mm - t);
   return v => { const x = Math.min(Math.max((v - lo) / (1 - lo), 0), 1); return ((mb - 1) * x) / ((2 * mb - 1) * x - mb); };
}

// Per-star figures in the ghost zone (22-40 px at the reference scale, background from 62-72 px).
// Real images have no ground truth, so these describe what the correction did; they are meant to
// be compared between versions, and judged with the contact sheet:
//   removed: light removed, % of the light above background in the zone
//   uneven:  spread of the removed light over 8 sectors, (max - min) / mean in %; 0 = round ghost
//   colour:  largest change of a channel share of the zone (points) between raw and corrected
function starMetrics(raw, ghost, nc, W, H, s, sc) {
   const rIn = 22 * sc, rOut = 40 * sc, rBg0 = 62 * sc, rBg1 = 72 * sc, R = Math.ceil(rBg1);
   const x0 = Math.round(s.x), y0 = Math.round(s.y), bg = [];
   for (let c = 0; c < nc; ++c) bg.push([]);
   const secG = new Float64Array(8), zr = new Float64Array(nc), zg = new Float64Array(nc);
   for (let y = y0 - R; y <= y0 + R; ++y) for (let x = x0 - R; x <= x0 + R; ++x) {
      if (x < 0 || y < 0 || x >= W || y >= H) continue;
      const r = Math.hypot(x - s.x, y - s.y), q = y * W + x;
      if (r >= rBg0 && r <= rBg1) for (let c = 0; c < nc; ++c) bg[c].push(raw[c][q]);
   }
   const bgc = bg.map(median);
   for (let y = y0 - R; y <= y0 + R; ++y) for (let x = x0 - R; x <= x0 + R; ++x) {
      if (x < 0 || y < 0 || x >= W || y >= H) continue;
      const r = Math.hypot(x - s.x, y - s.y);
      if (r < rIn || r > rOut) continue;
      let a = Math.atan2(y - s.y, x - s.x); if (a < 0) a += 2 * Math.PI;
      const k = Math.min(Math.floor(a / (Math.PI / 4)), 7), q = y * W + x;
      for (let c = 0; c < nc; ++c) {
         const v = raw[c][q] - bgc[c], g = ghost[c][q];
         zr[c] += v; zg[c] += g; secG[k] += g;
      }
   }
   const sum = a => a.reduce((p, v) => p + v, 0), totR = sum(zr), totG = sum(zg);
   const removed = totR > 0 ? 100 * totG / totR : NaN;
   const mg = totG / 8;
   const uneven = mg > 0 ? 100 * (Math.max(...secG) - Math.min(...secG)) / mg : NaN;
   let colour = NaN;
   if (nc === 3 && totR > 0 && totR - totG > 0)
      colour = 100 * Math.max(...[0, 1, 2].map(c => Math.abs((zr[c] - zg[c]) / (totR - totG) - zr[c] / totR)));
   return { removed, uneven, colour };
}

// Contact sheet: one row per pair of stars, each star as raw | corrected. The image is flipped
// vertically so that it appears as in PixInsight (FITS rows are stored bottom-up).
function contactSheet(raw, cor, nc, W, H, stars, half, stf) {
   const tile = 2 * half, gap = 4, pairW = 2 * tile + gap, cols = 2;
   const rows = Math.ceil(stars.length / cols), w = cols * pairW + (cols - 1) * 3 * gap, h = rows * (tile + 3 * gap);
   const out = new Uint8Array(w * h * 3);
   stars.forEach((s, i) => {
      const ox0 = (i % cols) * (pairW + 3 * gap), oy = Math.floor(i / cols) * (tile + 3 * gap);
      const cx = Math.min(Math.max(Math.round(s.x), half), W - half), cy = Math.min(Math.max(Math.round(s.y), half), H - half);
      [raw, cor].forEach((img, k) => {
         const ox = ox0 + k * (tile + gap);
         for (let j = 0; j < tile; ++j) for (let i2 = 0; i2 < tile; ++i2) {
            const x = cx - half + i2, y = cy + half - 1 - j, q = y * W + x, o = ((oy + j) * w + ox + i2) * 3;
            if (x < 0 || y < 0 || x >= W || y >= H) continue;
            for (let c = 0; c < 3; ++c) out[o + c] = Math.round(255 * stf(img[nc === 3 ? c : 0][q]));
         }
      });
   });
   return encodePng(out, w, h);
}

// ---------------------------------------------------------------- main
function listInputs(args) {
   const files = [];
   for (const a of args) {
      const st = fs.statSync(a);
      if (st.isDirectory()) for (const f of fs.readdirSync(a).sort()) { if (/\.(fits?|fts)$/i.test(f)) files.push(path.join(a, f)); }
      else files.push(a);
   }
   return files;
}

function main() {
   const argv = process.argv.slice(2), inputs = [];
   let outRoot = "bench-results", nStars = 8, all = false;
   for (let i = 0; i < argv.length; ++i) {
      if (argv[i] === "--out") outRoot = argv[++i];
      else if (argv[i] === "--stars") nStars = Math.max(1, parseInt(argv[++i], 10) || 8);
      else if (argv[i] === "--all") all = true;
      else inputs.push(argv[i]);
   }
   if (!inputs.length) { console.log("usage: node bench/run_bench.js <file or folder> [...] [--out dir] [--stars n] [--all]"); process.exit(1); }
   const GB = loadCore(), { version } = source();
   const outDir = path.join(outRoot, "v" + version);
   fs.mkdirSync(outDir, { recursive: true });
   const files = listInputs(inputs);
   const report = ["# GhostBuster " + version + " - benchmark on real stacks", "",
      "Automated mode, " + (all ? "all detected stars corrected" : "default selection (white stars)") + ".",
      "Figures in the ghost zone (22-40 px at the reference scale). Real stacks have no ground truth: these",
      "figures describe what the correction did and are meant for comparison between versions; the contact",
      "sheets are the reference.", "",
      "- **removed**: light removed, % of the light above the local background",
      "- **uneven**: spread of the removed light over 8 sectors, (max - min) / mean in %; 0 = round ghost, large = lopsided",
      "- **colour**: largest change of a channel share between raw and corrected (points)", ""];
   const results = { version, images: [] };
   for (const f of files) {
      const name = path.basename(f).replace(/\.(fits?|fts)$/i, "");
      const slug = name.replace(/[^A-Za-z0-9._-]+/g, "_");
      process.stdout.write(name + ": reading... ");
      let im;
      try { im = readFits(f); } catch (e) { console.log("skipped (" + e.message + ")"); report.push("## " + name, "", "Skipped: " + e.message, ""); continue; }
      const auto = GB.autoScale(im.kw, im.W, im.H), p = GB.defaults();
      p.scale = auto.scale;
      process.stdout.write("analysing... ");
      const dry = GB.process(im.chans, im.W, im.H, p, null, () => {}, { dryRun: true });
      const key = s => Math.round(s.x) + "," + Math.round(s.y);
      const pq = GB.scaled(p), chosen = {};
      for (const s of dry.stars) if (all || (s.status === "ok" && s.conf >= pq.confLo && s.area >= pq.minCoreArea)) chosen[key(s)] = s;
      process.stdout.write("correcting... ");
      const t0 = Date.now(), log = [];
      const res = GB.process(im.chans, im.W, im.H, p, s => !!chosen[key(s)], m => log.push(m));
      const dt = (Date.now() - t0) / 1000;
      const cor = im.chans.map((a, c) => { const o = new Float32Array(a.length); for (let i = 0; i < a.length; ++i) o[i] = Math.min(Math.max(a[i] - res.ghost[c][i], 0), 1); return o; });
      const done = Object.values(chosen).sort((a, b) => (b.rs || 0) - (a.rs || 0));
      const rows = done.map(s => {
         return Object.assign({ x: Math.round(s.x), y: Math.round(s.y), conf: s.conf },
                              starMetrics(im.chans, res.ghost, im.nc, im.W, im.H, s, pq.scale));
      });
      const shown = done.slice(0, nStars);
      const png = contactSheet(im.chans, cor, im.nc, im.W, im.H, shown, Math.round(90 * pq.scale), makeStf(im.chans, im.W, im.H));
      fs.writeFileSync(path.join(outDir, slug + ".png"), png);
      const f1 = v => isFinite(v) ? v.toFixed(1) : "-";
      report.push("## " + name, "",
         `${im.W} x ${im.H}, ${im.nc} channel(s); scale ${auto.scale.toFixed(2)} (${auto.source})${im.note ? "; " + im.note : ""}; ` +
         `${dry.stars.length} candidate stars, ${done.length} corrected in ${dt.toFixed(1)} s.`, "",
         `![${name}](${slug}.png)`, "",
         `Contact sheet: the ${shown.length} largest corrected stars, raw | corrected, two per row, in the order of the table.`, "",
         "| # | star (x, y) | confidence | removed % | uneven % | colour pts |", "|---|---|---|---|---|---|");
      rows.forEach((r, i) => report.push(`| ${i + 1} | (${r.x}, ${r.y}) | ${r.conf.toFixed(1)} | ${f1(r.removed)} | ${f1(r.uneven)} | ${f1(r.colour)} |`));
      report.push("");
      results.images.push({ name, width: im.W, height: im.H, channels: im.nc, scale: auto.scale, scaleSource: auto.source,
                            candidates: dry.stars.length, corrected: done.length, seconds: dt, stars: rows });
      console.log(`${done.length} stars, ${dt.toFixed(1)} s`);
   }
   fs.writeFileSync(path.join(outDir, "report.md"), report.join("\n") + "\n");
   fs.writeFileSync(path.join(outDir, "results.json"), JSON.stringify(results, null, 1) + "\n");
   console.log("Results in " + outDir);
}

if (require.main === module) main();
module.exports = { readFits, encodePng, crc32, starMetrics };
