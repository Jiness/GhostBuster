// Test suite for the GhostBuster processing core and dialog logic.
// Usage: node tests/run_tests.js        (Node.js >= 14, no dependencies)
"use strict";
const assert = require("assert");
const { loadCore } = require("./load.js");
const { makeField } = require("./synth.js");

const GB = loadCore();
let passed = 0, failed = 0;
function test(name, fn) {
   const t0 = Date.now();
   try { fn(); ++passed; console.log("  ok    " + name + " (" + (Date.now() - t0) + " ms)"); }
   catch (e) { ++failed; console.log("  FAIL  " + name + "\n        " + (e && e.message)); }
}
const quiet = () => {};
function checksum(chans) { let s = 0; for (const a of chans) for (let i = 0; i < a.length; i += 7) s += a[i] * (i % 13 + 1); return s; }

// Annulus error between the injected ghost and the estimated ghost, relative to the injected one.
function ghostError(f, est, star) {
   let truth = 0, err = 0;
   for (let y = Math.round(star.y) - 46; y <= Math.round(star.y) + 46; ++y)
      for (let x = Math.round(star.x) - 46; x <= Math.round(star.x) + 46; ++x) {
         const d = Math.hypot(x - star.x, y - star.y);
         if (d < 18 || d > 44 || x < 0 || y < 0 || x >= f.W || y >= f.H) continue;
         const q = y * f.W + x;
         for (let c = 0; c < f.nc; ++c) { truth += f.ghost[c][q]; err += Math.abs(f.ghost[c][q] - est[c][q]); }
      }
   return err / truth;
}

console.log("Numerical helpers");
test("median: odd and even lengths", () => {
   assert.strictEqual(GB.median([5, 1, 3], 3), 3);
   assert.strictEqual(GB.median([4, 1, 3, 2], 4), 2.5);
});
test("median: non-finite values are ignored", () => {
   assert.strictEqual(GB.median([NaN, 2, Infinity, 1, 3], 5), 2);
   assert.strictEqual(GB.median([NaN, NaN], 2), 0);
});
test("select: k-th smallest", () => {
   const a = new Float64Array([9, 3, 7, 1, 5]);
   assert.strictEqual(GB.select(a, 5, 0), 1);
   assert.strictEqual(GB.select(new Float64Array([9, 3, 7, 1, 5]), 5, 4), 9);
});
test("solve: 3x3 normal equations", () => {
   // A^T A x = A^T y for x = (1, -2, 3)
   const ata = [4, 1, 0, 1, 3, 1, 0, 1, 2], x = [1, -2, 3];
   const aty = [0, 1, 2].map(i => ata[3 * i] * x[0] + ata[3 * i + 1] * x[1] + ata[3 * i + 2] * x[2]);
   const s = GB.solve(ata, aty, 3);
   for (let i = 0; i < 3; ++i) assert.ok(Math.abs(s[i] - x[i]) < 1e-9);
});
test("lsqRobust: line fit ignores outliers", () => {
   const n = 200, c0 = new Float64Array(n), c1 = new Float64Array(n), y = new Float64Array(n), ok = new Uint8Array(n);
   for (let i = 0; i < n; ++i) { c0[i] = 1; c1[i] = i / n; y[i] = 2 + 3 * c1[i] + (i % 10 == 0 ? 50 : 0); ok[i] = 1; }
   const r = GB.lsqRobust([c0, c1], y, n, ok, 5, 2.5, 4);
   assert.ok(Math.abs(r.coef[0] - 2) < 1e-3 && Math.abs(r.coef[1] - 3) < 1e-3, r.coef[0] + " " + r.coef[1]);
});

console.log("Settings");
test("validate: out-of-range and non-finite values are clamped", () => {
   const p = GB.validate({ scale: NaN, annIn: -5, annOut: 3, rOut: 1, confLo: 2, confHi: 2, nHarm: 99, iters: 0.2 });
   assert.strictEqual(p.scale, 1);
   assert.ok(p.annIn >= 2 && p.annOut >= p.annIn + 4 && p.rOut >= p.annOut + 10);
   assert.ok(p.confHi > p.confLo);
   assert.ok(p.nHarm <= 6 && p.iters >= 1);
});
test("validate: defaults are unchanged", () => {
   const d = GB.defaults(), v = GB.validate(d);
   for (const k in d) assert.deepStrictEqual(v[k], d[k], k);
});
test("autoScale: astrometric solution of a native Draco stack", () => {
   const kw = [["CD1_1", "-0.000274181350922"], ["CD1_2", "-0.000296601211037"], ["CD2_1", "0.000296579170808"], ["CD2_2", "-0.000274363011826"]]
      .map(([name, value]) => ({ name, value }));
   const a = GB.autoScale(kw, 4096, 3072);
   assert.strictEqual(a.scale, 1); assert.ok(a.fromMetadata);
});
test("autoScale: 2x drizzle from the astrometric solution", () => {
   const kw = [["CD1_1", "-0.000137"], ["CD1_2", "-0.000148"], ["CD2_1", "0.000148"], ["CD2_2", "-0.000137"]].map(([name, value]) => ({ name, value }));
   assert.strictEqual(GB.autoScale(kw, 8150, 6100).scale, 2);
});
test("autoScale: pixel size and focal length", () => {
   const a = GB.autoScale([{ name: "XPIXSZ", value: "2.394" }, { name: "FOCALLEN", value: "340." }], 4096, 3072);
   assert.strictEqual(a.scale, 1); assert.ok(a.fromMetadata);
});
test("autoScale: no metadata falls back to the image size and says so", () => {
   const a = GB.autoScale([], 8150, 6100);
   assert.strictEqual(a.scale, 2); assert.strictEqual(a.fromMetadata, false);
   assert.strictEqual(GB.autoScale([], 4000, 3000).scale, 1);
});
test("autoScale: a cropped image keeps the scale of its pixel size", () => {
   const kw = [{ name: "XPIXSZ", value: "2.394" }, { name: "FOCALLEN", value: "340." }];
   assert.strictEqual(GB.autoScale(kw, 1200, 900).scale, 1);
   assert.strictEqual(GB.autoScale([], 1200, 900).scale, 1);
});
test("autoScale: drizzled image whose header kept the original pixel size", () => {
   const kw = [{ name: "XPIXSZ", value: "2.394" }, { name: "FOCALLEN", value: "340." }];
   assert.strictEqual(GB.autoScale(kw, 8150, 6100).scale, 2);
});

console.log("Ghost removal on synthetic fields");
for (const seed of [1, 2]) {
   const f = makeField({ seed });
   const before = checksum(f.chans);
   let res;
   test("seed " + seed + ": all ghosted stars detected and corrected", () => {
      res = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet);
      assert.strictEqual(res.done, f.stars.filter(s => s.ghost).length);
   });
   test("seed " + seed + ": ghost estimate error below 60% overall and 80% per star", () => {
      let tot = 0;
      const ghosted = f.stars.filter(s => s.ghost);
      for (const s of ghosted) { const e = ghostError(f, res.ghost, s); tot += e; assert.ok(e < 0.8, "star " + s.x.toFixed(0) + ": " + e.toFixed(2)); }
      assert.ok(tot / ghosted.length < 0.6, "overall " + (tot / ghosted.length).toFixed(2));
   });
   test("seed " + seed + ": ghost colour learned (red > green, blue)", () => {
      let r = 0, g = 0, b = 0;
      for (let q = 0; q < f.W * f.H; ++q) { r += res.ghost[0][q]; g += res.ghost[1][q]; b += res.ghost[2][q]; }
      assert.ok(r > 1.2 * g && r > 1.2 * b, [r, g, b].map(v => v.toFixed(1)).join("/"));
   });
   test("seed " + seed + ": input arrays are not modified", () => {
      assert.strictEqual(checksum(f.chans), before);
   });
}

console.log("Colour and coloured nebula");
// Colour shares (R, G, B) of an image sum over the ghost disk of a star.
function shares(chans, f, star) {
   const v = [0, 0, 0];
   for (let y = Math.round(star.y) - 44; y <= Math.round(star.y) + 44; ++y)
      for (let x = Math.round(star.x) - 44; x <= Math.round(star.x) + 44; ++x) {
         const d = Math.hypot(x - star.x, y - star.y);
         if (d < 14 || d > 40 || x < 0 || y < 0 || x >= f.W || y >= f.H) continue;
         for (let c = 0; c < 3; ++c) v[c] += chans[c][y * f.W + x];
      }
   const s = v[0] + v[1] + v[2];
   return v.map(a => a / s);
}
{
   const f = makeField({ seed: 4, bright: 6 });
   const res = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet);
   test("ghost colour of each star within 0.06 of the injected one", () => {
      for (const s of f.stars.filter(t => t.ghost)) {
         const est = shares(res.ghost, f, s), tru = shares(f.ghost, f, s);
         for (let c = 0; c < 3; ++c) assert.ok(Math.abs(est[c] - tru[c]) < 0.06, "star " + s.x.toFixed(0) + ": " + est.map(v => v.toFixed(3)) + " vs " + tru.map(v => v.toFixed(3)));
      }
   });
}
for (const amp of [0.01, 0.03]) {
   // red emission patch next to one ghosted star, from about the ghost level to 3 times brighter
   const f = makeField({ seed: 5, bright: 6 });
   const s = f.stars.find(t => t.ghost);
   for (let y = 0; y < f.H; ++y) for (let x = 0; x < f.W; ++x) {
      const d = Math.hypot(x - s.x - 25, y - s.y + 15), v = amp * Math.exp(-d * d / (2 * 60 * 60));
      const q = y * f.W + x;
      f.chans[0][q] = Math.min(1, f.chans[0][q] + v); f.clean[0][q] += v;
   }
   const res = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet);
   test("red nebula " + amp + ": no sector of the corrected disk over-subtracted by 60% of the ghost", () => {
      // residual = corrected - clean, by 45 degree sector, compared with the injected ghost
      for (let sec = 0; sec < 8; ++sec) {
         const rsd = [0, 0, 0], gh = [0, 0, 0];
         for (let y = Math.round(s.y) - 44; y <= Math.round(s.y) + 44; ++y)
            for (let x = Math.round(s.x) - 44; x <= Math.round(s.x) + 44; ++x) {
               const d = Math.hypot(x - s.x, y - s.y);
               let a = Math.atan2(y - s.y, x - s.x); if (a < 0) a += 2 * Math.PI;
               if (d < 18 || d > 38 || Math.floor(a / (Math.PI / 4)) !== sec) continue;
               const q = y * f.W + x;
               for (let c = 0; c < 3; ++c) { rsd[c] += f.chans[c][q] - res.ghost[c][q] - f.clean[c][q]; gh[c] += f.ghost[c][q]; }
            }
         for (let c = 0; c < 3; ++c) assert.ok(rsd[c] > -0.6 * gh[c], "sector " + sec + " channel " + c + ": residual " + (rsd[c] / gh[c]).toFixed(2) + " of the ghost");
      }
   });
}

console.log("Benchmark tool");
{
   const fs = require("fs"), os = require("os"), path = require("path"), zlib = require("zlib");
   const bench = require("../bench/run_bench.js");
   // minimal FITS writer for the tests: 3 planes, given BITPIX
   function fitsFile(W, H, bitpix, values, extra) {
      const cards = ["SIMPLE  =                    T", `BITPIX  = ${String(bitpix).padStart(20)}`, "NAXIS   =                    3",
                     `NAXIS1  = ${String(W).padStart(20)}`, `NAXIS2  = ${String(H).padStart(20)}`, "NAXIS3  =                    3"].concat(extra || [], ["END"]);
      let hdr = cards.map(c => c.padEnd(80)).join("");
      hdr = hdr.padEnd(Math.ceil(hdr.length / 2880) * 2880);
      const bytes = Math.abs(bitpix) / 8, data = Buffer.alloc(Math.ceil(W * H * 3 * bytes / 2880) * 2880);
      values.forEach((v, i) => bitpix === 16 ? data.writeInt16BE(v, i * 2) : data.writeFloatBE(v, i * 4));
      const f = path.join(os.tmpdir(), "gb_test_" + process.pid + "_" + bitpix + ".fits");
      fs.writeFileSync(f, Buffer.concat([Buffer.from(hdr, "latin1"), data]));
      return f;
   }
   test("FITS reader: unsigned 16-bit (BZERO 32768) mapped to [0, 1], keywords read", () => {
      const raw = [-32768, 0, 32767, -32768, 0, 32767, 0, 0, 0, 0, 0, 0];
      const f = fitsFile(2, 2, 16, raw, ["BZERO   =                32768", "BSCALE  =                    1", "XPIXSZ  =                2.394"]);
      const im = bench.readFits(f); fs.unlinkSync(f);
      assert.strictEqual(im.W, 2); assert.strictEqual(im.nc, 3);
      assert.strictEqual(im.chans[0][0], 0); assert.strictEqual(im.chans[0][2], 1);
      assert.ok(Math.abs(im.chans[0][1] - 32768 / 65535) < 1e-6);
      assert.strictEqual(im.kw.find(k => k.name === "XPIXSZ").value, "2.394");
   });
   test("FITS reader: 32-bit float above 1 rescaled from 16-bit ADU", () => {
      const f = fitsFile(2, 1, -32, [0, 65535, 100, 200, 300, 400]);
      const im = bench.readFits(f); fs.unlinkSync(f);
      assert.strictEqual(im.chans[0][1], 1); assert.ok(im.note.length > 0);
   });
   test("PNG writer: valid signature, chunk CRCs and pixel data", () => {
      const rgb = new Uint8Array([255, 0, 0, 0, 255, 0, 0, 0, 255, 10, 20, 30]);
      const png = bench.encodePng(rgb, 2, 2);
      assert.deepStrictEqual(Array.from(png.slice(0, 8)), [137, 80, 78, 71, 13, 10, 26, 10]);
      let pos = 8, idat = [];
      while (pos < png.length) {
         const len = png.readUInt32BE(pos), type = png.toString("latin1", pos + 4, pos + 8);
         assert.strictEqual(png.readUInt32BE(pos + 8 + len), bench.crc32(png.slice(pos + 4, pos + 8 + len)), type + " CRC");
         if (type === "IDAT") idat.push(png.slice(pos + 8, pos + 8 + len));
         pos += 12 + len;
      }
      const raw = zlib.inflateSync(Buffer.concat(idat));
      assert.deepStrictEqual(Array.from(raw), [0, 255, 0, 0, 0, 255, 0, 0, 0, 0, 255, 10, 20, 30]);
   });
}

console.log("Analysis and manual stars");
{
   const f = makeField({ seed: 3, bright: 5 });
   test("dry run classifies every candidate and allocates no output", () => {
      const r = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet, { dryRun: true });
      assert.ok(r.stars.length >= 5);
      for (const s of r.stars) assert.ok(["ok", "gate", "noghost"].indexOf(s.status) >= 0);
      assert.strictEqual(r.ghost.length, 0);
   });
   test("invalid manual stars are ignored (NaN, out of bounds, null)", () => {
      const extra = [{ x: NaN, y: 10 }, { x: -50, y: 10 }, { x: 10, y: f.H + 100 }, null, { x: f.W, y: 0 }];
      const r = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet, { dryRun: true, extra });
      assert.strictEqual(r.stars.filter(s => s.manual).length, 0);
   });
   test("valid manual star is snapped to the nearest bright pixel and caller data is untouched", () => {
      const faint = f.stars.find(s => !s.ghost && s.peak > 0.15 && s.x > 50 && s.y > 50 && s.x < f.W - 50 && s.y < f.H - 50);
      const ex = { x: faint.x + 4, y: faint.y - 3 };
      const r = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet, { dryRun: true, extra: [ex] });
      const m = r.stars.filter(s => s.manual);
      assert.strictEqual(m.length, 1);
      assert.ok(Math.hypot(m[0].x - faint.x, m[0].y - faint.y) <= 2);
      assert.deepStrictEqual(Object.keys(ex).sort(), ["x", "y"]);
   });
   test("selectFn excludes stars", () => {
      const r = GB.process(f.chans, f.W, f.H, GB.defaults(), () => false, quiet);
      assert.strictEqual(r.done, 0);
      for (const g of r.ghost) for (let q = 0; q < g.length; q += 101) assert.strictEqual(g[q], 0);
   });
}

console.log("Robustness");
const fin = (r) => r.ghost.every(g => { for (let i = 0; i < g.length; ++i) if (!isFinite(g[i])) return false; return true; });
test("mono image", () => {
   const f = makeField({ seed: 4, nc: 1, bright: 4, W: 800, H: 600 });
   const r = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet);
   assert.ok(r.done >= 1 && fin(r));
});
test("image containing NaN and Infinity", () => {
   const f = makeField({ seed: 5, bright: 3, W: 700, H: 500 });
   for (const a of f.chans) for (let i = 0; i < a.length; i += 977) a[i] = (i % 2) ? NaN : Infinity;
   const r = GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet);
   assert.ok(fin(r));
});
test("tiny image (40 x 30)", () => {
   const f = makeField({ seed: 6, bright: 0, faint: 3, W: 40, H: 30 });
   assert.ok(fin(GB.process(f.chans, f.W, f.H, GB.defaults(), null, quiet)));
});
test("fully saturated image is not taken for a star", () => {
   const W = 600, H = 400, ch = [0, 1, 2].map(() => new Float32Array(W * H).fill(1));
   const r = GB.process(ch, W, H, GB.defaults(), null, quiet);
   assert.strictEqual(r.total, 0);
});
test("black image", () => {
   const W = 600, H = 400, ch = [0, 1, 2].map(() => new Float32Array(W * H));
   assert.strictEqual(GB.process(ch, W, H, GB.defaults(), null, quiet).total, 0);
});
test("absurd settings are clamped instead of failing", () => {
   const f = makeField({ seed: 7, bright: 2, W: 600, H: 450 });
   const r = GB.process(f.chans, f.W, f.H, { scale: 0, annIn: 1e9, rStep: 0, thStepDeg: -1, confLo: 5, confHi: 1 }, null, quiet);
   assert.ok(fin(r));
});

console.log("\n" + passed + " passed, " + failed + " failed");
process.exitCode = failed ? 1 : 0;
