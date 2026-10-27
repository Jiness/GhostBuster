// Deterministic synthetic star fields with known reflection ghosts.
"use strict";

function rng(seed) {
   let s = seed >>> 0;
   const next = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
   return {
      next,
      gauss() { const u = Math.max(next(), 1e-12), v = next(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
   };
}

// Returns { W, H, chans, clean, ghost, stars }:
//   chans: observed image (star + ghost + background + noise), values in [0, 1]
//   clean: same image without the ghosts, ghost: the injected ghosts only
function makeField(opts) {
   const o = Object.assign({ W: 1024, H: 768, nc: 3, seed: 1, bright: 8, faint: 150, ghostRatio: 0.0003,
                             ghostColour: [1.5, 0.8, 0.8], background: 0.06, noise: 0.0015 }, opts || {});
   const { W, H, nc } = o, N = W * H, r = rng(o.seed);
   const clean = [], ghost = [];
   for (let c = 0; c < nc; ++c) { clean.push(new Float32Array(N).fill(o.background)); ghost.push(new Float32Array(N)); }
   const stars = [];
   const addStar = (x, y, peak, colour, withGhost) => {
      const R = withGhost ? 110 : 12;
      const gdx = withGhost ? (r.next() - 0.5) * 6 : 0, gdy = withGhost ? (r.next() - 0.5) * 6 : 0;
      for (let yy = Math.max(0, Math.floor(y - R)); yy <= Math.min(H - 1, Math.ceil(y + R)); ++yy)
         for (let xx = Math.max(0, Math.floor(x - R)); xx <= Math.min(W - 1, Math.ceil(x + R)); ++xx) {
            const d = Math.hypot(xx - x, yy - y);
            // Moffat-like core plus a power-law wing
            const core = peak * Math.pow(1 + (d / 2.2) * (d / 2.2), -2.5) + peak * 0.0006 * Math.pow(Math.max(d, 3) / 20, -2.6);
            let gh = 0;
            if (withGhost) {
               const dg = Math.hypot(xx - x - gdx, yy - y - gdy);
               gh = peak * o.ghostRatio * 0.5 * (1 - Math.tanh((dg - 33) / 1.2)) * (0.7 + 0.3 * dg / 33);
            }
            const q = yy * W + xx;
            for (let c = 0; c < nc; ++c) {
               const col = nc == 1 ? 1 : colour[c];
               clean[c][q] += core * col;
               ghost[c][q] += gh * col * (nc == 1 ? 1 : o.ghostColour[c]) / 3;
            }
         }
      stars.push({ x, y, peak, ghost: withGhost, gdx, gdy });
   };
   // bright saturated stars with ghosts, kept apart from each other and from the edges
   const placed = [], margin = Math.min(130, W / 4, H / 4), minSep = Math.min(220, Math.max(W, H) / 3);
   for (let tries = 0; placed.length < o.bright && tries < 10000; ++tries) {
      const x = margin + r.next() * (W - 2 * margin), y = margin + r.next() * (H - 2 * margin);
      if (placed.every(p => Math.hypot(p[0] - x, p[1] - y) > minSep)) placed.push([x, y]);
   }
   for (const [x, y] of placed) addStar(x, y, 30 + r.next() * 30, [0.8 + 0.2 * r.next(), 1, 0.9 + 0.2 * r.next()], true);
   for (let i = 0; i < o.faint; ++i) addStar(r.next() * W, r.next() * H, 0.02 + r.next() * 0.3, [1, 1, 1], false);
   const chans = [];
   for (let c = 0; c < nc; ++c) {
      const a = new Float32Array(N);
      for (let q = 0; q < N; ++q) {
         a[q] = Math.min(1, Math.max(0, clean[c][q] + ghost[c][q] + o.noise * r.gauss()));
         clean[c][q] = Math.min(1, clean[c][q]);
      }
      chans.push(a);
   }
   return { W, H, nc, chans, clean, ghost, stars };
}

module.exports = { makeField, rng };
