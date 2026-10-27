// ****************************************************************************
// GhostBuster.js - PixInsight script
// https://github.com/Jiness/GhostBuster
//
// Copyright (C) 2026 Jiness
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version. This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details: <https://www.gnu.org/licenses/>.
//
// Removes the reflection ghost rings / disks that appear around bright,
// saturated stars on Dwarflab Draco duo-band stacks (and on other instruments
// with internal reflections).
//
// METHOD
//   1. Candidate stars: saturated cores.
//   2. Learning, on stars in dark sky: radial shape of the ghost, typical
//      ghost/star intensity ratio and ghost/star colour ratio. A built-in
//      Draco profile is used when the image has too few such stars.
//   3. For each star: model of the smooth stellar halo and of the local
//      background (quadratic surface); the excess in an annulus is the ghost
//      signal, with neighbouring stars masked. It is fitted with the learned
//      shape (radius adjusted per star, limited lopsidedness), then refined by
//      a bounded free correction. A fit above the physical ceiling (ratio
//      cap) has absorbed nebula and falls back to the typical ghost.
//   4. Ghost colour: split chosen so that the corrected disk keeps the colour
//      of the star light around it, within bounds around the predicted colour
//      (star colour x learned colour ratio).
//   5. Confidence = ghost amplitude / residual structure. The freedom of the
//      shape scales with it: weak ghosts in structured nebula stay regular
//      rings.
//
// USAGE
//   * Linear data only, after background extraction, before deconvolution
//     and before star removal.
//   * Automated mode reads the pixel scale from the header (astrometric
//     solution, else pixel size and focal length, else image size) and
//     derives every setting. Manual mode exposes the settings.
//   * The dialog shows the image with a circle around each candidate star:
//     white = corrected, orange = uncertain (corrected only if selected),
//     grey = excluded, green = added manually.
//   * Outputs: <id>_deghosted and, optionally, <id>_ghost (removed signal).
//     The source image is not modified.
//
// Version 0.1.2
// ****************************************************************************

#feature-id    Utilities > GhostBuster
#feature-info  Removes reflection ghost rings around bright saturated stars \
               on linear images. Run before deconvolution and star removal.

#include <pjsr/NumericControl.jsh>
#include <pjsr/StdButton.jsh>
#include <pjsr/StdIcon.jsh>
#include <pjsr/UndoFlag.jsh>
#include <pjsr/ColorSpace.jsh>
#include <pjsr/SampleType.jsh>

#define TITLE   "GhostBuster"
#define VERSION "0.1.2"
// PixInsight version as a string (the predefined macros expand to bare numbers)
#ifdef __PI_MAJOR__
#define GB_PI_VERSION (__PI_MAJOR__ + "." + __PI_MINOR__ + "." + __PI_RELEASE__)
#else
#define GB_PI_VERSION "unknown"
#endif

// ===================== Core (no PixInsight dependencies) =====================
var GB = {};

GB.defaults = function () {
   return {
      satThreshold : 0.90,  // fraction of the brightest luminance value considered saturated
      minArea      : 12,    // minimum saturated area (px) for a star to be treated
      maxArea      : 5000,  // larger saturated regions are not stars (clipped nebula core, artefacts)
      scale        : 1.0,   // multiplies all radii (2.0 for a 2x drizzled image)
      rOut         : 90,    // halo fit radius (px)
      annIn        : 16,    // ghost annulus inner radius (px)
      annOut       : 46,    // ghost annulus outer radius (px)
      gateK        : 8.0,   // skip a star if halo-fit residual > gateK x image noise
      nHarm        : 3,     // angular harmonics for the ghost estimate in sparse fields
      maxStars     : 200,   // candidates kept (largest saturated cores first), bounds the analysis time
      maxLop       : 0.6,   // max lopsidedness of the ghost (fraction of its mean level)
      ratioCap     : 3.0,   // ghost may not exceed ratioCap x the typical ghost/star ratio
      freeBound    : 1.0,   // free correction allowed around the fitted shape (fraction of it)
      colLo        : 0.4,   // local ghost colour may deviate from the predicted one by these factors
      colHi        : 2.5,
      bgSpread     : 5.0,   // background variation around the star (x ghost level) above which the
                            // colour is fitted on the modelled residual instead of the star light
      confLo       : 1.5,   // below: regular round ring (weak ghost in nebula)
      confHi       : 3.0,   // above: full shape freedom (strong clear ghost)
      haloPows     : [2, 3, 4],
      iters        : 5,
      clipHi       : 2.5,
      clipLo       : 4.0,
      rStep        : 0.5,
      thStepDeg    : 2.0,
      edgeTaper    : 3.0,
      coreFactor   : 1.8,
      coreMin      : 10,
      satDilate    : 3
   };
};

// k-th smallest (quickselect, in place on b[0..n-1])
GB.select = function (b, n, k) {
   var lo = 0, hi = n - 1;
   while (hi > lo) {
      var piv = b[(lo + hi) >> 1], i = lo, j = hi;
      while (i <= j) {
         while (b[i] < piv) ++i;
         while (b[j] > piv) --j;
         if (i <= j) { var t = b[i]; b[i] = b[j]; b[j] = t; ++i; --j; }
      }
      if (k <= j) hi = j; else if (k >= i) lo = i; else return b[k];
   }
   return b[k];
};
GB.median = function (a, n) {
   // NaN / infinite values are ignored
   var b = new Float64Array(n), m = 0;
   for (var i = 0; i < n; ++i) { var v = a[i]; if (isFinite(v)) b[m++] = v; }
   if (m <= 0) return 0;
   if (m & 1) return GB.select(b, m, (m - 1) >> 1);
   var hi = GB.select(b, m, m >> 1), lo = -Infinity;
   for (var i = 0; i < (m >> 1); ++i) if (b[i] > lo) lo = b[i];
   return 0.5 * (lo + hi);
};

// Solve (A^T A) x = A^T y with partial pivoting. ata: m*m, aty: m
GB.solve = function (ata, aty, m) {
   var M = new Float64Array(m * (m + 1));
   var maxd = 0;
   for (var i = 0; i < m; ++i) if (ata[i * m + i] > maxd) maxd = ata[i * m + i];
   for (var i = 0; i < m; ++i) {
      for (var j = 0; j < m; ++j) M[i * (m + 1) + j] = ata[i * m + j];
      M[i * (m + 1) + i] += 1e-12 * maxd + 1e-300;
      M[i * (m + 1) + m] = aty[i];
   }
   var w = m + 1;
   for (var c = 0; c < m; ++c) {
      var p = c, best = Math.abs(M[c * w + c]);
      for (var r = c + 1; r < m; ++r) { var v = Math.abs(M[r * w + c]); if (v > best) { best = v; p = r; } }
      if (p != c) for (var k = 0; k < w; ++k) { var t = M[c * w + k]; M[c * w + k] = M[p * w + k]; M[p * w + k] = t; }
      var d = M[c * w + c];
      if (d == 0) continue;
      for (var r = c + 1; r < m; ++r) {
         var f = M[r * w + c] / d;
         if (f != 0) for (var k = c; k < w; ++k) M[r * w + k] -= f * M[c * w + k];
      }
   }
   var x = new Float64Array(m);
   for (var i = m - 1; i >= 0; --i) {
      var s = M[i * w + m];
      for (var k = i + 1; k < m; ++k) s -= M[i * w + k] * x[k];
      x[i] = (M[i * w + i] != 0) ? s / M[i * w + i] : 0;
   }
   return x;
};

// Robust least squares. cols: array of m Float64Array(n); y: Float64Array(n); ok: Uint8Array(n) initial usable rows
GB.lsqRobust = function (cols, y, n, ok, iters, clipHi, clipLo) {
   var m = cols.length, w = new Uint8Array(n), res = new Float64Array(n), tmp = new Float64Array(n);
   for (var i = 0; i < n; ++i) w[i] = ok[i];
   var coef = new Float64Array(m), sigma = 0;
   for (var it = 0; it < iters; ++it) {
      var ata = new Float64Array(m * m), aty = new Float64Array(m), cnt = 0;
      for (var i = 0; i < n; ++i) {
         if (!w[i]) continue;
         ++cnt;
         for (var a = 0; a < m; ++a) {
            var va = cols[a][i];
            aty[a] += va * y[i];
            for (var b = a; b < m; ++b) ata[a * m + b] += va * cols[b][i];
         }
      }
      if (cnt < m) break;
      for (var a = 0; a < m; ++a) for (var b = 0; b < a; ++b) ata[a * m + b] = ata[b * m + a];
      coef = GB.solve(ata, aty, m);
      var k = 0;
      for (var i = 0; i < n; ++i) {
         var s = 0;
         for (var a = 0; a < m; ++a) s += coef[a] * cols[a][i];
         res[i] = y[i] - s;
         if (w[i]) tmp[k++] = res[i];
      }
      var med = GB.median(tmp, k);
      for (var i = 0; i < k; ++i) tmp[i] = Math.abs(tmp[i] - med);
      sigma = 1.4826 * GB.median(tmp, k) + 1e-12;
      // clip around the median residual; keep the previous selection if everything would be rejected
      var kept = 0;
      for (var i = 0; i < n; ++i) {
         var d = res[i] - med;
         tmp[i] = ok[i] && d < clipHi * sigma && d > -clipLo * sigma ? 1 : 0;
         kept += tmp[i];
      }
      if (kept < m) break;
      for (var i = 0; i < n; ++i) w[i] = tmp[i];
   }
   return { coef: coef, sigma: sigma };
};

// chans: array of Float32Array (width*height), values normalized [0,1]
GB.detect = function (chans, W, H, p) {
   var nc = chans.length, N = W * H, L = new Float32Array(N), Lmax = 0;
   for (var i = 0; i < N; ++i) {
      var s = 0;
      for (var c = 0; c < nc; ++c) s += chans[c][i];
      L[i] = s / nc;
      if (L[i] > Lmax) Lmax = L[i];
   }
   var thr = p.satThreshold * Lmax, m = new Uint8Array(N);
   for (var i = 0; i < N; ++i) if (L[i] > thr) m[i] = 1;
   var lab = new Int32Array(N), stars = [], stack = new Int32Array(N), nl = 0;
   for (var i = 0; i < N; ++i) {
      if (!m[i] || lab[i]) continue;
      ++nl; var sp = 0; stack[sp++] = i; lab[i] = nl;
      var area = 0, sx = 0, sy = 0;
      while (sp > 0) {
         var q = stack[--sp], qx = q % W, qy = (q - qx) / W;
         ++area; sx += qx; sy += qy;
         var nb = [q - 1, q + 1, q - W, q + W], ok = [qx > 0, qx < W - 1, qy > 0, qy < H - 1];
         for (var k = 0; k < 4; ++k) if (ok[k] && m[nb[k]] && !lab[nb[k]]) { lab[nb[k]] = nl; stack[sp++] = nb[k]; }
      }
      if (area >= p.minArea && area <= p.maxArea) stars.push({ x: sx / area, y: sy / area, rs: Math.sqrt(area / Math.PI), area: area });
   }
   stars.sort(function (a, b) { return b.area - a.area; });
   // saturation mask dilated (diamond)
   var sat = new Uint8Array(N), D = p.satDilate;
   for (var i = 0; i < N; ++i) {
      if (!m[i]) continue;
      var x = i % W, y = (i - x) / W;
      for (var dy = -D; dy <= D; ++dy) for (var dx = -D; dx <= D; ++dx) {
         if (Math.abs(dx) + Math.abs(dy) > D) continue;
         var xx = x + dx, yy = y + dy;
         if (xx >= 0 && yy >= 0 && xx < W && yy < H) sat[yy * W + xx] = 1;
      }
   }
   return { stars: stars, sat: sat };
};

GB.noise = function (ch, W, H) {
   var a = [];
   for (var y = 0; y < H; y += 4) for (var x = 0; x < W; x += 4) a.push(ch[y * W + x]);
   var n = a.length, med = GB.median(a, n);
   for (var i = 0; i < n; ++i) a[i] = Math.abs(a[i] - med);
   return 1.4826 * GB.median(a, n);
};

// Returns a copy of p with every numeric setting finite and inside a safe range
// (array sizes, loop counts, interpolation steps and divisors depend on them).
GB.validate = function (p) {
   var d = GB.defaults(), q = {};
   for (var k in d) q[k] = d[k];
   for (var k in p) q[k] = p[k];
   function num(key, lo, hi) {
      var v = Number(q[key]);
      if (!isFinite(v)) v = d[key];
      q[key] = Math.min(Math.max(v, lo), hi);
   }
   num("scale", 0.25, 8);         num("satThreshold", 0.5, 1);
   num("minArea", 1, 1e6);        num("maxArea", 1, 1e7);
   num("annIn", 2, 500);          num("annOut", 6, 1000);
   num("rOut", 16, 2000);         num("gateK", 0.1, 1e3);
   num("nHarm", 1, 6);            num("maxLop", 0, 1);
   num("maxStars", 1, 5000);
   num("ratioCap", 0.1, 100);     num("freeBound", 0, 10);
   num("colLo", 0.01, 1);         num("colHi", 1, 100);
   num("bgSpread", 0, 1000);
   num("confLo", 0, 100);         num("confHi", 0, 1000);
   num("iters", 1, 20);           num("clipHi", 0.5, 20);
   num("clipLo", 0.5, 20);        num("rStep", 0.1, 5);
   num("thStepDeg", 0.5, 15);     num("edgeTaper", 0.5, 50);
   num("coreFactor", 1, 10);      num("coreMin", 1, 200);
   num("satDilate", 0, 20);
   if (q.maxArea < q.minArea) q.maxArea = q.minArea;
   if (q.annOut < q.annIn + 4) q.annOut = q.annIn + 4;
   if (q.rOut < q.annOut + 10) q.rOut = q.annOut + 10;
   if (q.confHi < q.confLo + 0.1) q.confHi = q.confLo + 0.1;
   q.iters = Math.round(q.iters);
   if (!(q.haloPows instanceof Array) || q.haloPows.length == 0) q.haloPows = d.haloPows;
   return q;
};

GB.scaled = function (p0) {
   var p = GB.validate(p0);
   var q = {}; for (var k in p) q[k] = p[k];
   var s = p.scale;
   q.rOut = p.rOut * s; q.annIn = p.annIn * s; q.annOut = p.annOut * s; q.coreMin = p.coreMin * s;
   q.edgeTaper = p.edgeTaper * s; q.rStep = p.rStep * s; q.minArea = Math.max(1, Math.round(p.minArea * s * s));
   q.satDilate = Math.max(1, Math.round(p.satDilate * s));
   q.maxArea = Math.round(p.maxArea * s * s);
   q.nHarm = Math.max(1, Math.round(p.nHarm));
   return q;
};

// Per-radius robust Fourier fit of the polar map. pol: Float64Array(nth*nr) row-major [t*nr+j], NaN = missing
GB.fourierRing = function (pol, nth, nr, thStep, nh) {
   var ncMax = 1 + 2 * nh, F = [];
   for (var c = 0; c < ncMax; ++c) F.push(new Float64Array(nth));
   for (var t = 0; t < nth; ++t) {
      var th = t * thStep; F[0][t] = 1;
      for (var m = 1; m <= nh; ++m) { F[2 * m - 1][t] = Math.cos(m * th); F[2 * m][t] = Math.sin(m * th); }
   }
   var out = new Float64Array(nth * nr);
   for (var j = 0; j < nr; ++j) {
      var j0 = Math.max(j - 1, 0), j1 = Math.min(j + 1, nr - 1), nn = nth * (j1 - j0 + 1);
      // angular coverage: fraction of angles with at least one valid sample
      var covered = 0;
      for (var t = 0; t < nth; ++t) {
         for (var jj = j0; jj <= j1; ++jj) if (!isNaN(pol[t * nr + jj])) { ++covered; break; }
      }
      var cov = covered / nth;
      if (cov < 0.15) continue;                       // too little data: leave 0
      var nhEff = cov > 0.6 ? nh : (cov > 0.35 ? Math.min(nh, 1) : 0);
      var nc = 1 + 2 * nhEff;
      var cols = [], y = new Float64Array(nn), ok = new Uint8Array(nn), vals = [];
      for (var c = 0; c < nc; ++c) cols.push(new Float64Array(nn));
      var k = 0;
      for (var t = 0; t < nth; ++t) for (var jj = j0; jj <= j1; ++jj, ++k) {
         var v = pol[t * nr + jj];
         y[k] = isNaN(v) ? 0 : v; ok[k] = isNaN(v) ? 0 : 1;
         if (!isNaN(v)) vals.push(v);
         for (var c = 0; c < nc; ++c) cols[c][k] = F[c][t];
      }
      if (vals.length < 3 * nc) continue;
      // ceiling: the model may not exceed the 90th percentile of what was measured at this radius
      var b = new Float64Array(vals.length);
      for (var i = 0; i < vals.length; ++i) b[i] = vals[i];
      var cap = GB.select(b, vals.length, Math.floor(0.9 * (vals.length - 1)));
      var r = GB.lsqRobust(cols, y, nn, ok, 6, 2.5, 4.0);
      for (var t = 0; t < nth; ++t) {
         var s = 0;
         for (var c = 0; c < nc; ++c) s += r.coef[c] * F[c][t];
         out[t * nr + j] = Math.min(s, Math.max(cap, 0));
      }
   }
   return out;
};

GB.addHalo = function (cols, cx, cy, x0, y0, w, h, pows) {
   var n = w * h, hs = [];
   for (var k = 0; k < pows.length; ++k) { hs.push(new Float64Array(n)); hs.push(new Float64Array(n)); hs.push(new Float64Array(n)); }
   for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
      var dx = x0 + i - cx, dy = y0 + j - cy, rr = Math.max(Math.sqrt(dx * dx + dy * dy), 1), q = j * w + i;
      for (var k = 0; k < pows.length; ++k) {
         var hv = Math.pow(rr / 20.0, -pows[k]);
         hs[3 * k][q] = hv; hs[3 * k + 1][q] = hv * dx / rr; hs[3 * k + 2][q] = hv * dy / rr;
      }
   }
   for (var k = 0; k < hs.length; ++k) cols.push(hs[k]);
};

// Mask compact (point-like) sources in a window so they do not bias the ghost estimate.
// L: Float64Array(w*h) luminance; sigL: luminance noise; scale: image scale.
GB.compactMask = function (L, w, h, sigL, scale, centerX, centerY, minDist) {
   var n = w * h, B = Math.max(5, Math.round(10 * scale));
   // box blur via integral image
   var I = new Float64Array((w + 1) * (h + 1));
   for (var j = 0; j < h; ++j) {
      var rs = 0;
      for (var i = 0; i < w; ++i) { rs += L[j * w + i]; I[(j + 1) * (w + 1) + i + 1] = I[j * (w + 1) + i + 1] + rs; }
   }
   var hp = new Float64Array(n);
   for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
      var a0 = Math.max(i - B, 0), a1 = Math.min(i + B + 1, w), b0 = Math.max(j - B, 0), b1 = Math.min(j + B + 1, h);
      var sum = I[b1 * (w + 1) + a1] - I[b0 * (w + 1) + a1] - I[b1 * (w + 1) + a0] + I[b0 * (w + 1) + a0];
      hp[j * w + i] = L[j * w + i] - sum / ((a1 - a0) * (b1 - b0));
   }
   var mask = new Uint8Array(n), thr = 10 * sigL, D = Math.max(2, Math.round(3 * scale)), found = 0;
   var dirs = [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, -1], [1, -1], [-1, 1]];
   for (var j = D; j < h - D; ++j) for (var i = D; i < w - D; ++i) {
      var q = j * w + i, v = hp[q];
      if (v < thr) continue;
      var isMax = true;
      for (var dj = -1; dj <= 1 && isMax; ++dj) for (var di = -1; di <= 1; ++di)
         if ((di || dj) && hp[q + dj * w + di] > v) { isMax = false; break; }
      if (!isMax) continue;
      if (Math.sqrt((i - centerX) * (i - centerX) + (j - centerY) * (j - centerY)) < minDist) continue;
      // compactness: brightness must fall in every direction (a ring rim does not, along the rim)
      var worst = 0;
      for (var k = 0; k < 8; ++k) {
         var f = (k < 4) ? 1 : 0.7071;
         var di2 = Math.round(dirs[k][0] * D * f), dj2 = Math.round(dirs[k][1] * D * f);
         var vv = hp[(j + dj2) * w + i + di2];
         if (vv > worst) worst = vv;
      }
      if (worst > 0.6 * v) continue;
      var rm = Math.min(scale * (3 + 2 * Math.log(v / thr) / Math.LN2), 18 * scale);
      var rmi = Math.ceil(rm);
      for (var dj = -rmi; dj <= rmi; ++dj) for (var di = -rmi; di <= rmi; ++di) {
         var ii = i + di, jj = j + dj;
         if (ii < 0 || jj < 0 || ii >= w || jj >= h) continue;
         if (di * di + dj * dj <= rm * rm) mask[jj * w + ii] = 1;
      }
      ++found;
   }
   return { mask: mask, count: found };
};

// Analyses one star: halo/background fit, neighbour masking, polar excess per channel.
// Returns { ok: false, gateRatio } when the surroundings are too structured (unless ignoreGate),
// otherwise { ok: true, gateRatio, pols, polL, tap, starCol, masked, ... }.
GB.analyseStar = function (ctx, s, work, ignoreGate) {
   var p = ctx.p, W = ctx.W, H = ctx.H, nc = ctx.nc, R = p.rOut, stars = ctx.stars, excl = ctx.excl;
   var x0 = Math.floor(Math.max(s.x - R, 0)), x1 = Math.floor(Math.min(s.x + R + 1, W));
   var y0 = Math.floor(Math.max(s.y - R, 0)), y1 = Math.floor(Math.min(s.y + R + 1, H));
   var w = x1 - x0, h = y1 - y0, n = w * h;
   var nbrs = [];
   for (var oi = 0; oi < stars.length; ++oi) {
      var o = stars[oi];
      if (o === s) continue;
      if (Math.sqrt((o.x - s.x) * (o.x - s.x) + (o.y - s.y) * (o.y - s.y)) <= R + p.annOut) nbrs.push(o);
   }
   var cols = [];
   for (var k = 0; k < 6; ++k) cols.push(new Float64Array(n));
   var r = new Float64Array(n), ex = new Uint8Array(n), ann = new Uint8Array(n), fit = new Uint8Array(n);
   for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
      var q = j * w + i, dx = x0 + i - s.x, dy = y0 + j - s.y, u = dx / 100, v = dy / 100;
      cols[0][q] = 1; cols[1][q] = u; cols[2][q] = v; cols[3][q] = u * u; cols[4][q] = u * v; cols[5][q] = v * v;
      r[q] = Math.sqrt(dx * dx + dy * dy);
      ex[q] = excl[(y0 + j) * W + x0 + i];
      for (var oi = 0; oi < nbrs.length; ++oi) {
         var ox = x0 + i - nbrs[oi].x, oy = y0 + j - nbrs[oi].y;
         if (Math.sqrt(ox * ox + oy * oy) < p.annIn) ex[q] = 1;
      }
      ann[q] = (r[q] >= p.annIn && r[q] <= p.annOut) ? 1 : 0;
      fit[q] = (r[q] <= R && !ex[q] && !ann[q]) ? 1 : 0;
   }
   GB.addHalo(cols, s.x, s.y, x0, y0, w, h, p.haloPows);
   for (var oi = 0; oi < nbrs.length; ++oi) GB.addHalo(cols, nbrs[oi].x, nbrs[oi].y, x0, y0, w, h, p.haloPows);
   var m = cols.length;
   var Lw = new Float64Array(n);
   for (var q = 0; q < n; ++q) {
      var xq = q % w, yq = (q - xq) / w, vv = 0;
      for (var c = 0; c < nc; ++c) vv += work[c][(y0 + yq) * W + x0 + xq];
      Lw[q] = vv;
   }
   var gate = GB.lsqRobust(cols, Lw, n, fit, p.iters, p.clipHi, p.clipLo);
   var res = { x0: x0, y0: y0, w: w, h: h, r: r, gateRatio: gate.sigma / ctx.sig };
   if (gate.sigma > p.gateK * ctx.sig && !ignoreGate) return res;
   var cm = GB.compactMask(Lw, w, h, ctx.sigL, p.scale, s.x - x0, s.y - y0, p.annIn * 0.75);
   res.masked = cm.count;
   var exE = new Uint8Array(n);
   for (var q = 0; q < n; ++q) exE[q] = ex[q] | cm.mask[q];
   var nth = ctx.nth, nr = ctx.nr, thStep = ctx.thStep, npol = nth * nr;
   var PX = new Float64Array(npol), PY = new Float64Array(npol), bad = new Uint8Array(npol);
   for (var t = 0; t < nth; ++t) for (var jr = 0; jr < nr; ++jr) {
      var k = t * nr + jr, rad = p.annIn + jr * p.rStep, th = t * thStep;
      var px = s.x - x0 + rad * Math.cos(th), py = s.y - y0 + rad * Math.sin(th);
      PX[k] = px; PY[k] = py;
      if (px < 0 || py < 0 || px > w - 1 || py > h - 1) { bad[k] = 1; continue; }
      var ix = Math.floor(px), iy = Math.floor(py), fx = px - ix, fy = py - iy;
      var ix1 = Math.min(ix + 1, w - 1), iy1 = Math.min(iy + 1, h - 1);
      var e = (1 - fx) * (1 - fy) * exE[iy * w + ix] + fx * (1 - fy) * exE[iy * w + ix1] + (1 - fx) * fy * exE[iy1 * w + ix] + fx * fy * exE[iy1 * w + ix1];
      if (e > 0.01) bad[k] = 1;
   }
   // rIn0: just outside the star core (star colour and star-light reference)
   var rIn0 = Math.max(p.coreMin, p.coreFactor * s.rs) + 1;
   var tap = new Float64Array(n), T = p.edgeTaper, Ti = Math.ceil(T) + 1;
   for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
      var q = j * w + i;
      if (r[q] < p.annIn || r[q] > p.annOut) continue;
      var t1 = Math.min(Math.max(Math.min(r[q] - p.annIn, p.annOut - r[q]) / T, 0), 1);
      var dmin = 1e9;
      if (ex[q]) dmin = 0;
      else for (var dj = -Ti; dj <= Ti; ++dj) for (var di = -Ti; di <= Ti; ++di) {
         var ii = i + di, jj = j + dj;
         if (ii < 0 || jj < 0 || ii >= w || jj >= h) continue;
         if (ex[jj * w + ii]) { var dd = Math.sqrt(di * di + dj * dj); if (dd < dmin) dmin = dd; }
      }
      tap[q] = t1 * Math.min(Math.max(dmin / T, 0), 1);
   }
   var pols = [], polL = new Float64Array(npol);
   for (var c = 0; c < nc; ++c) {
      var yc = new Float64Array(n);
      for (var q = 0; q < n; ++q) { var xq = q % w, yq = (q - xq) / w; yc[q] = work[c][(y0 + yq) * W + x0 + xq]; }
      var fr = GB.lsqRobust(cols, yc, n, fit, p.iters, p.clipHi, p.clipLo);
      var E = new Float64Array(n);
      for (var q = 0; q < n; ++q) {
         var mv = 0;
         for (var a = 0; a < m; ++a) mv += fr.coef[a] * cols[a][q];
         E[q] = yc[q] - mv;
      }
      var pol = new Float64Array(npol);
      for (var k = 0; k < npol; ++k) {
         if (bad[k]) { pol[k] = NaN; continue; }
         var px = PX[k], py = PY[k], ix = Math.floor(px), iy = Math.floor(py), fx = px - ix, fy = py - iy;
         var ix1 = Math.min(ix + 1, w - 1), iy1 = Math.min(iy + 1, h - 1);
         pol[k] = (1 - fx) * (1 - fy) * E[iy * w + ix] + fx * (1 - fy) * E[iy * w + ix1] + (1 - fx) * fy * E[iy1 * w + ix] + fx * fy * E[iy1 * w + ix1];
      }
      pols.push(pol);
      for (var k = 0; k < npol; ++k) polL[k] = (c == 0 ? 0 : polL[k]) + pol[k];
   }
   // star colour: inner wing of the star (dominates any nebula there) minus the background of the outer window
   var starCol = [];
   for (var c = 0; c < nc; ++c) {
      var inner = [], outer = [];
      for (var q = 0; q < n; ++q) {
         var xq = q % w, yq = (q - xq) / w, gi = (y0 + yq) * W + x0 + xq;
         if (r[q] >= rIn0 && r[q] < p.annIn && !ctx.sat[gi]) inner.push(work[c][gi]);
         else if (r[q] >= R - 20 && r[q] <= R) outer.push(work[c][gi]);
      }
      starCol.push(inner.length > 10 && outer.length > 10 ? GB.median(inner, inner.length) - GB.median(outer, outer.length) : 0);
   }
   res.pols = pols; res.polL = polL; res.tap = tap; res.exE = exE; res.rIn0 = rIn0; res.starCol = starCol; res.ok = true;
   return res;
};

// Ghost colour split chosen so that the corrected disk keeps the colour of the star light just
// outside the core (measured above the local background at the edge of the window). A wrong
// split shows as a coloured ring after correction, or as a hole in a coloured nebula. The
// luminance of the ghost is unchanged; each channel stays within [lo, hi] x its predicted share.
GB.continuityColour = function (an, s, work, shape, kc, ctx, lo, hi) {
   var p = ctx.p, W = ctx.W, nc = ctx.nc, nth = ctx.nth, nr = ctx.nr, thStep = ctx.thStep;
   var w = an.w, h = an.h, r = an.r, ex = an.exE, tap = an.tap, R = p.rOut, band = 4;
   var vin = [], vbg = [], sec = [];
   for (var c = 0; c < nc; ++c) { vin.push([]); vbg.push([]); sec.push([[], [], [], [], [], [], [], []]); }
   for (var q = 0; q < w * h; ++q) {
      if (ex[q]) continue;
      var gi = (an.y0 + Math.floor(q / w)) * W + an.x0 + q % w;
      if (ctx.sat[gi]) continue;
      var dst = r[q] >= an.rIn0 && r[q] < an.rIn0 + band ? vin : r[q] >= R - 10 && r[q] <= R ? vbg : null;
      if (dst) for (var c = 0; c < nc; ++c) dst[c].push(work[c][gi]);
      if (dst === vbg) {
         var th8 = Math.atan2(an.y0 + Math.floor(q / w) - s.y, an.x0 + q % w - s.x); if (th8 < 0) th8 += 2 * Math.PI;
         var k8 = Math.min(Math.floor(th8 / (Math.PI / 4)), 7);
         for (var c = 0; c < nc; ++c) sec[c][k8].push(work[c][gi]);
      }
   }
   if (vin[0].length < 20 || vbg[0].length < 20) return false;
   // structured background (nebula gradient across the window): the star-light reference is not
   // reliable, the caller keeps its own colour estimate
   var gmean = 0, gcnt = 0;
   for (var k = 0; k < nth * nr; ++k) if (shape[k] > 0) { gmean += shape[k]; ++gcnt; }
   gmean = gcnt ? gmean / gcnt : 0;
   for (var c = 0; c < nc; ++c) {
      var lo8 = Infinity, hi8 = -Infinity;
      for (var k8 = 0; k8 < 8; ++k8) {
         if (sec[c][k8].length < 10) continue;
         var m8 = GB.median(sec[c][k8], sec[c][k8].length);
         if (m8 < lo8) lo8 = m8; if (m8 > hi8) hi8 = m8;
      }
      if (hi8 - lo8 > p.bgSpread * kc[c] * gmean) return false;
   }
   var bg = [], f = [], sf = 0;
   for (var c = 0; c < nc; ++c) {
      bg.push(GB.median(vbg[c], vbg[c].length));
      f.push(GB.median(vin[c], vin[c].length) - bg[c]);
      if (f[c] < 0) return false;
      sf += f[c];
   }
   if (sf <= 0) return false;
   for (var c = 0; c < nc; ++c) f[c] /= sf;
   // least squares on kc: (pixel - background) - kc * ghost = f * (luminance above background - ghost)
   var num = [], den = 0;
   for (var c = 0; c < nc; ++c) num.push(0);
   for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
      var q = j * w + i;
      if (tap[q] < 0.5) continue;
      var th = Math.atan2(an.y0 + j - s.y, an.x0 + i - s.x); if (th < 0) th += 2 * Math.PI;
      var ri = Math.max((r[q] - p.annIn) / p.rStep, 0);
      if (ri > nr - 1) continue;
      var v = shape[(Math.min(Math.round(th / thStep), nth - 1) % nth) * nr + Math.round(ri)];
      if (!(v > 0)) continue;
      var gi = (an.y0 + j) * W + an.x0 + i, DL = 0;
      for (var c = 0; c < nc; ++c) DL += work[c][gi] - bg[c];
      for (var c = 0; c < nc; ++c) num[c] += v * (work[c][gi] - bg[c] - f[c] * (DL - v));
      den += v * v;
   }
   if (den <= 0) return false;
   var kl = [], ks = 0;
   for (var c = 0; c < nc; ++c) { var k = Math.min(Math.max(num[c] / den, lo * kc[c]), hi * kc[c]); kl.push(k); ks += k; }
   if (ks <= 0) return false;
   for (var c = 0; c < nc; ++c) kc[c] = kl[c] / ks;
   return true;
};

// Robust azimuthal profile (median over angle of valid samples) per radius index
GB.azProfile = function (pol, nth, nr) {
   var out = new Float64Array(nr), buf = new Float64Array(nth);
   for (var j = 0; j < nr; ++j) {
      var k = 0;
      for (var t = 0; t < nth; ++t) { var v = pol[t * nr + j]; if (!isNaN(v)) buf[k++] = v; }
      out[j] = k >= nth * 0.3 ? GB.median(buf, k) : NaN;
   }
   return out;
};

// Fit the ghost of one star with a FIXED radial shape T: g(t,j) = T[j] * (a + bc cos + bs sin)
GB.fitTemplate = function (polL, T, nth, nr, thStep, maxLop) {
   var cols = [new Float64Array(nth * nr), new Float64Array(nth * nr), new Float64Array(nth * nr)];
   var y = new Float64Array(nth * nr), ok = new Uint8Array(nth * nr), cnt = 0;
   for (var t = 0; t < nth; ++t) {
      var ct = Math.cos(t * thStep), st = Math.sin(t * thStep);
      for (var j = 0; j < nr; ++j) {
         var k = t * nr + j, v = polL[k];
         cols[0][k] = T[j]; cols[1][k] = T[j] * ct; cols[2][k] = T[j] * st;
         if (!isNaN(v) && T[j] > 0) { y[k] = v; ok[k] = 1; ++cnt; }
      }
   }
   if (cnt < 50) return null;
   var fr = GB.lsqRobust(cols, y, nth * nr, ok, 6, 2.5, 2.5), f = fr.coef;
   var a = Math.max(f[0], 0), bc = f[1], bs = f[2], lop = Math.sqrt(bc * bc + bs * bs);
   if (lop > maxLop * a) { var sc = lop > 0 ? maxLop * a / lop : 0; bc *= sc; bs *= sc; }
   return { a: a, bc: bc, bs: bs, sigma: fr.sigma };
};

// Template stretched radially by alpha (ghost edge slightly larger/smaller for this star)
GB.stretchTemplate = function (T, nr, annIn, rStep, alpha) {
   var out = new Float64Array(nr);
   for (var j = 0; j < nr; ++j) {
      var r = (annIn + j * rStep) / alpha, x = (r - annIn) / rStep;
      if (x <= 0) { out[j] = T[0]; continue; }
      if (x >= nr - 1) { out[j] = T[nr - 1]; continue; }
      var i0 = Math.floor(x), f = x - i0;
      out[j] = (1 - f) * T[i0] + f * T[i0 + 1];
   }
   return out;
};

// Main entry. chans: array of Float32Array (normalized). selectFn(star) -> bool (optional). log(msg)
// Returns { ghost: array of Float32Array, done, skipped, total, templateSource, stars }
// opts (optional): { dryRun: bool, extra: [{x, y}] manual stars, forceFn: function(star) -> bool }
//   dryRun: classify only (star.status = 'gate' | 'noghost' | 'ok', star.conf), nothing subtracted
//   extra:  stars the user added; snapped to the brightest pixel within 12 px
//   forceFn: stars to correct even if the nebula safety gate would skip them
GB.process = function (chans, W, H, params, selectFn, log, opts) {
   opts = opts || {};
   var p = GB.scaled(params), nc = chans.length, N = W * H;
   log = log || function () {};
   // opts.progress(stage, i, n): called before each star; returning true aborts the run
   var tick = function (stage, i, n) {
      if (opts.progress && opts.progress(stage, i, n)) throw new Error("Aborted by the user.");
   };
   tick("Detecting saturated stars", 0, 1);
   var det = GB.detect(chans, W, H, p), stars = det.stars, sat = det.sat, nDet = stars.length;
   // the largest saturated cores come first: beyond maxStars, the rest is left out
   if (nDet > p.maxStars) stars = stars.slice(0, p.maxStars);
   if (opts.extra) {
      for (var ei = 0; ei < opts.extra.length; ++ei) {
         var ex0 = opts.extra[ei];
         if (!ex0 || !isFinite(ex0.x) || !isFinite(ex0.y)) continue;
         var bx = Math.round(ex0.x), by = Math.round(ex0.y);
         if (bx < 0 || by < 0 || bx >= W || by >= H) continue;
         // snap to the brightest pixel nearby (skipped when the caller already snapped)
         var SN = ex0.snapped ? 0 : Math.round(12 * p.scale), bv = -1, mx = bx, my = by;
         for (var yy = Math.max(0, by - SN); yy <= Math.min(H - 1, by + SN); ++yy)
            for (var xx = Math.max(0, bx - SN); xx <= Math.min(W - 1, bx + SN); ++xx) {
               var lv = 0; for (var c = 0; c < nc; ++c) lv += chans[c][yy * W + xx];
               if (lv > bv) { bv = lv; mx = xx; my = yy; }
            }
         var dup = false;
         for (var si = 0; si < stars.length; ++si)
            if (Math.abs(stars[si].x - mx) < 4 && Math.abs(stars[si].y - my) < 4) { dup = true; break; }
         if (!dup) stars.push({ x: mx, y: my, rs: 2, area: 0, manual: true });
      }
   }
   var sig = 0, sig2 = 0;
   for (var c = 0; c < nc; ++c) { var sc = GB.noise(chans[c], W, H); sig += sc; sig2 += sc * sc; }
   var excl = new Uint8Array(N);
   for (var i = 0; i < N; ++i) excl[i] = sat[i];
   for (var si = 0; si < stars.length; ++si) {
      var s = stars[si], rc = Math.max(p.coreMin, p.coreFactor * s.rs);
      for (var y = Math.max(0, Math.floor(s.y - rc)); y <= Math.min(H - 1, Math.ceil(s.y + rc)); ++y)
         for (var x = Math.max(0, Math.floor(s.x - rc)); x <= Math.min(W - 1, Math.ceil(s.x + rc)); ++x)
            if (Math.sqrt((x - s.x) * (x - s.x) + (y - s.y) * (y - s.y)) < rc) excl[y * W + x] = 1;
   }
   var ctx = { p: p, W: W, H: H, nc: nc, stars: stars, excl: excl, sat: sat, sig: sig, sigL: Math.sqrt(sig2),
               thStep: p.thStepDeg * Math.PI / 180, nth: Math.round(360 / p.thStepDeg),
               nr: Math.floor((p.annOut - p.annIn) / p.rStep + 1e-6) + 1 };
   var nth = ctx.nth, nr = ctx.nr, thStep = ctx.thStep;
   log("Detected " + nDet + " saturated stars" + (nDet > p.maxStars ? "; the " + p.maxStars + " largest are kept." : "."));

   // ---- pass 1: learn the ghost radial shape on clean (dark-sky) stars ----
   var profs = [];
   for (var si = 0; si < stars.length && profs.length < 25; ++si) {
      if (stars[si].manual) continue;
      tick("Learning the ghost shape", si, stars.length);
      var an = GB.analyseStar(ctx, stars[si], chans);
      if (!an.ok || an.gateRatio > 3) continue;
      var pr = GB.azProfile(an.polL, nth, nr), mx = 0, nanCount = 0;
      for (var j = 0; j < nr; ++j) { if (isNaN(pr[j])) ++nanCount; else if (pr[j] > mx) mx = pr[j]; }
      if (nanCount > nr * 0.2 || mx <= 6 * ctx.sigL / Math.sqrt(nth)) continue;
      for (var j = 0; j < nr; ++j) pr[j] = isNaN(pr[j]) ? NaN : Math.max(pr[j], 0) / mx;
      profs.push(pr);
   }
   var T = new Float64Array(nr), src;
   if (profs.length >= 3) {
      var buf = [];
      for (var j = 0; j < nr; ++j) {
         buf.length = 0;
         for (var k = 0; k < profs.length; ++k) if (!isNaN(profs[k][j])) buf.push(profs[k][j]);
         T[j] = buf.length ? GB.median(buf, buf.length) : 0;
      }
      var Ts = new Float64Array(nr);
      for (var j = 0; j < nr; ++j) {
         var acc = 0, cn = 0;
         for (var d = -2; d <= 2; ++d) { var jj = Math.min(Math.max(j + d, 0), nr - 1); acc += T[jj]; ++cn; }
         Ts[j] = acc / cn;
      }
      T = Ts;
      src = "learned from " + profs.length + " clean stars";
   } else {
      for (var j = 0; j < nr; ++j) T[j] = GB.TEMPLATE[Math.min(j, GB.TEMPLATE.length - 1)];
      src = "built-in Draco profile (too few clean stars: " + profs.length + ")";
   }
   var tmax = 0;
   for (var j = 0; j < nr; ++j) if (T[j] > tmax) tmax = T[j];
   for (var j = 0; j < nr; ++j) T[j] = tmax > 0 ? Math.max(T[j], 0) / tmax : 0;
   log("Ghost shape: " + src + ".");
   // ---- pass 1b: ghost-to-star ratio on clean stars (physical ceiling for stars in nebula) ----
   var ratios = [], tr = [];
   for (var c = 0; c < nc; ++c) tr.push([]);
   for (var si = 0; si < stars.length && ratios.length < 25; ++si) {
      if (stars[si].manual) continue;
      tick("Learning the ghost intensity", si, stars.length);
      var an = GB.analyseStar(ctx, stars[si], chans);
      if (!an.ok || an.gateRatio > 3) continue;
      var hL = 0; for (var c = 0; c < nc; ++c) hL += an.starCol[c];
      var gf = GB.fitTemplate(an.polL, T, nth, nr, thStep, p.maxLop);
      if (!gf || hL <= 0) continue;
      ratios.push(gf.a / hL);
      // colour transfer: ghost colour relative to the star colour, per channel
      for (var c = 0; c < nc; ++c) {
         var gc = GB.fitTemplate(an.pols[c], T, nth, nr, thStep, p.maxLop);
         if (gc && an.starCol[c] > 0 && gf.a > 0) tr[c].push((gc.a / gf.a) / (an.starCol[c] / hL));
      }
   }
   var ctr = [];
   for (var c = 0; c < nc; ++c) ctr.push(tr[c].length >= 3 ? Math.max(GB.median(tr[c], tr[c].length), 0.05) : 1.0);
   log("Ghost colour vs star colour: " + ctr.map(function (v) { return v.toFixed(2); }).join(" / ") + ".");
   var rho = ratios.length >= 3 ? GB.median(ratios, ratios.length) : GB.RATIO;
   log("Ghost/star ratio: " + rho.toFixed(3) + (ratios.length >= 3 ? " (from " + ratios.length + " clean stars)" : " (built-in)") + ".");

   // ---- pass 2: fit and subtract the ghost of every selected star ----
   // dry-run: nothing is subtracted, so the input is read directly and no output is allocated
   var work = [], ghost = [];
   for (var c = 0; c < nc; ++c) {
      work.push(opts.dryRun ? chans[c] : new Float32Array(chans[c]));
      if (!opts.dryRun) ghost.push(new Float32Array(N));
   }
   var nDone = 0, nSkip = 0;
   for (var si = 0; si < stars.length; ++si) {
      var s = stars[si];
      if (selectFn && !selectFn(s)) continue;
      var tag = "(" + Math.round(s.x) + ", " + Math.round(s.y) + ")";
      tick(opts.dryRun ? "Analysing stars" : "Correcting stars", si, stars.length);
      var force = !!(s.manual || (opts.forceFn && opts.forceFn(s)));
      var an = GB.analyseStar(ctx, s, work, force);
      if (!an.ok && opts.dryRun) { s.status = "gate"; s.conf = 0; continue; }
      if (!an.ok) {
         log("SKIP " + tag + ": residual " + Math.round(an.gateRatio) + "x noise (bright structured nebula)");
         ++nSkip; continue;
      }
      // predicted ghost colour: star colour x learned colour transfer
      var kc = [], ksum = 0;
      for (var c = 0; c < nc; ++c) { var v0 = Math.max(an.starCol[c], 0) * ctr[c]; kc.push(v0); ksum += v0; }
      for (var c = 0; c < nc; ++c) kc[c] = ksum > 0 ? kc[c] / ksum : 1 / nc;
      var polG = an.polL;
      // radial stretch of the learned shape for this star (best of a small grid)
      var gf = null, Ts = T;
      for (var al = 0.85; al <= 1.151; al += 0.025) {
         var Ta = GB.stretchTemplate(T, nr, p.annIn, p.rStep, al);
         var g1 = GB.fitTemplate(polG, Ta, nth, nr, thStep, p.maxLop);
         if (g1 && (!gf || g1.sigma < gf.sigma)) { gf = g1; Ts = Ta; }
      }
      if (!gf || gf.a <= 0) {
         if (opts.dryRun) { s.status = "noghost"; s.conf = 0; continue; }
         log("OK   " + tag + ": no ghost detected"); ++nDone; continue;
      }
      // confidence: ghost amplitude vs residual structure around it. Weak ghost in structured
      // nebula -> regular round ring; strong clear ghost -> full shape freedom.
      var conf = gf.sigma > 0 ? gf.a / gf.sigma : 0;
      if (opts.dryRun) { s.status = "ok"; s.conf = conf; continue; }
      var fr0 = Math.min(Math.max((conf - p.confLo) / (p.confHi - p.confLo), 0), 1);
      // physical ceiling: a fit above it has absorbed surrounding signal (nebula, neighbour),
      // so the fit is discarded and the ghost falls back to the typical ghost of a star of this
      // brightness: a round ring, without free correction
      var hL = 0; for (var c = 0; c < nc; ++c) hL += an.starCol[c];
      var capA = p.ratioCap * rho * Math.max(hL, 0), capped = false;
      if (gf.a > capA) { gf.a = rho * Math.max(hL, 0); gf.bc = 0; gf.bs = 0; fr0 = 0; capped = true; }
      var lopNow = Math.sqrt(gf.bc * gf.bc + gf.bs * gf.bs), lopMax = p.maxLop * fr0 * gf.a;
      if (lopNow > lopMax) { var sl = lopNow > 0 ? lopMax / lopNow : 0; gf.bc *= sl; gf.bs *= sl; }
      var freeNow = p.freeBound * fr0;
      var shape = new Float64Array(nth * nr);
      for (var t = 0; t < nth; ++t) {
         var ang = gf.a + gf.bc * Math.cos(t * thStep) + gf.bs * Math.sin(t * thStep);
         for (var j = 0; j < nr; ++j) shape[t * nr + j] = Ts[j] * Math.max(ang, 0);
      }
      // bounded free correction: follows the real ghost (sharp edge, offset) but never
      // departs from the fitted shape by more than +/- p.freeBound of it
      if (freeNow > 0) {
         var resid = new Float64Array(nth * nr);
         for (var k = 0; k < nth * nr; ++k) resid[k] = polG[k] - shape[k];
         var dlt = GB.fourierRing(resid, nth, nr, thStep, p.nHarm);
         for (var k = 0; k < nth * nr; ++k) {
            var lim = freeNow * shape[k], dv = Math.min(Math.max(dlt[k], -lim), lim);
            shape[k] = Math.max(shape[k] + dv, 0);
         }
      }
      if (nc > 1 && !GB.continuityColour(an, s, work, shape, kc, ctx, p.colLo, p.colHi)) {
         // fallback: per-channel amplitude on the final shape, on the data with background and
         // halos removed, kept within [colLo, colHi] x predicted colour
         var kl = [], ksum2 = 0;
         for (var c = 0; c < nc; ++c) {
            var okk = new Uint8Array(nth * nr), yv = new Float64Array(nth * nr), cnt = 0;
            for (var k = 0; k < nth * nr; ++k) {
               if (isNaN(an.pols[c][k]) || shape[k] <= 0) continue;
               okk[k] = 1; yv[k] = an.pols[c][k]; ++cnt;
            }
            var kk = cnt > 50 ? GB.lsqRobust([shape], yv, nth * nr, okk, 5, 2.5, 2.5).coef[0] : kc[c];
            kl.push(Math.min(Math.max(kk, p.colLo * kc[c]), p.colHi * kc[c])); ksum2 += kl[c];
         }
         if (ksum2 > 0) for (var c = 0; c < nc; ++c) kc[c] = kl[c] / ksum2;
      }
      var peaks = [], x0 = an.x0, y0 = an.y0, w = an.w, h = an.h, r = an.r, tap = an.tap;
      for (var c = 0; c < nc; ++c) {
         var pk = 0, kk = kc[c];
         for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
            var q = j * w + i;
            if (tap[q] <= 0) continue;
            var dx = x0 + i - s.x, dy = y0 + j - s.y;
            var th = Math.atan2(dy, dx); if (th < 0) th += 2 * Math.PI;
            var ti = th / thStep, ri = (r[q] - p.annIn) / p.rStep;
            ti = Math.min(Math.max(ti, 0), nth); ri = Math.min(Math.max(ri, 0), nr - 1);
            var it0 = Math.floor(ti), ir0 = Math.floor(ri), ft = ti - it0, fr2 = ri - ir0;
            var it1 = Math.min(it0 + 1, nth), ir1 = Math.min(ir0 + 1, nr - 1);
            var ta = (it0 % nth) * nr, tb = (it1 % nth) * nr;
            var v = (1 - ft) * (1 - fr2) * shape[ta + ir0] + ft * (1 - fr2) * shape[tb + ir0] + (1 - ft) * fr2 * shape[ta + ir1] + ft * fr2 * shape[tb + ir1];
            var g = kk * v * tap[q];
            if (g > 0) {
               var gi = (y0 + j) * W + x0 + i;
               ghost[c][gi] += g; work[c][gi] -= g;
               if (g > pk) pk = g;
            }
         }
         peaks.push(pk);
      }
      ++nDone;
      log("OK   " + tag + ": ghost removed (confidence " + conf.toFixed(1) + "; " + an.masked + " neighbour stars masked" + (capped ? ", capped" : "") + "), peak " +
          peaks.map(function (v) { return v.toExponential(2); }).join(" / "));
   }
   return { ghost: ghost, done: nDone, skipped: nSkip, total: stars.length, templateSource: src, stars: stars };
};

// Fallback values used when the image has fewer than 3 clean stars.
// RATIO: ghost amplitude (sum of channels, at the profile maximum) divided by the star level
//   (sum of channels) measured on the inner wing; median of the values measured on four
//   native Draco duo-band stacks (0.016 to 0.064).
// TEMPLATE: radial ghost profile, normalised to 1 at its maximum, index j = radius annIn + j * rStep,
//   measured on a native Draco duo-band stack.
GB.RATIO = 0.03;
GB.TEMPLATE = [0.584, 0.605, 0.623, 0.635, 0.624, 0.607, 0.592, 0.585, 0.585, 0.604, 0.620, 0.650, 0.684, 0.729, 0.765, 0.794, 0.814, 0.826, 0.829, 0.847, 0.877, 0.913, 0.948, 0.985, 1.000, 0.998, 0.989, 0.967, 0.932, 0.894, 0.842, 0.777, 0.710, 0.648, 0.582, 0.523, 0.465, 0.408, 0.345, 0.288, 0.238, 0.200, 0.169, 0.145, 0.124, 0.107, 0.086, 0.068, 0.055, 0.051, 0.046, 0.045, 0.044, 0.042, 0.038, 0.035, 0.030, 0.026, 0.023, 0.018, 0.016];
// ===================== end of core =====================

// ---- Automated mode: scale from image metadata ----
// Native Draco (12 MP binned mode): 2.394 um pixels, 340 mm focal length -> 1.452 arcsec/px.
GB.REF_ARCSEC_PER_PX = 1.452;
GB.REF_WIDTH = 4096;

// kw: array of {name, value} (strings). Returns { scale, source, fromMetadata }
GB.autoScale = function (kw, W, H) {
   var k = {};
   for (var i = 0; i < kw.length; ++i) {
      var v = String(kw[i].value).replace(/'/g, "").trim();
      k[String(kw[i].name).trim().toUpperCase()] = v;
   }
   function num(n) { if (!k.hasOwnProperty(n)) return NaN; var x = parseFloat(k[n]); return isFinite(x) ? x : NaN; }
   var ps = NaN, src = "";
   var cd11 = num("CD1_1"), cd21 = num("CD2_1"), cd12 = num("CD1_2"), cd22 = num("CD2_2");
   if (isFinite(cd11) && isFinite(cd21) && isFinite(cd12) && isFinite(cd22)) {
      ps = 3600 * Math.sqrt(Math.abs(cd11 * cd22 - cd12 * cd21));
      src = "astrometric solution";
   } else if (isFinite(num("CDELT1"))) {
      ps = 3600 * Math.abs(num("CDELT1")); src = "astrometric solution";
   } else if (isFinite(num("XPIXSZ")) && isFinite(num("FOCALLEN")) && num("FOCALLEN") > 0) {
      ps = 206.265 * num("XPIXSZ") / num("FOCALLEN"); src = "pixel size and focal length";
   }
   // The image size is only used to detect an enlarged (drizzled / upsampled) image, which can keep
   // the original XPIXSZ in its header. A smaller image is most often a crop, so it never lowers the scale.
   var ratio = Math.max(W, H * 4 / 3) / GB.REF_WIDTH;
   var sizeScale = ratio > 1.3 ? Math.round(ratio * 2) / 2 : 1.0;
   var scale;
   var astrometric = src.indexOf("astrometric") == 0;
   if (isFinite(ps) && ps > 0.05 && ps < 30 && (astrometric || sizeScale == 1.0)) {
      scale = GB.REF_ARCSEC_PER_PX / ps;
      src += " (" + ps.toFixed(2) + "\"/px)";
   } else if (isFinite(ps) && ps > 0.05 && ps < 30) {
      scale = sizeScale;
      src = "image size (" + W + " x " + H + "), larger than the pixel size in the header suggests";
   } else {
      scale = sizeScale;
      src = "image size (" + W + " x " + H + ")";
   }
   scale = Math.min(Math.max(scale, 0.25), 4);
   return { scale: Math.round(scale * 100) / 100, source: src, fromMetadata: src.indexOf("image size") != 0 };
};

// ===================== PixInsight interface =====================

function GB_readChannels(img) {
   var W = img.width, H = img.height, nc = img.numberOfChannels, ch = [];
   for (var c = 0; c < nc; ++c) {
      var a = new Float32Array(W * H);
      try {
         img.getSamples(a, new Rect(0, 0, W, H), c);
      } catch (e) {   // older PixInsight: plain Array only
         var tmp = [];
         img.getSamples(tmp, new Rect(0, 0, W, H), c);
         for (var i = 0; i < W * H; ++i) a[i] = tmp[i];
      }
      // non-finite pixels (NaN / Inf) are replaced by 0 in this working copy only
      var bad = 0;
      for (var i = 0; i < W * H; ++i) if (!isFinite(a[i])) { a[i] = 0; ++bad; }
      if (bad > 0)
         console.warningln("Channel " + c + ": " + bad + " non-finite pixel(s) treated as 0 " +
                           "(the source image is not modified).");
      ch.push(a);
   }
   return ch;
}

// HTML-escape text inserted in rich-text labels
function GB_esc(s) {
   return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}

function GB_validId(id) {
   var s = id.replace(/[^A-Za-z0-9_]/g, "_");
   if (!/^[A-Za-z_]/.test(s)) s = "_" + s;
   return s;
}

function GB_makeWindow(srcWin, id, chans, W, H) {
   var nc = chans.length;
   var win = new ImageWindow(W, H, nc, 32, true, nc >= 3, GB_validId(id));
   var v = win.mainView;
   v.beginProcess(UndoFlag_NoSwapFile);
   for (var c = 0; c < nc; ++c) {
      try {
         v.image.setSamples(chans[c], new Rect(0, 0, W, H), c);
      } catch (e) {
         var tmp = new Array(W * H);
         for (var i = 0; i < W * H; ++i) tmp[i] = chans[c][i];
         v.image.setSamples(tmp, new Rect(0, 0, W, H), c);
      }
   }
   v.endProcess();
   try { win.keywords = srcWin.keywords; }
   catch (e) { console.warningln("Could not copy FITS keywords to " + win.mainView.id + "."); }
   try { win.copyAstrometricSolution(srcWin); }
   catch (e) { /* not available in older PixInsight versions, or no solution: optional */ }
   return win;
}

function GB_keywords(win) {
   var out = [];
   try { var kw = win.keywords; for (var i = 0; i < kw.length; ++i) out.push({ name: kw[i].name, value: kw[i].value }); }
   catch (e) { /* no keywords: the scale falls back to the image size */ }
   return out;
}

// ---- preview data: block-averaged, auto-stretched (same formula as PixInsight's STF AutoStretch) ----
function GB_mtf(m, x) { return x <= 0 ? 0 : (x >= 1 ? 1 : ((m - 1) * x) / ((2 * m - 1) * x - m)); }

function GB_makePreview(chans, W, H, maxW, maxH) {
   var f = Math.max(1, Math.ceil(Math.max(W / maxW, H / maxH)));
   var w = Math.floor(W / f), h = Math.floor(H / f), nc = chans.length, small = [], stf = [];
   for (var c = 0; c < nc; ++c) {
      var a = new Float32Array(w * h), src = chans[c];
      for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
         var s = 0;
         for (var dj = 0; dj < f; ++dj) { var row = (j * f + dj) * W + i * f; for (var di = 0; di < f; ++di) s += src[row + di]; }
         a[j * w + i] = s / (f * f);
      }
      // auto stretch: shadows at median - 2.8 MAD, median mapped to 0.25
      var med = GB.median(a, a.length), dev = new Float64Array(a.length);
      for (var i = 0; i < a.length; ++i) dev[i] = Math.abs(a[i] - med);
      var mad = 1.4826 * GB.median(dev, dev.length), c0 = Math.min(Math.max(med - 2.8 * mad, 0), 0.99);
      var m = med > c0 ? GB_mtf(0.25, (med - c0) / (1 - c0)) : 0.5;
      for (var i = 0; i < a.length; ++i) a[i] = GB_mtf(m, (a[i] - c0) / (1 - c0));
      small.push(a);
      stf.push({ c0: c0, m: m });
   }
   return { f: f, w: w, h: h, small: small, stf: stf };
}

// Build a Bitmap from 3 Float32Arrays in [0,1] (Image.render when available, else setPixel)
function GB_toBitmap(rgb, w, h) {
   var bmp = null;
   try {
      var img = new Image(w, h, 3, ColorSpace_RGB, 32, SampleType_Real);
      for (var c = 0; c < 3; ++c) {
         try { img.setSamples(rgb[c], new Rect(0, 0, w, h), c); }
         catch (e1) { var t = new Array(w * h); for (var i = 0; i < w * h; ++i) t[i] = rgb[c][i]; img.setSamples(t, new Rect(0, 0, w, h), c); }
      }
      bmp = img.render();
      try { img.free(); } catch (e2) { /* freed by the garbage collector anyway */ }
   } catch (e) { bmp = null; }
   if (!bmp) {
      bmp = new Bitmap(w, h);
      for (var j = 0; j < h; ++j) for (var i = 0; i < w; ++i) {
         var q = j * w + i;
         var r = Math.round(255 * rgb[0][q]), g = Math.round(255 * rgb[1][q]), b = Math.round(255 * rgb[2][q]);
         bmp.setPixel(i, j, (0xff000000 | (r << 16) | (g << 8) | b) >>> 0);
      }
   }
   return bmp;
}

// brightest pixel (sum of channels) within rad px of (x, y)
function GB_snap(chans, W, H, x, y, rad) {
   var bx = Math.round(x), by = Math.round(y), best = -1, sx = bx, sy = by, nc = chans.length;
   for (var yy = Math.max(0, by - rad); yy <= Math.min(H - 1, by + rad); ++yy)
      for (var xx = Math.max(0, bx - rad); xx <= Math.min(W - 1, bx + rad); ++xx) {
         var v = 0; for (var c = 0; c < nc; ++c) v += chans[c][yy * W + xx];
         if (v > best) { best = v; sx = xx; sy = yy; }
      }
   return { x: sx, y: sy };
}

function GB_key(x, y) { return Math.round(x) + "," + Math.round(y); }

function GB_Dialog(win, params, chans) {
   this.__base__ = Dialog;
   this.__base__();
   var self = this;
   var img0 = win.mainView.image, W = img0.width, H = img0.height;
   this.params = params;
   this.chans = chans;
   this.auto = GB.autoScale(GB_keywords(win), W, H);
   this.params.manual = false;
   this.params.scale = this.auto.scale;
   this.marks = [];        // { x, y, cls: 'white'|'orange'|'manual', on, status, conf }
   this.override = {};     // key -> true/false chosen by the user
   this.extra = [];        // user-added stars { x, y } (image coordinates, snapped)

   this.effectiveParams = function () {
      if (self.params.manual) return self.params;
      var q = GB.defaults();
      q.scale = self.auto.scale;
      return q;
   };

   // ---------- left: zoomable, pannable, clickable preview ----------
   this.prev = GB_makePreview(chans, W, H, 880, 660);
   this.cw = 880; this.ch = 660;
   this.zFit = Math.max(W / this.cw, H / this.ch);
   this.zoom = this.zFit;
   this.vx = (W - this.cw * this.zoom) / 2;    // image coordinate of the canvas top-left corner
   this.vy = (H - this.ch * this.zoom) / 2;
   this.drag = null;                            // { x, y, vx, vy, moved }
   this.dragDx = 0; this.dragDy = 0;

   this.clampView = function () {
      var z = self.zoom;
      self.vx = Math.min(Math.max(self.vx, -self.cw * z / 2), W - self.cw * z / 2);
      self.vy = Math.min(Math.max(self.vy, -self.ch * z / 2), H - self.ch * z / 2);
   };

   // half-resolution stretched level, for intermediate zooms (computed once)
   this.half = null;
   this.makeHalf = function () {
      var w2 = Math.floor(W / 2), h2 = Math.floor(H / 2), nc = chans.length, lv = [];
      for (var c = 0; c < nc; ++c) {
         var a = new Float32Array(w2 * h2), src = chans[c], st = self.prev.stf[c];
         for (var j = 0; j < h2; ++j) for (var i = 0; i < w2; ++i) {
            var q0 = 2 * j * W + 2 * i;
            a[j * w2 + i] = GB_mtf(st.m, ((src[q0] + src[q0 + 1] + src[q0 + W] + src[q0 + W + 1]) / 4 - st.c0) / (1 - st.c0));
         }
         lv.push(a);
      }
      self.half = { w: w2, h: h2, data: lv };
   };

   // render the visible region from the best level: small preview (zoomed out),
   // half resolution (intermediate), full resolution with the same stretch (zoomed in)
   this.renderView = function () {
      var cw = self.cw, ch = self.ch, z = self.zoom, nc = chans.length, rgb = [];
      for (var c = 0; c < 3; ++c) rgb.push(new Float32Array(cw * ch));
      var P = self.prev, level = z >= P.f * 0.75 ? 0 : (z >= 1.5 ? 1 : 2);
      if (level == 1 && !self.half) self.makeHalf();
      for (var j = 0; j < ch; ++j) for (var i = 0; i < cw; ++i) {
         var q = j * cw + i, X = self.vx + (i + 0.5) * z, Y = self.vy + (j + 0.5) * z;
         if (X < 0 || Y < 0 || X >= W || Y >= H) { rgb[0][q] = rgb[1][q] = rgb[2][q] = 0.12; continue; }
         for (var c = 0; c < 3; ++c) {
            var cc = Math.min(c, nc - 1), v;
            if (level == 0) {
               var sx = Math.min(P.w - 1, Math.floor(X / P.f)), sy = Math.min(P.h - 1, Math.floor(Y / P.f));
               v = P.small[cc][sy * P.w + sx];
            } else if (level == 1) {
               var hx = Math.min(self.half.w - 1, Math.floor(X / 2)), hy = Math.min(self.half.h - 1, Math.floor(Y / 2));
               v = self.half.data[cc][hy * self.half.w + hx];
            } else {
               var st = P.stf[cc];
               v = GB_mtf(st.m, (chans[cc][Math.floor(Y) * W + Math.floor(X)] - st.c0) / (1 - st.c0));
            }
            rgb[c][q] = v;
         }
      }
      self.viewBmp = GB_toBitmap(rgb, cw, ch);
   };

   this.canvas = new Control(this);
   this.canvas.setFixedSize(this.cw, this.ch);
   this.canvas.toolTip = "<p>Mouse wheel: zoom. Drag: move around. Click a circle: include / exclude. Click a star without a circle: add it.</p>";
   this.markRadius = function () { return Math.min(Math.max(8, 46 * self.auto.scale / self.zoom), 200); };
   this.toScreen = function (x, y) { return { x: (x - self.vx) / self.zoom + self.dragDx, y: (y - self.vy) / self.zoom + self.dragDy }; };
   this.canvas.onPaint = function (x0, y0, x1, y1) {
      var g = new Graphics(this);
      g.fillRect(0, 0, self.cw, self.ch, new Brush(0xff1e1e1e));
      g.drawBitmap(self.dragDx, self.dragDy, self.viewBmp);
      g.antialiasing = true;
      var rad = self.markRadius();
      for (var i = 0; i < self.marks.length; ++i) {
         var mk = self.marks[i], col, wid = 1.5, sp = self.toScreen(mk.x, mk.y);
         if (sp.x < -rad || sp.y < -rad || sp.x > self.cw + rad || sp.y > self.ch + rad) continue;
         if (mk.cls == "manual") col = 0xff40ff60;
         else if (mk.on) col = 0xd0ffffff;
         else if (mk.cls == "orange") { col = 0xffff9020; wid = 2; }
         else col = 0x90808080;
         // dark outline first so white circles stay visible on bright stars
         g.pen = new Pen(0xa0000000, wid + 2);
         g.drawCircle(sp.x, sp.y, rad);
         g.pen = new Pen(col, wid);
         g.drawCircle(sp.x, sp.y, rad);
      }
      g.end();
   };

   this.setZoom = function (nz, sx, sy) {      // zoom around screen point (sx, sy)
      nz = Math.min(Math.max(nz, 0.25), self.zFit);
      var X = self.vx + sx * self.zoom, Y = self.vy + sy * self.zoom;
      self.zoom = nz;
      self.vx = X - sx * nz; self.vy = Y - sy * nz;
      if (nz >= self.zFit - 1e-9) { self.vx = (W - self.cw * nz) / 2; self.vy = (H - self.ch * nz) / 2; }
      self.clampView();
      self.renderView();
      self.updateZoomLabel();
      self.canvas.update();
   };
   this.canvas.onMouseWheel = function (x, y, delta, buttonState, modifiers) {
      self.setZoom(delta > 0 ? self.zoom / 1.5 : self.zoom * 1.5, x, y);
   };

   this.clickAt = function (x, y) {
      var rad = self.markRadius(), best = -1, bd = 1e9;
      for (var i = 0; i < self.marks.length; ++i) {
         var sp = self.toScreen(self.marks[i].x, self.marks[i].y), dx = sp.x - x, dy = sp.y - y, d = Math.sqrt(dx * dx + dy * dy);
         if (d < bd) { bd = d; best = i; }
      }
      if (best >= 0 && bd <= rad + 4) {
         var mk = self.marks[best];
         if (mk.cls == "manual") {
            for (var k = 0; k < self.extra.length; ++k)
               if (self.extra[k].x == mk.x && self.extra[k].y == mk.y) { self.extra.splice(k, 1); break; }
            self.marks.splice(best, 1);
         } else {
            mk.on = !mk.on;
            self.override[GB_key(mk.x, mk.y)] = mk.on;
         }
      } else {
         var X = self.vx + (x + 0.5) * self.zoom, Y = self.vy + (y + 0.5) * self.zoom;
         if (X < 0 || Y < 0 || X >= W || Y >= H) return;
         var pt = GB_snap(self.chans, W, H, X, Y, Math.round(12 * self.auto.scale) + Math.ceil(self.zoom));
         self.extra.push({ x: pt.x, y: pt.y });
         self.marks.push({ x: pt.x, y: pt.y, cls: "manual", on: true, status: "manual", conf: 0 });
      }
      self.updateCounts();
      self.canvas.update();
   };
   this.canvas.onMousePress = function (x, y, button, buttonState, modifiers) {
      self.drag = { x: x, y: y, vx: self.vx, vy: self.vy, moved: false };
   };
   this.canvas.onMouseMove = function (x, y, buttonState, modifiers) {
      if (!self.drag) return;
      var dx = x - self.drag.x, dy = y - self.drag.y;
      if (!self.drag.moved && Math.sqrt(dx * dx + dy * dy) < 4) return;
      self.drag.moved = true;
      self.dragDx = dx; self.dragDy = dy;          // live: shift the current bitmap
      self.canvas.update();
   };
   this.canvas.onMouseRelease = function (x, y, button, buttonState, modifiers) {
      if (!self.drag) return;
      var d = self.drag;
      self.drag = null;
      if (!d.moved) { self.clickAt(x, y); return; }
      self.vx = d.vx - (x - d.x) * self.zoom;
      self.vy = d.vy - (y - d.y) * self.zoom;
      self.dragDx = 0; self.dragDy = 0;
      self.clampView();
      self.renderView();
      self.canvas.update();
   };

   this.zoomInButton = new PushButton(this);
   this.zoomInButton.text = "+";
   this.zoomInButton.toolTip = "Zoom in";
   this.zoomInButton.onClick = function () { self.setZoom(self.zoom / 1.5, self.cw / 2, self.ch / 2); };
   this.zoomOutButton = new PushButton(this);
   this.zoomOutButton.text = "-";
   this.zoomOutButton.toolTip = "Zoom out";
   this.zoomOutButton.onClick = function () { self.setZoom(self.zoom * 1.5, self.cw / 2, self.ch / 2); };
   this.fitButton = new PushButton(this);
   this.fitButton.text = "Fit";
   this.fitButton.toolTip = "Show the whole image";
   this.fitButton.onClick = function () { self.setZoom(self.zFit, self.cw / 2, self.ch / 2); };
   this.zoomLabel = new Label(this);
   this.updateZoomLabel = function () {
      self.zoomLabel.text = self.zoom >= 1 ? "Zoom 1:" + self.zoom.toFixed(1) : "Zoom " + (1 / self.zoom).toFixed(1) + ":1";
   };
   this.zoomBar = new HorizontalSizer;
   this.zoomBar.spacing = 6;
   this.zoomBar.add(this.zoomOutButton);
   this.zoomBar.add(this.zoomInButton);
   this.zoomBar.add(this.fitButton);
   this.zoomBar.add(this.zoomLabel);
   this.zoomBar.addStretch();
   this.left = new VerticalSizer;
   this.left.spacing = 6;
   this.left.add(this.canvas);
   this.left.add(this.zoomBar);
   this.renderView();
   this.updateZoomLabel();

   this.analyse = function () {
      console.writeln("<end><cbr>Analysing stars (nothing is modified yet)...");
      processEvents();
      var t0 = Date.now();
      var res = GB.process(self.chans, W, H, self.effectiveParams(), null, function () {}, { dryRun: true, progress: GB_progress() });
      var p = self.effectiveParams(), marks = [];
      for (var i = 0; i < res.stars.length; ++i) {
         var s = res.stars[i];
         if (s.manual) continue;
         var cls = (s.status == "ok" && s.conf >= p.confLo) ? "white" : "orange";
         var key = GB_key(s.x, s.y);
         var on = (key in self.override) ? self.override[key] : (cls == "white");
         marks.push({ x: s.x, y: s.y, cls: cls, on: on, status: s.status, conf: s.conf || 0 });
      }
      for (var k = 0; k < self.extra.length; ++k)
         marks.push({ x: self.extra[k].x, y: self.extra[k].y, cls: "manual", on: true, status: "manual", conf: 0 });
      self.marks = marks;
      console.writeln(format("Analysis done in %.1f s: %d candidate stars.", (Date.now() - t0) / 1000, marks.length));
   };

   // ---------- right column ----------
   this.helpLabel = new Label(this);
   this.helpLabel.wordWrapping = true;
   this.helpLabel.useRichText = true;
   this.helpLabel.minWidth = 340;
   this.helpLabel.maxWidth = 380;
   this.helpLabel.text = "<p><b>" + TITLE + " v" + VERSION + "</b></p>" +
      "<p>Removes reflection ghost rings around bright saturated stars. Run it on <b>linear</b> data, " +
      "before BlurXTerminator and StarXTerminator. The original image is not modified.</p>" +
      "<p>Target: <b>" + GB_esc(win.mainView.id) + "</b> (" + W + " x " + H + ")</p>";

   this.legendLabel = new Label(this);
   this.legendLabel.wordWrapping = true;
   this.legendLabel.useRichText = true;
   this.legendLabel.maxWidth = 380;
   this.legendLabel.text =
      "<p><b>White</b> circle: will be corrected.<br>" +
      "<b><span style='color:#e07010'>Orange</span></b>: uncertain (nebula or weak ghost), NOT corrected until you click it.<br>" +
      "<b><span style='color:#808080'>Grey</span></b>: excluded.<br>" +
      "<b><span style='color:#20a040'>Green</span></b>: star you added.</p>" +
      "<p>Click a circle to switch it on/off. Click a star without a circle to add it; click a green circle to remove it.</p>" +
      "<p>Mouse wheel or +/- to zoom, drag to move around, Fit to see the whole image.</p>";

   this.countLabel = new Label(this);
   this.countLabel.useRichText = true;
   this.updateCounts = function () {
      var nOn = 0, nOr = 0, nMan = 0, nOff = 0;
      for (var i = 0; i < self.marks.length; ++i) {
         var mk = self.marks[i];
         if (mk.cls == "manual") ++nMan;
         else if (mk.on) ++nOn;
         else if (mk.cls == "orange") ++nOr;
         else ++nOff;
      }
      self.countLabel.text = "<b>" + (nOn + nMan) + "</b> star(s) will be corrected (" + nMan + " added). " +
         nOr + " uncertain, " + nOff + " excluded.";
   };

   this.autoLabel = new Label(this);
   this.autoLabel.wordWrapping = true;
   this.autoLabel.useRichText = true;
   this.autoLabel.maxWidth = 380;
   this.autoLabel.text = "<p><b>Automated mode</b>: image scale " + this.auto.scale.toFixed(2) +
      " (from " + GB_esc(this.auto.source) + ")." +
      (this.auto.fromMetadata ? "" : " <b>Estimated:</b> no pixel scale in the metadata; check the scale in Manual mode " +
                                     "if the image was drizzled, binned or resampled.") + "</p>";

   this.reanalyseButton = new PushButton(this);
   this.reanalyseButton.text = "Re-analyse";
   this.reanalyseButton.toolTip = "<p>Run the star analysis again (needed after changing manual settings). Your clicks are kept.</p>";
   this.reanalyseButton.onClick = function () {
      try { self.analyse(); self.updateCounts(); self.canvas.update(); }
      catch (e) { GB_reportError(e); }
   };

   var labelWidth = this.font.width("Ghost outer radius (px):M");
   function numeric(label, tip, lo, hi, prec, key) {
      var nc = new NumericControl(self);
      nc.label.text = label;
      nc.label.minWidth = labelWidth;
      nc.setRange(lo, hi);
      nc.slider.setRange(0, 1000);
      nc.slider.minWidth = 120;
      nc.setPrecision(prec);
      nc.setValue(self.params[key]);
      nc.toolTip = tip;
      nc.onValueUpdated = function (value) { self.params[key] = value; };
      return nc;
   }
   this.satControl = numeric("Saturation level:", "<p>Fraction of the brightest value above which a pixel counts as saturated.</p>", 0.50, 1.00, 2, "satThreshold");
   this.minAreaControl = numeric("Min. core area (px):", "<p>Minimum saturated area for a star to be a candidate.</p>", 1, 500, 0, "minArea");
   this.scaleControl = numeric("Scale factor:", "<p>Multiplies every radius. 1.0 native Draco, 2.0 if drizzled 2x.</p>", 0.25, 4.00, 2, "scale");
   this.annInControl = numeric("Ghost inner radius (px):", "<p>Inner radius of the ghost zone (at scale 1).</p>", 4, 200, 0, "annIn");
   this.annOutControl = numeric("Ghost outer radius (px):", "<p>Outer radius of the ghost zone (at scale 1).</p>", 8, 400, 0, "annOut");
   this.rOutControl = numeric("Halo fit radius (px):", "<p>Radius used to model the stellar halo (at scale 1).</p>", 20, 600, 0, "rOut");
   this.gateControl = numeric("Nebula safety gate:", "<p>Stars in too structured surroundings are marked uncertain.</p>", 2, 100, 1, "gateK");
   this.lopControl = numeric("Ghost lopsidedness:", "<p>How much brighter the ghost may be on one side (0 = round).</p>", 0, 1, 2, "maxLop");
   this.scaleControl.setValue(this.auto.scale);

   this.manualCheck = new CheckBox(this);
   this.manualCheck.text = "Manual mode (show advanced settings)";
   this.manualCheck.checked = false;
   this.manualCheck.onCheck = function (checked) {
      self.params.manual = checked;
      self.paramBox.visible = checked;
      self.autoLabel.visible = !checked;
      self.adjustToContents();
   };

   this.ghostCheck = new CheckBox(this);
   this.ghostCheck.text = "Create ghost model image";
   this.ghostCheck.checked = this.params.makeGhost;
   this.ghostCheck.onCheck = function (checked) { self.params.makeGhost = checked; };

   this.paramBox = new GroupBox(this);
   this.paramBox.title = "Manual settings";
   this.paramBox.sizer = new VerticalSizer;
   this.paramBox.sizer.margin = 6;
   this.paramBox.sizer.spacing = 4;
   var ctl = [this.satControl, this.minAreaControl, this.scaleControl, this.annInControl,
              this.annOutControl, this.rOutControl, this.gateControl, this.lopControl];
   for (var i = 0; i < ctl.length; ++i) this.paramBox.sizer.add(ctl[i]);
   this.paramBox.visible = false;

   this.okButton = new PushButton(this);
   this.okButton.text = "Execute";
   this.okButton.onClick = function () {
      if (self.params.manual && self.params.annOut <= self.params.annIn + 4) {
         (new MessageBox("Ghost outer radius must be larger than inner radius (+4 px).", TITLE, StdIcon_Error, StdButton_Ok)).execute();
         return;
      }
      if (self.params.manual && self.params.rOut <= self.params.annOut + 10) {
         (new MessageBox("Halo fit radius must exceed the ghost outer radius by at least 10 px.", TITLE, StdIcon_Error, StdButton_Ok)).execute();
         return;
      }
      self.ok();
   };
   this.cancelButton = new PushButton(this);
   this.cancelButton.text = "Cancel";
   this.cancelButton.onClick = function () { self.cancel(); };

   this.buttons = new HorizontalSizer;
   this.buttons.spacing = 6;
   this.buttons.add(this.reanalyseButton);
   this.buttons.addStretch();
   this.buttons.add(this.okButton);
   this.buttons.add(this.cancelButton);

   this.right = new VerticalSizer;
   this.right.spacing = 8;
   this.right.add(this.helpLabel);
   this.right.add(this.legendLabel);
   this.right.add(this.countLabel);
   this.right.add(this.autoLabel);
   this.right.add(this.ghostCheck);
   this.right.add(this.manualCheck);
   this.right.add(this.paramBox);
   this.right.addStretch();
   this.right.add(this.buttons);

   this.sizer = new HorizontalSizer;
   this.sizer.margin = 8;
   this.sizer.spacing = 10;
   this.sizer.add(this.left);
   this.sizer.add(this.right);

   this.windowTitle = TITLE + " v" + VERSION;

   this.analyse();
   this.updateCounts();
   this.adjustToContents();
}
GB_Dialog.prototype = new Dialog;

// Progress callback for GB.process: writes the stage and count at most every 1.5 s, keeps
// PixInsight responsive and lets the console Abort button stop the run.
function GB_progress() {
   var last = "", lastT = 0;
   return function (stage, i, n) {
      var now = Date.now();
      if (stage != last || now - lastT > 1500 || i == n - 1) {
         console.writeln(format("  %s: %d / %d", stage, i + 1, n));
         last = stage; lastT = now;
      }
      processEvents();
      try { return !!console.abortRequested; } catch (e) { return false; }
   };
}

// Environment summary written at start, to be pasted in bug reports.
// The platform is inferred from the form of the home directory, which is not written out.
function GB_environment() {
   var ver = "unknown", plat = "unknown platform";
   try { ver = String(typeof GB_PI_VERSION != "undefined" ? GB_PI_VERSION : "unknown"); } catch (e) { /* keep "unknown" */ }
   try {
      var home = String(File.homeDirectory);
      plat = /^[A-Za-z]:/.test(home) ? "Windows" : (home.indexOf("/Users/") == 0 ? "macOS" : "Linux/Unix");
   } catch (e) { /* keep "unknown platform" */ }
   return "PixInsight " + ver + ", " + plat;
}

function GB_reportError(e) {
   var msg = e && e.message ? e.message : String(e);
   if (msg == "Aborted by the user.") { console.warningln("<end><cbr>" + TITLE + ": aborted by the user."); return; }
   console.criticalln("<end><cbr>" + TITLE + " v" + VERSION + " stopped: " + msg +
                      (e && e.lineNumber ? " (line " + e.lineNumber + ")" : ""));
   if (e && e.stack) console.writeln(String(e.stack));
   console.writeln("If this looks like a bug, please open an issue at github.com/Jiness/GhostBuster/issues " +
                   "and paste the whole console output of this run.");
}

function main() {
   console.show();
   console.abortEnabled = true;
   try { GB_run(); }
   catch (e) { GB_reportError(e); }
   finally { console.abortEnabled = false; }
}

function GB_run() {
   console.show();
   var win = ImageWindow.activeWindow;
   if (win.isNull) {
      (new MessageBox("Open the linear image to process first.", TITLE, StdIcon_Error, StdButton_Ok)).execute();
      return;
   }
   var img = win.mainView.image, W = img.width, H = img.height;
   console.writeln("<end><cbr><br><b>" + TITLE + " v" + VERSION + "</b>");
   console.writeln("Target: " + win.mainView.id);
   console.writeln(GB_environment() + format("; image %d x %d, %d channel(s)", W, H, img.numberOfChannels));
   var auto0 = GB.autoScale(GB_keywords(win), W, H);
   console.writeln(format("Pixel scale factor %.2f (from %s).", auto0.scale, auto0.source));
   // peak memory: per pixel and channel, 4 float32 buffers (input copy, working copy, ghost, output)
   // and a half-resolution preview level; per pixel, detection buffers (float32 luminance,
   // two int32 arrays, three byte masks)
   var nch = img.numberOfChannels, estMB = Math.round(W * H * (nch * (16 + 1) + 4 + 8 + 3) / 1048576);
   if (estMB > 16384) {
      (new MessageBox("This image is too large for this script (about " + estMB + " MB of memory needed).",
                      TITLE, StdIcon_Error, StdButton_Ok)).execute();
      return;
   }
   if (estMB > 4096) {
      var mbm = new MessageBox("This image is large (" + W + " x " + H + "): about " + estMB +
         " MB of memory will be needed.\nContinue?", TITLE, StdIcon_Warning, StdButton_Yes, StdButton_No);
      if (mbm.execute() != StdButton_Yes) return;
   }
   var chans = GB_readChannels(img);
   // Linearity check on a sample spread over the whole image. Heuristic only: a linear image
   // has most pixels near the background level, so a median above 0.25 suggests a stretched
   // image. Very bright fields or images with a large pedestal can trigger it wrongly, and a
   // lightly stretched image can pass, hence a warning rather than an error.
   var sample = [], step = Math.max(1, Math.floor(W * H / 500000));
   for (var i = 0; i < W * H; i += step) sample.push(chans[0][i]);
   var med = GB.median(sample, sample.length);
   console.writeln(format("Image median %.4f (linear data is usually below 0.05).", med));
   if (med > 0.25) {
      var mb = new MessageBox("The image median is " + med.toFixed(3) + ". It looks stretched.\n" +
         "This script is designed for LINEAR data. Continue anyway?", TITLE, StdIcon_Warning, StdButton_Yes, StdButton_No);
      if (mb.execute() != StdButton_Yes) return;
   }

   var params = GB.defaults();
   params.makeGhost = true;
   var dlg = new GB_Dialog(win, params, chans);
   if (!dlg.execute()) return;

   var run = dlg.effectiveParams();
   run.makeGhost = dlg.params.makeGhost;
   console.writeln(dlg.params.manual ? "Manual mode." : "Automated mode: scale " + run.scale.toFixed(2) + " (from " + dlg.auto.source + ")");
   if (!dlg.params.manual && !dlg.auto.fromMetadata)
      console.warningln("No pixel scale in the image metadata: the scale was estimated from the image size. " +
                        "If the image was drizzled, binned or resampled, set the scale in Manual mode.");

   var sel = {}, forced = {};
   for (var i = 0; i < dlg.marks.length; ++i) {
      var mk = dlg.marks[i];
      if (mk.cls == "manual") continue;
      sel[GB_key(mk.x, mk.y)] = mk.on;
      if (mk.on && mk.status == "gate") forced[GB_key(mk.x, mk.y)] = true;
   }
   var extra = [];
   for (var k = 0; k < dlg.extra.length; ++k) extra.push({ x: dlg.extra[k].x, y: dlg.extra[k].y, snapped: true });
   var selectFn = function (s) { return s.manual ? true : !!sel[GB_key(s.x, s.y)]; };
   var forceFn = function (s) { return !!forced[GB_key(s.x, s.y)]; };

   var t0 = Date.now();
   var res = GB.process(chans, W, H, run, selectFn, function (msg) { console.writeln(msg); processEvents(); },
                         { extra: extra, forceFn: forceFn, progress: GB_progress() });

   var nc = chans.length, out = [];
   for (var c = 0; c < nc; ++c) {
      var a = new Float32Array(W * H), src = chans[c], g = res.ghost[c];
      for (var i = 0; i < W * H; ++i) { var v = src[i] - g[i]; a[i] = v < 0 ? 0 : (v > 1 ? 1 : v); }
      out.push(a);
   }
   var w1 = GB_makeWindow(win, win.mainView.id + "_deghosted", out, W, H);
   w1.show();
   if (run.makeGhost) {
      var w2 = GB_makeWindow(win, win.mainView.id + "_ghost", res.ghost, W, H);
      w2.show();
   }
   console.writeln("Ghost shape: " + res.templateSource);
   console.writeln(format("<br>Done in %.1f s: %d star(s) corrected.", (Date.now() - t0) / 1000, res.done));
}

main();
