// Runs the complete script (dialog included) against a minimal simulation of the PixInsight
// scripting API, on a synthetic image, and simulates user actions in the star selection preview.
// Usage: node tests/ui_test.js
"use strict";
const assert = require("assert");
const { source } = require("./load.js");
const { makeField } = require("./synth.js");

const f = makeField({ seed: 11, bright: 6, W: 1200, H: 900 });
const { W, H } = f;
const srcWrites = [], logs = [], outputs = [];
const print = (s) => process.stdout.write(s + "\n");
const g = global;

// ---- minimal PJSR simulation ----
g.console = { writeln: (m) => logs.push(String(m)), warningln: (m) => logs.push("WARNING " + m), show: () => {} };
g.processEvents = () => {};
g.format = (fmt, ...a) => { let i = 0; return fmt.replace(/%(\.\d+)?[fd]/g, (m, p) => { const v = a[i++]; return p ? Number(v).toFixed(+p.slice(1)) : String(Math.round(v)); }); };
for (const k of ["StdIcon_Error", "StdIcon_Warning", "StdButton_Ok", "StdButton_Yes", "StdButton_No", "ColorSpace_RGB", "SampleType_Real", "UndoFlag_NoSwapFile"]) g[k] = k;
g.Rect = function (x0, y0, x1, y1) { Object.assign(this, { x0, y0, x1, y1 }); };
g.Pen = function (c, w) { this.color = c; this.width = w; };
g.Brush = function (c) { this.color = c; };
g.Bitmap = function (w, h) { this.width = w; this.height = h; this.setPixel = () => {}; };
g.Image = function (w, h) { this.setSamples = (a) => { assert.strictEqual(a.length, w * h); }; this.render = () => new g.Bitmap(w, h); this.free = () => {}; };
let circles = 0;
g.Graphics = function () { this.drawBitmap = () => {}; this.fillRect = () => {}; this.drawCircle = () => { ++circles; }; this.end = () => {}; };
function Widget() { this.visible = true; }
for (const k of ["Label", "PushButton", "CheckBox", "GroupBox"]) g[k] = function () { Widget.call(this); };
g.Control = function () { Widget.call(this); this.setFixedSize = () => {}; this.update = () => {}; };
g.NumericControl = function () { Widget.call(this); this.label = {}; this.slider = { setRange: () => {} }; this.setRange = () => {}; this.setPrecision = () => {}; this.setValue = () => {}; };
g.VerticalSizer = g.HorizontalSizer = function () { this.add = (x) => assert.ok(x !== undefined, "sizer.add(undefined)"); this.addStretch = () => {}; };
g.MessageBox = function () { this.execute = () => g.StdButton_Yes; };
let simulateUser = null;
g.Dialog = function () { this.font = { width: (s) => 7 * s.length }; };
g.Dialog.prototype.adjustToContents = function () {};
g.Dialog.prototype.ok = function () {};
g.Dialog.prototype.cancel = function () {};
g.Dialog.prototype.execute = function () { simulateUser(this); return true; };
g.ImageWindow = function (w, h, nc, bits, isFloat, isColor, id) {
   const data = [];
   this.mainView = { id, beginProcess: () => {}, endProcess: () => {}, image: { setSamples: (a, r, c) => { data[c] = Float32Array.from(a); } } };
   this.data = data; this.show = () => {}; this.copyAstrometricSolution = () => {};
   outputs.push(this);
};
g.ImageWindow.activeWindow = {
   isNull: false,
   keywords: [{ name: "XPIXSZ", value: "2.394" }, { name: "FOCALLEN", value: "340." }],
   mainView: { id: "Synthetic_Field", image: {
      width: W, height: H, numberOfChannels: 3,
      getSamples: (a, r, c) => a.set(f.chans[c]),
      setSamples: () => srcWrites.push("source image written") } }
};

// ---- load and run ----
const { code, version, title } = source();
g.VERSION = version; g.TITLE = title;
const before = f.chans.map(a => a.slice());
let state = {};
simulateUser = (dlg) => {
   dlg.canvas.onPaint();
   state.circlesFirstPaint = circles;
   const marks = dlg.marks.slice();
   const white = marks.find(m => m.cls == "white"), orange = marks.find(m => m.cls == "orange");
   assert.ok(white, "at least one star marked for correction");
   const scr = (m) => dlg.toScreen(m.x, m.y);
   const click = (p) => { dlg.canvas.onMousePress(p.x, p.y, 1, 0, 0); dlg.canvas.onMouseRelease(p.x, p.y, 1, 0, 0); };
   // zoom in three notches on the first white star, then exclude it with a click
   for (let k = 0; k < 3; ++k) { const p = scr(white); dlg.canvas.onMouseWheel(p.x, p.y, 120, 0, 0); }
   state.zoomed = dlg.zoom < dlg.zFit;
   click(scr(white));
   state.whiteExcluded = { x: white.x, y: white.y, on: white.on };
   // a drag must pan the view without toggling anything
   const vx0 = dlg.vx, on0 = white.on;
   dlg.canvas.onMousePress(300, 300, 1, 0, 0); dlg.canvas.onMouseMove(360, 330, 1, 0); dlg.canvas.onMouseRelease(400, 350, 1, 0, 0);
   state.panned = Math.abs(dlg.vx - (vx0 - 100 * dlg.zoom)) < 1e-6 && white.on === on0;
   // confirm an orange star if there is one
   dlg.fitButton.onClick();
   if (orange) { click(scr(orange)); state.orangeConfirmed = { x: orange.x, y: orange.y, on: orange.on }; }
   // add a star on a faint field star
   const faint = f.stars.find(s => !s.ghost && s.peak > 0.15 && s.x > 60 && s.y > 60 && s.x < W - 60 && s.y < H - 60 &&
                                   marks.every(m => Math.hypot(m.x - s.x, m.y - s.y) > 120));
   dlg.vx = faint.x - dlg.cw * dlg.zoom / 2; dlg.vy = faint.y - dlg.ch * dlg.zoom / 2; dlg.clampView(); dlg.renderView();
   click(dlg.toScreen(faint.x, faint.y));
   state.added = dlg.marks.filter(m => m.cls == "manual").length;
};
new Function(code)();

// ---- checks ----
let failed = 0;
function check(name, cond) { print((cond ? "  ok    " : "  FAIL  ") + name); if (!cond) ++failed; }
const ran = (x, y) => logs.some(l => l.startsWith("OK   (" + Math.round(x) + ", " + Math.round(y) + ")"));
check("candidates drawn in the preview", state.circlesFirstPaint > 0);
check("mouse wheel zooms in", state.zoomed);
check("click excludes a white star", state.whiteExcluded.on === false);
check("drag pans without toggling", state.panned);
check("excluded star is not processed", !ran(state.whiteExcluded.x, state.whiteExcluded.y));
if (state.orangeConfirmed) check("confirmed orange star is processed", state.orangeConfirmed.on && ran(state.orangeConfirmed.x, state.orangeConfirmed.y));
check("manually added star is listed", state.added === 1);
check("two output windows created", outputs.map(o => o.mainView.id).join(",") === "Synthetic_Field_deghosted,Synthetic_Field_ghost");
check("source image never written", srcWrites.length === 0);
check("source data unchanged", f.chans.every((a, c) => a.every((v, i) => v === before[c][i])));
const out = outputs[0].data, gh = outputs[1].data;
let consistent = true;
for (let c = 0; c < 3 && consistent; ++c) for (let i = 0; i < W * H; i += 37) {
   const expect = Math.min(1, Math.max(0, f.chans[c][i] - gh[c][i]));
   if (Math.abs(out[c][i] - expect) > 1e-6) { consistent = false; break; }
}
check("deghosted = source - ghost (clipped to [0, 1])", consistent);
check("script reported completion", logs.some(l => /Done in/.test(l)));
print("\n" + (failed ? failed + " check(s) failed" : "all checks passed"));
process.exitCode = failed ? 1 : 0;
