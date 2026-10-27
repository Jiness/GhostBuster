// Loads GhostBuster.js in Node.js: removes the PixInsight preprocessor directives and
// returns either the processing core only (no PixInsight dependency) or the full script.
"use strict";
const fs = require("fs");
const path = require("path");

function source() {
   const src = fs.readFileSync(path.join(__dirname, "..", "GhostBuster.js"), "utf8");
   const version = (src.match(/#define VERSION "([^"]+)"/) || [])[1];
   const title = (src.match(/#define TITLE\s+"([^"]+)"/) || [])[1];
   const lines = src.split("\n");
   const out = [];
   let continuation = false;
   for (const l of lines) {
      if (continuation) { continuation = /\\\s*$/.test(l); continue; }
      if (l.startsWith("#")) { continuation = /\\\s*$/.test(l); continue; }
      out.push(l);
   }
   return { code: out.join("\n"), version, title };
}

function loadCore() {
   const { code } = source();
   const end = code.indexOf("// ===================== PixInsight interface");
   if (end < 0) throw new Error("interface marker not found");
   return new Function(code.slice(0, end) + "\nreturn GB;")();
}

module.exports = { source, loadCore };
