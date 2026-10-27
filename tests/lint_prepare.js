// Writes build/GhostBuster.lint.js: GhostBuster.js without the PixInsight preprocessor
// directives, so that standard JavaScript tools (ESLint, CodeQL) can parse it.
"use strict";
const fs = require("fs");
const path = require("path");
const { source } = require("./load.js");

const { code, version, title } = source();
const out = path.join(__dirname, "..", "build");
fs.mkdirSync(out, { recursive: true });
fs.writeFileSync(path.join(out, "GhostBuster.lint.js"),
   "var TITLE = " + JSON.stringify(title) + ", VERSION = " + JSON.stringify(version) + ";\n" + code);
console.log("build/GhostBuster.lint.js written (version " + version + ")");
