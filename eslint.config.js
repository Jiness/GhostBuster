// ESLint configuration (ESLint 9, flat config) for the parsable copy produced by tests/lint_prepare.js.
// PixInsight scripts run on an ECMAScript 5 engine, so newer syntax is reported as an error.
// ES5 requires a parameter in every catch clause, hence caughtErrors: "none".
"use strict";
module.exports = [
  {
    files: ["build/GhostBuster.lint.js"],
    languageOptions: {
      ecmaVersion: 5,
      sourceType: "script",
      globals: {
      Bitmap: "readonly",
      Brush: "readonly",
      CheckBox: "readonly",
      ColorSpace_RGB: "readonly",
      Control: "readonly",
      Dialog: "readonly",
      Float32Array: "readonly",
      Float64Array: "readonly",
      Graphics: "readonly",
      GroupBox: "readonly",
      HorizontalSizer: "readonly",
      Image: "readonly",
      ImageWindow: "readonly",
      Int32Array: "readonly",
      Label: "readonly",
      MessageBox: "readonly",
      NumericControl: "readonly",
      Pen: "readonly",
      PushButton: "readonly",
      Rect: "readonly",
      SampleType_Real: "readonly",
      StdButton_No: "readonly",
      StdButton_Ok: "readonly",
      StdButton_Yes: "readonly",
      StdIcon_Error: "readonly",
      StdIcon_Warning: "readonly",
      Uint8Array: "readonly",
      UndoFlag_NoSwapFile: "readonly",
      VerticalSizer: "readonly",
      console: "readonly",
      format: "readonly",
      isFinite: "readonly",
      processEvents: "readonly"
      }
    },
    rules: {
      "no-unused-vars": ["warn", {"vars": "all", "args": "none", "caughtErrors": "none"}],
      "no-undef": "error",
      "no-redeclare": "off",
      "no-unreachable": "error",
      "no-dupe-keys": "error",
      "no-empty": "warn",
      "no-self-assign": "error",
      "no-constant-condition": "warn",
      "no-shadow": "off",
      "eqeqeq": "off"
    }
  }
];
