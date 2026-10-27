# Changelog

## 0.1.0 (first public pre-release)

Initial public version: ghost-ring removal for Dwarflab Draco duo-band stacks, interactive star
selection with zoom and pan, automated scale detection, manual mode.

### Pre-release review

| # | Finding | Resolution |
|---|---|---|
| 1 | Personal references | Copyright holder is the repository owner name only; no social-media links. |
| 2 | Repository URLs | Aligned with `Jiness/GhostBuster`. |
| 3 | `.gitignore` | Extended to credential files, keys and development artefacts. New repository: no history to clean. |
| 4 | PNG metadata | Checked: no text, EXIF or time chunks in `docs/*.png`. |
| 5 | Manual star coordinates | Non-finite, null and out-of-bounds coordinates are ignored; snapping no longer writes into caller objects. |
| 6, 13 | Parameter validation | New `GB.validate`: every setting is finite and clamped; `confHi > confLo`, `annOut > annIn`, `rOut > annOut` enforced. |
| 7 | Memory | Peak estimate covers all buffers; warning above 4 GB, refusal above 16 GB. |
| 8 | Non-finite pixels | Count reported in the console; only the script's working copy is changed. |
| 9 | Fallback ghost/star ratio | Was 0.5, inconsistent with measurements; now 0.03 (median of 0.016 to 0.064 measured on four Draco stacks), definition documented. |
| 10 | Scale from image size | Console warning and dialog note when the scale is estimated. |
| 11 | Dry-run allocations | Analysis no longer copies the image nor allocates output buffers. |
| 12 | Linearity check | Heuristic and its limits documented in the code. |
| 14, 15 | Tests | `tests/run_tests.js` (31 tests on synthetic fields with known ghosts) and `tests/ui_test.js` (full script and dialog against a simulated PixInsight API). |
| 16, 17 | Includes, dead code | Only used PixInsight headers kept; ESLint clean with zero warnings. |
| 18, 19 | Comments, personal data | Neutral technical comments; no names, e-mail addresses, local paths or credentials. |

### Defects found by the new tests

- Robust least squares clipped residuals around zero instead of their median, which could reject
  every point and keep a biased first fit. Clipping is now centred on the median residual.
- Automated scale: an image smaller than a native frame (typically a crop) was treated as binned and
  its radii halved. The image size now only detects enlarged (drizzled) images.
- Saturated regions larger than `maxArea` are no longer treated as stars.
