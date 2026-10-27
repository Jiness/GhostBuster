# Changelog

## 0.1.3

- Author credit in the dialog.
- CI: tests also run on macOS with Apple Silicon (M1), Node.js 24 instead of 20 (end of life),
  GitHub actions at their current major versions.
- Dependabot: minor and patch updates of the GitHub actions are merged automatically once the
  required checks pass; major updates stay for manual review.

## 0.1.2

Diagnostics, after a report of PixInsight appearing frozen on a Mac M1:

- Errors are caught and written to the Process Console (message, line, call stack) with a link
  to the issue tracker, instead of leaving the script stopped silently.
- At start, the console shows the PixInsight version and platform, the image size and channels,
  the pixel scale factor with its source, and the image median: everything needed in a bug report.
- Progress is written during analysis and correction (stage and star count, at most every
  1.5 s); PixInsight stays responsive between stars and the console Abort button stops the run.
- At most 200 candidate stars are analysed (largest saturated cores first, `maxStars`), which
  bounds the analysis time on images with very many saturated areas.
- No more "reference to undefined property" warning when the image has no astrometric solution.
- Processing results are unchanged (checked pixel by pixel on two stacks). 39 tests.

## 0.1.1

Fixes found on IC 1805 (Heart Nebula):

- **Dark hole next to a star in bright red nebula.** When the fitted ghost exceeded the physical
  ceiling (ratio cap), it was only scaled down and kept its free shape and colour, so it still
  followed the nebula. A fit above the ceiling is now discarded: the star gets the typical ghost
  for its brightness, as a round ring, without free correction.
- **Green ring left around some stars.** The ghost colour came from a colour ratio learned on the
  image, which picks up the colour of the nebula (red x1.3 to x1.6 depending on the image). The
  colour split is now chosen so that the corrected disk keeps the colour of the star light just
  outside the core, within the same bounds around the predicted colour. Where the background
  varies strongly around the star, the previous per-channel fit is kept.
- Tests: ghost colour per star, and a star next to a red emission patch (no sector of the
  corrected disk over-subtracted). 37 tests.
- New `bench/run_bench.js`: runs the automated mode on real FITS stacks and writes contact sheets
  and per-star figures to `bench-results/v<version>/`, for comparison between versions.

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
