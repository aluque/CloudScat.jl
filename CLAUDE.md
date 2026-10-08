# CloudScat.jl — agent notes

> `CLAUDE.md` and `AGENTS.md` are mirrored copies. Keep them identical: when
> editing one, apply the same change to the other.

## What this repo is

CloudScat is a Julia Monte Carlo code that simulates light from an optical
source (lightning) inside or above a thundercloud, as seen by satellite-borne
observers. It models Rayleigh (molecular) and Mie (droplet) scattering, plus
optional absorption. Published in GMD: https://doi.org/10.5194/gmd-13-5549-2020

Layout:

- `src/` — the Julia package (`CloudScat.jl` is the entry point; `MieScatter/`
  is a vendored Mie solver by O. Wilkman).
- `samples/` — heavily commented example run scripts (`sample.jl` is the
  simplest). `samples/gmd/` reproduces the paper figures.
- `util/` — Python/matplotlib scripts to plot output (`plot_photo.py`,
  `plot_image.py`).
- `viewer/` — a static browser app for viewing output files (see below).
- `data/` — refractive index tables used by the Mie code.

Running a simulation: `julia --color=yes samples/sample.jl` (or
`include("sample.jl"); run()` from the REPL). It writes an `.h5` file named
after the script.

## Output file format (the viewer's input)

Written by `save()` in `src/CloudScat.jl`. Everything is SI (m, s, m⁻³, …).

- `/parameters` — a group holding **attributes only**, one for each field of
  `Params` (`src/parameters.jl`): `zmin`, `source_a`, `source_b`, `N`, `λ`,
  `weight_min`, `σray`, `max_iter`, `nair`, `νray_ground`, `H`, `νraymax`.
  The `source` field is expanded as `source.<field>`. The field names depend on
  the source type, e.g. `source.a`/`source.b` for segment sources,
  `source.r0` for point sources. Attribute names use Unicode (`λ`, `σray`, `ν…`).
- `/obsNNNNN` (`obs00001`, `obs00002`, … 1-based, zero-padded to 5 digits), one
  group per observer:
  - attributes: `altitude` (= r[3], m), `shift` (= r[1], m), `delay` (s,
    light travel time from the source centroid).
  - datasets: `t` (bin centres, s), `timeline` (photons / m² / s per source
    photon), `image` (2-D, pixels × pixels), `depth` (2-D, same shape as
    `image`). All are deflate-compressed with shuffle.
- `/cloudgrid` (optional, only when the `savecloudgrid` callback is used):
  `density` (3-D), `x`, `y`, `z`.

Julia is column-major and HDF5/JS are row-major, so `image` is **transposed**
relative to the Julia array (`util/plot_image.py` notes this too).

If you change `save()` or `Params`, check that the viewer still works. Its
`PARAM_FORMAT` table in `viewer/app.js` is keyed by `Params` field names.

## The viewer (`viewer/`)

A single-page, browser-only HDF5 viewer. It is published via GitHub Pages at
https://aluque.github.io/CloudScat.jl/viewer/ and served directly from the
`viewer/` folder on `master`. There is no build step, no bundler, no
`package.json`, and no tests.

- `viewer/index.html` holds the markup plus all CSS in one inline `<style>`.
  It loads its dependencies from the jsDelivr CDN:
  - `pako@2.2.0` (zlib, needed by jsfive for deflate; jsfive requires
    `pako ^2`, and pako 3.x has no `dist/pako.min.js`)
  - `jsfive@0.4.2` (pure-JS HDF5 reader; `hdf5.File`)
  - `plotly.js-dist@2.27.1`
- `viewer/app.js` holds all the logic in one IIFE, written in **ES5 style**
  (`var`, `function`, no modules, no arrow functions). Keep that style for
  consistency.

### How app.js is organized

Sections are separated by `// ── Name ───` banner comments:

- **State**: a single `state` object (`file`, `fname`, `obsKeys`,
  `currentObs`, `paramAttrs`, `imageLog`, `imageData`, `cmap`,
  `cmapReverse`, `cmin`, `cmax`, `autoLimits`, `lightcurveLog`,
  `currentTab`, `theme`).
- **Theme colours**: `THEMES.dark` / `THEMES.light` hold the Plotly colours.
  `tc()` returns the active set. The CSS side uses custom properties on
  `:root`, overridden by `body.theme-light`. A theme change must update
  **both**: the CSS vars and the `THEMES` entries used in Plotly layouts.
- **Colormaps**: `COLORMAPS` is a list of `{ name, stops }` entries, with
  stops given as `[t, '#hex']` (the matplotlib maps are sampled at 10
  points). `colormapStops(name, reverse)` feeds both the Plotly
  `colorscale` and the CSS gradient swatches (`cssGradient`).
- **File loading**: `FileReader` reads the file into an ArrayBuffer, then
  `new hdf5.File(ab, name)` and `parseFile()`. The file must contain at least
  one `obsNNNNN` group. Loading happens via the open buttons or drag-and-drop
  anywhere on the document.
- **Sidebar**: an observer `<select>` (shown only when there are 2 or more
  observers), an observer-attributes table, and a simulation-parameters table.
  `renderParamTable(attrs, spec)` first renders the keys in the spec, in spec
  order, then every other attribute with `formatAuto`. jsfive attrs may not be
  plain objects, so the code enumerates keys defensively and uses
  `scalarVal()` to unwrap 1-element arrays.
- **Format helpers**: `PARAM_FORMAT` and `OBS_FORMAT` are
  `{ key, label, convert, unit, decimals }` specs; `decimals: null` means
  `formatAuto`. Some `PARAM_FORMAT` keys (`radius`, `nscat`, `g`, `Qext`,
  `ω₀`) are not written by the current `save()`. They are harmless, and are
  skipped when absent.
- **Tabs**: "Camera Image" (`image`) and "Light Curve" (`lc`). The image
  tab shows the "Color scale" sidebar section (`#cscale-section`). The light
  curve tab shows a log-scale checkbox in the tab bar.
  `switchTab` calls `Plotly.Plots.resize` because hidden panes have zero size.
- **Image rendering**: `renderImage(preserveRange)` reads
  `obsNNNNN/image`, turns the flat row-major array into nested rows, and
  draws a `heatmap` with the y-axis reversed. In log mode it plots log10.
  With `preserveRange`, it keeps the current zoom (read from `_fullLayout`).
  - Pixels are always square, but the axes are **not** linked with
    `scaleanchor`. With linked axes, Plotly forces the box-zoom rectangle
    to the current aspect ratio while dragging, and `constrain` does not
    change that. Instead, `squareDomains(xspan, yspan)` computes centred
    axis `domain`s that give square pixels for the given ranges.
    `renderImage` uses it for the first draw, and `fitImageArea` (on
    `plotly_afterplot`) re-applies it after every zoom, reset and resize.
    Any box shape is shown exactly, with blank space around it.
  - `fitImageArea` also moves the colorbar to the image's right edge and
    matches its height (minimum length 0.5). Its tolerance checks stop the
    relayout/restyle from re-triggering itself.
- **Colour scale controls**: the Linear/Log segmented control, Min/Max
  inputs, *Fit to view*, *Auto*, the colormap grid, and *Reverse*.
  - `state.cmin`/`cmax` hold the limits in **physical units** in both
    scales (null means auto). `colorLimits()` merges them with the auto
    limits (full range when linear, 4 decades below the max when log),
    flags unusable inputs with `.invalid` plus a tooltip, and shows the
    auto values as placeholders.
  - `fitLimitsToView()` uses the pixels inside the current axis ranges. In
    log mode it keeps the min within 4 decades of the max, because MC
    images contain stray pixels around 1e-305.
  - The colour settings persist across observer switches.
- **Light curve rendering**: `renderLightCurve(preserveRange)` plots
  `timeline` against `(t - delay)` in ms. A log y-axis spans 4 decades
  below the max.
- **Error toast**: `showError(msg)` shows a toast that disappears after
  5 seconds.
- Re-renders use `Plotly.react`. Theme and colour-scale changes re-render
  with `preserveRange = true`.

Not yet visualized: the `depth` dataset and `/cloudgrid`.
`util/plot_image.py --field depth` is the reference for plotting depth.

### Developing / testing the viewer

- Open `viewer/index.html` directly in a browser, or serve the repo root
  (e.g. `python3 -m http.server` and go to `/viewer/`). Then drop in an `.h5`
  produced by a sample run.
- Failures are logged to the console with a `[CloudScat]` prefix
  (`console.warn`/`console.error`). Don't leave debug `console.log` calls
  in committed code.
- For headless checks, Chrome can be driven over the DevTools protocol from
  Node (which has a built-in `WebSocket`). Load a file by dispatching a
  `drop` event whose `DataTransfer` holds a `File`. Test zooms with real
  `Input.dispatchMouseEvent` drags, not `Plotly.relayout`, because a
  relayout bypasses the drag-time zoom-box behaviour. Mouse coordinates
  are page coordinates: add the plot div's `getBoundingClientRect()` to
  the axis `_offset`s.
- No `.h5` files are committed. Generate one with `samples/sample.jl`, which
  can be made quick by lowering `N`.
- Keep the viewer free of build steps and server-side code: it must keep
  working as static files on GitHub Pages. If you add a CDN dependency, pin
  its version.

- Use the CSS variables for colours, never hard-coded hex values, so that both
  themes work (e.g. the `.param-empty` class for placeholder table rows).

## Conventions

- The git main branch is `master`.
- Julia ≥ 1.11 (`Project.toml`). Python helpers use numpy/matplotlib/h5py.
- Keep edits in the style of the surrounding code. Commit only when asked.
