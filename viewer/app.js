// CloudScat Viewer — app.js
// Browser-only, no build step. Requires jsfive, pako, and Plotly via CDN.

(function () {
  'use strict';

  // ── State ────────────────────────────────────────────────────────────────────

  var state = {
    file: null,
    fname: '',
    obsKeys: [],
    currentObs: 0,
    paramAttrs: null,
    imageLog: false,
    imageData: null,   // { data, nrows, ncols } of the image currently shown
    cmap: 'Hot',
    cmapReverse: false,
    cmin: null,        // user colour limits in physical units; null = auto
    cmax: null,
    autoLimits: null,  // { min, max } used when cmin / cmax are null
    lightcurveLog: false,    // log intensity (y) axis
    lightcurveLogX: false,   // log time (x) axis
    currentTab: 'image',
    theme: 'dark',
  };

  // ── Theme colours ─────────────────────────────────────────────────────────────

  var THEMES = {
    dark: {
      paper:    '#1a1a2e',
      plotBg:   '#16213e',
      plotBgImg:'#1a1a2e',
      grid:     '#2a3a5e',
      zeroline: '#2a3a5e',
      axisText: '#8892a4',
      tickText: '#e0e0e0',
      line:     '#e94560',
    },
    light: {
      paper:    '#f0f2f5',
      plotBg:   '#ffffff',
      plotBgImg:'#f0f2f5',
      grid:     '#dce3ef',
      zeroline: '#b0bccc',
      axisText: '#5a6478',
      tickText: '#1a1a2e',
      line:     '#c0392b',
    },
  };

  function tc() { return THEMES[state.theme]; }

  // ── Colormaps ─────────────────────────────────────────────────────────────────

  // Evenly spaced stops; the matplotlib maps are sampled at 10 points.
  function evenStops(colors) {
    return colors.map(function (c, i) { return [i / (colors.length - 1), c]; });
  }

  var COLORMAPS = [
    { name: 'Hot',     stops: [[0, '#000000'], [0.3, '#e60000'], [0.6, '#ffd200'], [1, '#ffffff']] },
    { name: 'Inferno', stops: evenStops(['#000004', '#1b0c41', '#4a0c6b', '#781c6d', '#a52c60',
                                         '#cf4446', '#ed6925', '#fb9b06', '#f7d13d', '#fcffa4']) },
    { name: 'Magma',   stops: evenStops(['#000004', '#180f3d', '#440f76', '#721f81', '#9e2f7f',
                                         '#cd4071', '#f1605d', '#fd9668', '#feca8d', '#fcfdbf']) },
    { name: 'Plasma',  stops: evenStops(['#0d0887', '#46039f', '#7201a8', '#9c179e', '#bd3786',
                                         '#d8576b', '#ed7953', '#fb9f3a', '#fdca26', '#f0f921']) },
    { name: 'Viridis', stops: evenStops(['#440154', '#482878', '#3e4989', '#31688e', '#26828e',
                                         '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725']) },
    { name: 'Cividis', stops: evenStops(['#00224e', '#123570', '#3b496c', '#575d6d', '#707173',
                                         '#8a8678', '#a59c74', '#c3b369', '#e1cc55', '#fee838']) },
    { name: 'Greys',   stops: evenStops(['#000000', '#ffffff']) },
    { name: 'Jet',     stops: evenStops(['#000080', '#0000ff', '#0080ff', '#00ffff', '#80ff80',
                                         '#ffff00', '#ff8000', '#ff0000', '#800000']) },
  ];

  function colormapStops(name, reverse) {
    var cm = COLORMAPS.filter(function (m) { return m.name === name; })[0] || COLORMAPS[0];
    if (!reverse) return cm.stops;
    return cm.stops.map(function (s) { return [1 - s[0], s[1]]; }).reverse();
  }

  function cssGradient(stops) {
    return 'linear-gradient(to right, ' + stops.map(function (s) {
      return s[1] + ' ' + (s[0] * 100) + '%';
    }).join(', ') + ')';
  }

  // ── DOM refs ─────────────────────────────────────────────────────────────────

  var dropzone       = document.getElementById('dropzone');
  var dropInner      = document.getElementById('drop-inner');
  var dropOpenBtn    = document.getElementById('drop-open-btn');
  var openBtn        = document.getElementById('open-btn');
  var fileInput      = document.getElementById('file-input');
  var fnameLabel     = document.getElementById('fname-label');
  var obsSelectorWrap = document.getElementById('obs-selector-wrap');
  var obsSelector    = document.getElementById('obs-selector');
  var obsAttrsSection = document.getElementById('obs-attrs-section');
  var obsAttrsTable  = document.getElementById('obs-attrs-table');
  var paramsSection  = document.getElementById('params-section');
  var paramsTable    = document.getElementById('params-table');
  var tabBtns        = document.querySelectorAll('.tab-btn');
  var cscaleSection  = document.getElementById('cscale-section');
  var csScaleBtns    = document.querySelectorAll('#cs-scale button');
  var csMin          = document.getElementById('cs-min');
  var csMax          = document.getElementById('cs-max');
  var csFitBtn       = document.getElementById('cs-fit');
  var csAutoBtn      = document.getElementById('cs-auto');
  var cmapGrid       = document.getElementById('cmap-grid');
  var cmapReverseCb  = document.getElementById('cmap-reverse');
  var lcscaleSection = document.getElementById('lcscale-section');
  var lcScaleXBtns   = document.querySelectorAll('#lc-scale-x button');
  var lcScaleYBtns   = document.querySelectorAll('#lc-scale-y button');
  var plotImage      = document.getElementById('plot-image');
  var plotLc         = document.getElementById('plot-lc');
  var toast          = document.getElementById('toast');
  var themeBtn       = document.getElementById('theme-btn');

  // ── Init ──────────────────────────────────────────────────────────────────────

  document.addEventListener('DOMContentLoaded', function () {
    bindEvents();
  });

  function bindEvents() {
    // File open
    openBtn.addEventListener('click', function () { fileInput.click(); });
    dropOpenBtn.addEventListener('click', function () { fileInput.click(); });
    fileInput.addEventListener('change', function () {
      if (fileInput.files.length > 0) loadFile(fileInput.files[0]);
    });

    // Drag-and-drop on dropzone
    document.addEventListener('dragover', function (e) {
      e.preventDefault();
      dropInner.classList.add('drag-over');
    });
    document.addEventListener('dragleave', function (e) {
      if (!e.relatedTarget) dropInner.classList.remove('drag-over');
    });
    document.addEventListener('drop', function (e) {
      e.preventDefault();
      dropInner.classList.remove('drag-over');
      var files = e.dataTransfer.files;
      if (files.length > 0) loadFile(files[0]);
    });

    // Tabs
    tabBtns.forEach(function (btn) {
      btn.addEventListener('click', function () {
        switchTab(btn.dataset.tab);
      });
    });

    // Colour scale controls
    buildCmapGrid();
    bindScaleControl(csScaleBtns, function (log) {
      state.imageLog = log;
      if (state.file) renderImage(true);
    });
    csMin.addEventListener('change', function () { applyLimitInput(csMin, 'cmin'); });
    csMax.addEventListener('change', function () { applyLimitInput(csMax, 'cmax'); });
    [csMin, csMax].forEach(function (input) {
      input.addEventListener('keydown', function (e) {
        if (e.key === 'Enter') input.blur();   // blur fires 'change'
      });
    });
    csAutoBtn.addEventListener('click', function () {
      state.cmin = state.cmax = null;
      csMin.value = csMax.value = '';
      if (state.file) renderImage(true);
    });
    csFitBtn.addEventListener('click', fitLimitsToView);
    cmapReverseCb.addEventListener('change', function () {
      state.cmapReverse = cmapReverseCb.checked;
      updateCmapSwatches();
      if (state.file) renderImage(true);
    });

    // Light curve axis scales
    bindScaleControl(lcScaleXBtns, function (log) {
      state.lightcurveLogX = log;
      if (state.file) renderLightCurve(true);
    });
    bindScaleControl(lcScaleYBtns, function (log) {
      state.lightcurveLog = log;
      if (state.file) renderLightCurve(true);
    });

    // Observer selector
    obsSelector.addEventListener('change', function () {
      switchObserver(parseInt(obsSelector.value, 10));
    });

    // Theme toggle
    themeBtn.addEventListener('click', function () {
      state.theme = state.theme === 'dark' ? 'light' : 'dark';
      document.body.classList.toggle('theme-light', state.theme === 'light');
      themeBtn.textContent = state.theme === 'dark' ? '\u2600' : '\u263D';
      if (state.file) {
        renderImage(true);
        renderLightCurve(true);
      }
    });
  }

  // Linear/Log segmented control: onChange(isLog) runs after a click
  function bindScaleControl(buttons, onChange) {
    buttons.forEach(function (btn) {
      btn.addEventListener('click', function () {
        buttons.forEach(function (b) { b.classList.toggle('active', b === btn); });
        onChange(btn.dataset.scale === 'log');
      });
    });
  }

  // ── File loading ──────────────────────────────────────────────────────────────

  function loadFile(file) {
    state.fname = file.name;
    var reader = new FileReader();
    reader.onload = function (e) {
      try {
        var ab = e.target.result;
        var f = new hdf5.File(ab, file.name);
        state.file = f;
        parseFile(f);
      } catch (err) {
        showError('Failed to open file: ' + err.message);
      }
    };
    reader.onerror = function () {
      showError('Could not read file.');
    };
    reader.readAsArrayBuffer(file);
  }

  function parseFile(f) {
    // Detect observer groups
    var keys = f.keys;

    var obsKeys = keys.filter(function (k) { return /^obs\d{5}$/.test(k); });
    obsKeys.sort();

    if (obsKeys.length === 0) {
      showError('No observer groups (obs00001, …) found in this HDF5 file.');
      return;
    }

    var hasParams = keys.indexOf('parameters') !== -1;
    state.obsKeys = obsKeys;
    state.currentObs = 0;

    if (hasParams) {
      try {
        var pg = f.get('parameters');
        state.paramAttrs = pg ? pg.attrs : null;
      } catch (e) {
        console.warn('[CloudScat] could not read parameters attrs:', e);
        state.paramAttrs = null;
      }
    } else {
      state.paramAttrs = null;
    }

    fnameLabel.textContent = state.fname;
    populateSidebar();
    populateObsSelector();

    // Hide drop zone
    dropzone.classList.add('hidden');
    cscaleSection.style.display  = state.currentTab === 'image' ? 'block' : 'none';
    lcscaleSection.style.display = state.currentTab === 'lc' ? 'block' : 'none';

    switchObserver(0);
  }

  // ── Sidebar ───────────────────────────────────────────────────────────────────

  function populateObsSelector() {
    obsSelector.innerHTML = '';
    state.obsKeys.forEach(function (key, idx) {
      var opt = document.createElement('option');
      opt.value = idx;
      var attrs = {};
      try { attrs = state.file.get(key).attrs; } catch (e) {}
      var alt = attrs['altitude'] != null ? formatVal(attrs['altitude'] / 1000, 1) + ' km' : '?';
      var shift = attrs['shift'] != null ? formatVal(attrs['shift'] / 1000, 1) + ' km' : '?';
      opt.textContent = 'Observer ' + (idx + 1) + ' — alt ' + alt + ', shift ' + shift;
      obsSelector.appendChild(opt);
    });

    if (state.obsKeys.length > 1) {
      obsSelectorWrap.style.display = 'block';
    } else {
      obsSelectorWrap.style.display = 'none';
    }
  }

  function populateSidebar() {
    paramsTable.innerHTML = renderParamTable(state.paramAttrs, PARAM_FORMAT);
    paramsSection.style.display = 'block';
  }

  function populateObsAttrs(obsKey) {
    var attrs = null;
    try {
      var og = state.file.get(obsKey);
      attrs = og ? og.attrs : null;
    } catch (e) {
      console.warn('[CloudScat] could not read obs attrs:', e);
    }
    obsAttrsTable.innerHTML = renderParamTable(attrs, OBS_FORMAT);
    obsAttrsSection.style.display = 'block';
  }

  // Render a parameter table from an attrs object, using a format spec array.
  // formatSpec: array of { key, label, convert, unit, decimals }
  function renderParamTable(attrs, formatSpec) {
    if (!attrs) return '<tr><td colspan="2" class="param-empty">No attributes</td></tr>';

    var rows = '';
    var shownKeys = {};

    // First pass: render known keys in spec order
    formatSpec.forEach(function (spec) {
      var raw;
      try { raw = attrs[spec.key]; } catch (e) { return; }
      if (raw == null) return;
      shownKeys[spec.key] = true;
      var val = spec.convert ? spec.convert(scalarVal(raw)) : scalarVal(raw);
      var formatted = spec.decimals != null
        ? formatVal(val, spec.decimals)
        : formatAuto(val);
      rows += '<tr><td>' + esc(spec.label) + '</td><td>' + esc(formatted) + (spec.unit ? '&thinsp;' + esc(spec.unit) : '') + '</td></tr>';
    });

    // Second pass: all remaining keys not already shown
    // Try multiple enumeration methods since jsfive attrs may not be a plain object
    var allKeys = [];
    try { allKeys = allKeys.concat(Object.keys(attrs)); } catch (e) {}
    try {
      for (var k in attrs) {
        if (Object.prototype.hasOwnProperty.call(attrs, k) && allKeys.indexOf(k) === -1) {
          allKeys.push(k);
        }
      }
    } catch (e) {}

    allKeys.forEach(function (k) {
      if (shownKeys[k]) return;
      var raw;
      try { raw = attrs[k]; } catch (e) { return; }
      if (raw == null) return;
      rows += '<tr><td>' + esc(k) + '</td><td>' + esc(formatAuto(raw)) + '</td></tr>';
    });

    if (!rows) {
      rows = '<tr><td colspan="2" class="param-empty">No attributes found</td></tr>';
    }
    return rows;
  }

  // Extract scalar from a possible 1-element typed array
  function scalarVal(v) {
    if (v != null && typeof v === 'object' && v.length === 1) return v[0];
    return v;
  }

  // HTML-escape helper
  function esc(s) {
    return String(s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');
  }

  // ── Format helpers ────────────────────────────────────────────────────────────

  // Parameter formatting spec — keys match Julia Params struct field names
  var PARAM_FORMAT = [
    { key: 'λ',             label: 'Wavelength λ',        convert: function(v){ return v * 1e9; },    unit: 'nm',   decimals: 1 },
    { key: 'N',             label: 'Photons N',            convert: null,                             unit: '',     decimals: 0 },
    { key: 'H',             label: 'Scale height H',       convert: function(v){ return v / 1000; },  unit: 'km',   decimals: 2 },
    { key: 'σray',          label: 'σ_ray',                convert: null,                             unit: 'm²',   decimals: null },
    { key: 'nair',          label: 'Air density n_air',    convert: null,                             unit: 'm⁻³',  decimals: null },
    { key: 'νray_ground',   label: 'ν_ray (ground)',       convert: null,                             unit: 'm⁻¹',  decimals: null },
    { key: 'νraymax',       label: 'ν_ray (max)',          convert: null,                             unit: 'm⁻¹',  decimals: null },
    { key: 'zmin',          label: 'z_min',                convert: function(v){ return v / 1000; },  unit: 'km',   decimals: 2 },
    { key: 'weight_min',    label: 'Weight min',           convert: null,                             unit: '',     decimals: 4 },
    { key: 'max_iter',      label: 'Max iterations',       convert: null,                             unit: '',     decimals: 0 },
    // Mie parameters (present in some output files)
    { key: 'radius',        label: 'Droplet radius',       convert: function(v){ return v * 1e6; },   unit: 'µm',   decimals: 2 },
    { key: 'nscat',         label: 'Number density n',     convert: null,                             unit: 'm⁻³',  decimals: null },
    { key: 'g',             label: 'Asymmetry g',          convert: null,                             unit: '',     decimals: 4 },
    { key: 'Qext',          label: 'Q_ext',                convert: null,                             unit: '',     decimals: 4 },
    { key: 'ω₀',            label: 'Single scatter ω₀',   convert: null,                             unit: '',     decimals: 4 },
    { key: 'source_a',      label: 'Source A',             convert: null,                             unit: 'm',    decimals: null },
    { key: 'source_b',      label: 'Source B',             convert: null,                             unit: 'm',    decimals: null },
  ];

  var OBS_FORMAT = [
    { key: 'altitude', label: 'Altitude',   convert: function(v){ return v / 1000; }, unit: 'km', decimals: 2 },
    { key: 'shift',    label: 'Shift',      convert: function(v){ return v / 1000; }, unit: 'km', decimals: 2 },
    { key: 'delay',    label: 'Delay',      convert: function(v){ return v * 1000; }, unit: 'ms', decimals: 4 },
  ];

  function formatVal(v, decimals) {
    if (decimals === 0) return Math.round(v).toLocaleString();
    return Number(v).toFixed(decimals);
  }

  function formatAuto(v) {
    if (typeof v === 'number') {
      if (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-2 && v !== 0)) {
        return v.toExponential(3);
      }
      return Number(v.toPrecision(5)).toString();
    }
    return String(v);
  }

  // ── Observer switching ────────────────────────────────────────────────────────

  function switchObserver(idx) {
    state.currentObs = idx;
    obsSelector.value = idx;
    var obsKey = state.obsKeys[idx];
    populateObsAttrs(obsKey);
    renderImage();
    renderLightCurve();
  }

  // ── Tab switching ─────────────────────────────────────────────────────────────

  function switchTab(tab) {
    state.currentTab = tab;
    tabBtns.forEach(function (btn) {
      btn.classList.toggle('active', btn.dataset.tab === tab);
    });
    document.querySelectorAll('.plot-pane').forEach(function (pane) {
      pane.classList.toggle('active', pane.id === 'pane-' + tab);
    });

    // Show the controls that apply to this tab
    cscaleSection.style.display  = (tab === 'image' && state.file) ? 'block' : 'none';
    lcscaleSection.style.display = (tab === 'lc' && state.file) ? 'block' : 'none';

    // Relayout the active plot so Plotly fills the new dimensions
    if (tab === 'image') {
      Plotly.Plots.resize(plotImage);
    } else {
      Plotly.Plots.resize(plotLc);
    }
  }

  // ── Image rendering ───────────────────────────────────────────────────────────

  function renderImage(preserveRange) {
    var obsKey = state.obsKeys[state.currentObs];
    var ds;
    try {
      ds = state.file.get(obsKey + '/image');
    } catch (e) {
      showError('Could not load image dataset: ' + e.message);
      return;
    }

    var data = ds.value;       // Float64Array, row-major
    var shape = ds.shape;      // [npx, npx]
    if (!shape || shape.length < 2) {
      showError('Unexpected image shape.');
      return;
    }

    var nrows = shape[0];
    var ncols = shape[1];
    state.imageData = { data: data, nrows: nrows, ncols: ncols };

    var dataMin = Infinity;
    var dataMax = -Infinity;
    for (var i = 0; i < data.length; i++) {
      if (data[i] < dataMin) dataMin = data[i];
      if (data[i] > dataMax) dataMax = data[i];
    }
    var lim = colorLimits(dataMin, dataMax);

    // Build nested array (rows × cols), applying log scale if requested
    var z = [];
    for (var r = 0; r < nrows; r++) {
      var row = new Array(ncols);
      for (var c = 0; c < ncols; c++) {
        var v = data[r * ncols + c];
        row[c] = state.imageLog ? (v > 0 ? Math.log10(v) : NaN) : v;
      }
      z.push(row);
    }

    var savedX = null, savedY = null;
    if (preserveRange && plotImage._fullLayout) {
      savedX = plotImage._fullLayout.xaxis.range.slice();
      savedY = plotImage._fullLayout.yaxis.range.slice();
    }

    // Plot area for the ranges about to be shown, so that the first draw is
    // already square; fitImageArea() corrects it after zooms and resizes.
    var dom = squareDomains(savedX ? Math.abs(savedX[1] - savedX[0]) : ncols,
                            savedY ? Math.abs(savedY[1] - savedY[0]) : nrows);
    var cb = colorbarGeometry(dom.x, dom.y);
    var trace = {
      type: 'heatmap',
      z: z,
      colorscale: colormapStops(state.cmap, state.cmapReverse),
      zsmooth: false,
      showscale: true,
      colorbar: {
        thickness: 14,
        x: cb.x,
        y: cb.y,
        len: cb.len,
        tickfont: { color: tc().tickText, size: 11 },
        tickcolor: tc().tickText,
        outlinecolor: tc().grid,
        exponentformat: 'power',
        title: {
          text: state.imageLog ? 'log₁₀(ph / sr / m² / src ph)' : 'ph / sr / m² / src ph',
          font: { color: tc().axisText, size: 11 },
        },
      },
    };

    if (lim) {
      trace.zauto = false;
      trace.zmin = state.imageLog ? Math.log10(lim.min) : lim.min;
      trace.zmax = state.imageLog ? Math.log10(lim.max) : lim.max;
    }

    // The axes are deliberately not linked with scaleanchor: Plotly would then
    // force the zoom box to the current aspect ratio while dragging.  Square
    // pixels are kept instead by sizing the axis domains (squareDomains).
    var layout = {
      paper_bgcolor: tc().plotBgImg,
      plot_bgcolor:  tc().plotBgImg,
      margin: { t: 20, r: 90, b: 60, l: 90 },
      xaxis: {
        title: { text: 'px', font: { color: tc().axisText } },
        color: tc().axisText,
        gridcolor: tc().grid,
        zerolinecolor: tc().zeroline,
        exponentformat: 'power',
        domain: dom.x,
        range: savedX || undefined,
      },
      yaxis: {
        title: { text: 'py', font: { color: tc().axisText } },
        color: tc().axisText,
        gridcolor: tc().grid,
        zerolinecolor: tc().zeroline,
        exponentformat: 'power',
        domain: dom.y,
        // py grows along the camera's e2 axis (upwards for PointTo), so it
        // increases upwards on screen.
        autorange: savedY ? false : true,
        range: savedY || undefined,
      },
    };

    var config = { responsive: true, displaylogo: false };

    Plotly.react(plotImage, [trace], layout, config);

    if (!plotImage._fitHooked) {
      plotImage.on('plotly_afterplot', fitImageArea);
      plotImage._fitHooked = true;
    }
  }

  // Axis domains (fractions of the plot area) that show an xspan × yspan
  // data rectangle with square pixels, centred, as large as fits.
  function squareDomains(xspan, yspan) {
    var fl = plotImage._fullLayout;
    var w = fl && fl._size ? fl._size.w : plotImage.clientWidth - 180;
    var h = fl && fl._size ? fl._size.h : plotImage.clientHeight - 80;
    var scale = Math.min(w / xspan, h / yspan);
    var fx = Math.min(1, scale * xspan / w);
    var fy = Math.min(1, scale * yspan / h);
    return { x: [(1 - fx) / 2, (1 + fx) / 2], y: [(1 - fy) / 2, (1 + fy) / 2] };
  }

  // Colorbar placement that hugs the image area, in paper coordinates.  Keeps
  // a readable length when the image is a thin horizontal strip.
  function colorbarGeometry(xd, yd) {
    return { x: xd[1], y: (yd[0] + yd[1]) / 2, len: Math.max(yd[1] - yd[0], 0.5) };
  }

  // After every draw (box zoom, reset, resize) resize the plot area to keep
  // pixels square for the new ranges, then move the colorbar next to it.
  // The tolerance checks stop the relayout / restyle from re-triggering
  // themselves.
  function fitImageArea() {
    var fl = plotImage._fullLayout;
    if (!fl || !fl.xaxis || !fl._size) return;
    var xr = fl.xaxis.range, yr = fl.yaxis.range;
    var dom = squareDomains(Math.abs(xr[1] - xr[0]), Math.abs(yr[1] - yr[0]));
    var xd = fl.xaxis.domain, yd = fl.yaxis.domain;
    if (Math.abs(xd[0] - dom.x[0]) > 1e-4 || Math.abs(xd[1] - dom.x[1]) > 1e-4 ||
        Math.abs(yd[0] - dom.y[0]) > 1e-4 || Math.abs(yd[1] - dom.y[1]) > 1e-4) {
      // Pass the ranges too, so that an autoranged view is not recomputed
      Plotly.relayout(plotImage, {
        'xaxis.domain': dom.x, 'yaxis.domain': dom.y,
        'xaxis.range': xr.slice(), 'yaxis.range': yr.slice(),
      });
      return;   // the colorbar is synced on the next afterplot
    }

    var g = colorbarGeometry(xd, yd);
    var cb = plotImage.data && plotImage.data[0] && plotImage.data[0].colorbar;
    if (!cb) return;
    if (Math.abs((cb.x || 0) - g.x) < 1e-3 &&
        Math.abs((cb.y || 0) - g.y) < 1e-3 &&
        Math.abs((cb.len || 0) - g.len) < 1e-3) return;
    Plotly.restyle(plotImage, { 'colorbar.x': g.x, 'colorbar.y': g.y, 'colorbar.len': g.len }, [0]);
  }

  // ── Colour scale controls ─────────────────────────────────────────────────────

  // Effective colour limits in physical units, combining user values with the
  // automatic ones (full range, or 4 decades below the max in log scale).
  // Flags inputs whose values cannot be used.
  function colorLimits(dataMin, dataMax) {
    var autoMin, autoMax;
    if (state.imageLog) {
      autoMax = dataMax;
      autoMin = dataMax * 1e-4;
    } else {
      autoMin = dataMin;
      autoMax = dataMax;
    }
    var valid = state.imageLog ? dataMax > 0 : isFinite(dataMin) && isFinite(dataMax);
    state.autoLimits = valid ? { min: autoMin, max: autoMax } : null;

    var userMin = state.cmin, userMax = state.cmax;
    var minBad = userMin != null && state.imageLog && userMin <= 0;
    if (minBad) userMin = null;

    var lo = userMin != null ? userMin : autoMin;
    var hi = userMax != null ? userMax : autoMax;
    var orderBad = (userMin != null || userMax != null) && !(lo < hi);
    if (orderBad) {
      lo = autoMin;
      hi = autoMax;
    }

    flagInput(csMin, minBad ? 'Log scale needs Min > 0; using auto.'
                    : orderBad ? 'Min must be below Max; using auto.' : '');
    flagInput(csMax, orderBad ? 'Min must be below Max; using auto.' : '');
    csMin.placeholder = state.autoLimits ? 'auto (' + formatAuto(autoMin) + ')' : 'auto';
    csMax.placeholder = state.autoLimits ? 'auto (' + formatAuto(autoMax) + ')' : 'auto';

    if (!valid && userMin == null && userMax == null) return null;
    return { min: lo, max: hi };
  }

  function flagInput(input, msg) {
    input.classList.toggle('invalid', !!msg);
    input.title = msg;
  }

  // Read a Min/Max input into state[key]; empty means auto.
  function applyLimitInput(input, key) {
    var s = input.value.trim().replace('−', '-');
    if (s === '') {
      state[key] = null;
    } else {
      var v = Number(s);
      if (!isFinite(v)) {
        input.value = state[key] != null ? formatAuto(state[key]) : '';
        showError('Not a number: ' + s);
        return;
      }
      state[key] = v;
    }
    if (state.file) renderImage(true);
  }

  // Set the limits to the data range of the pixels currently in view.  In log
  // scale the min is kept within 4 decades of the max, as in auto mode, since
  // Monte Carlo images contain a few pixels with vanishingly small weights.
  function fitLimitsToView() {
    var img = state.imageData;
    var fl = plotImage._fullLayout;
    if (!img || !fl || !fl.xaxis) return;
    var xr = fl.xaxis.range, yr = fl.yaxis.range;
    var c0 = Math.max(0, Math.ceil(Math.min(xr[0], xr[1])));
    var c1 = Math.min(img.ncols - 1, Math.floor(Math.max(xr[0], xr[1])));
    var r0 = Math.max(0, Math.ceil(Math.min(yr[0], yr[1])));
    var r1 = Math.min(img.nrows - 1, Math.floor(Math.max(yr[0], yr[1])));

    var lo = Infinity, hi = -Infinity;
    for (var r = r0; r <= r1; r++) {
      for (var c = c0; c <= c1; c++) {
        var v = img.data[r * img.ncols + c];
        if (state.imageLog && !(v > 0)) continue;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
    }
    if (!isFinite(lo) || !isFinite(hi)) {
      showError(state.imageLog ? 'No positive pixels in view.' : 'No pixels in view.');
      return;
    }
    if (state.imageLog) lo = Math.max(lo, hi * 1e-4);
    state.cmin = lo;
    state.cmax = hi;
    csMin.value = formatAuto(lo);
    csMax.value = formatAuto(hi);
    renderImage(true);
  }

  function buildCmapGrid() {
    COLORMAPS.forEach(function (cm) {
      var btn = document.createElement('button');
      btn.className = 'cmap-btn' + (cm.name === state.cmap ? ' active' : '');
      btn.dataset.name = cm.name;
      btn.title = cm.name;
      btn.innerHTML = '<span class="swatch"></span><span class="cmap-name">' + esc(cm.name) + '</span>';
      btn.addEventListener('click', function () {
        state.cmap = cm.name;
        cmapGrid.querySelectorAll('.cmap-btn').forEach(function (b) {
          b.classList.toggle('active', b === btn);
        });
        if (state.file) renderImage(true);
      });
      cmapGrid.appendChild(btn);
    });
    updateCmapSwatches();
  }

  function updateCmapSwatches() {
    cmapGrid.querySelectorAll('.cmap-btn').forEach(function (btn) {
      btn.querySelector('.swatch').style.background =
        cssGradient(colormapStops(btn.dataset.name, state.cmapReverse));
    });
  }

  // ── Light curve rendering ─────────────────────────────────────────────────────

  function renderLightCurve(preserveRange) {
    var obsKey = state.obsKeys[state.currentObs];
    var dsT, dsTL;
    try {
      dsT  = state.file.get(obsKey + '/t');
      dsTL = state.file.get(obsKey + '/timeline');
    } catch (e) {
      showError('Could not load light curve datasets: ' + e.message);
      return;
    }

    var t  = dsT.value;   // Float64Array, seconds
    var tl = dsTL.value;  // Float64Array

    // Get propagation delay from observer attrs
    var delay = 0;
    try {
      var obsAttrs = state.file.get(obsKey).attrs;
      if (obsAttrs['delay'] != null) delay = obsAttrs['delay'];
    } catch (e) {}

    // Shift time and convert to ms
    var tMs = new Array(t.length);
    for (var i = 0; i < t.length; i++) {
      tMs[i] = (t[i] - delay) * 1000;
    }

    var trace = {
      type: 'scatter',
      mode: 'lines',
      x: tMs,
      y: Array.from(tl),
      line: { color: tc().line, width: 1.5 },
      name: obsKey,
    };

    var xType = state.lightcurveLogX ? 'log' : 'linear';

    var layout = {
      paper_bgcolor: tc().paper,
      plot_bgcolor:  tc().plotBg,
      margin: { t: 20, r: 50, b: 80, l: 140 },
      xaxis: {
        title: { text: 'Time (ms)', font: { color: tc().axisText }, standoff: 16 },
        color: tc().axisText,
        gridcolor: tc().grid,
        zerolinecolor: tc().zeroline,
        exponentformat: 'power',
        type: xType,
        // A range is only meaningful for the axis type it was taken from
        range: (preserveRange && plotLc._fullLayout && plotLc._fullLayout.xaxis.type === xType)
          ? plotLc._fullLayout.xaxis.range.slice() : undefined,
      },
      yaxis: (function() {
        var ax = {
          title: { text: 'photons / m² / s / source photon', font: { color: tc().axisText }, standoff: 28 },
          color: tc().axisText,
          gridcolor: tc().grid,
          zerolinecolor: tc().zeroline,
          exponentformat: 'power',
          type: state.lightcurveLog ? 'log' : 'linear',
        };
        if (state.lightcurveLog) {
          var maxVal = -Infinity;
          for (var i = 0; i < tl.length; i++) {
            if (tl[i] > maxVal) maxVal = tl[i];
          }
          if (maxVal > 0) {
            var logMax = Math.log10(maxVal);
            ax.range = [logMax - 4, logMax];
          }
        }
        return ax;
      })(),
      showlegend: false,
    };

    var config = { responsive: true, displaylogo: false };

    Plotly.react(plotLc, [trace], layout, config);
  }

  // ── Error toast ───────────────────────────────────────────────────────────────

  var toastTimer = null;

  function showError(msg) {
    console.error('[CloudScat Viewer]', msg);
    toast.textContent = msg;
    toast.classList.add('show');
    if (toastTimer) clearTimeout(toastTimer);
    toastTimer = setTimeout(function () {
      toast.classList.remove('show');
    }, 5000);
  }

  // ── Resize handling ───────────────────────────────────────────────────────────

  window.addEventListener('resize', function () {
    if (!state.file) return;
    Plotly.Plots.resize(plotImage);
    Plotly.Plots.resize(plotLc);
  });

})();
