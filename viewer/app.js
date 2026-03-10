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
    lightcurveLog: false,
    currentTab: 'image',
  };

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
  var logImageWrap   = document.getElementById('log-image-wrap');
  var logImageCb     = document.getElementById('log-image-cb');
  var logLcWrap      = document.getElementById('log-lc-wrap');
  var logLcCb        = document.getElementById('log-lc-cb');
  var plotImage      = document.getElementById('plot-image');
  var plotLc         = document.getElementById('plot-lc');
  var toast          = document.getElementById('toast');

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

    // Log scale toggles
    logImageCb.addEventListener('change', function () {
      state.imageLog = logImageCb.checked;
      if (state.file) renderImage(true);
    });
    logLcCb.addEventListener('change', function () {
      state.lightcurveLog = logLcCb.checked;
      if (state.file) renderLightCurve(true);
    });

    // Observer selector
    obsSelector.addEventListener('change', function () {
      switchObserver(parseInt(obsSelector.value, 10));
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
    console.log('[CloudScat] top-level keys:', keys);

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
        console.log('[CloudScat] parameters attrs:', state.paramAttrs);
        console.log('[CloudScat] parameters attrs keys:', state.paramAttrs ? Object.keys(state.paramAttrs) : 'null');
        if (state.paramAttrs) {
          var sampleKeys = ['N', 'λ', 'H', 'zmin'];
          sampleKeys.forEach(function(k) {
            var v = state.paramAttrs[k];
            console.log('[CloudScat] attrs[' + k + '] =', v, '| type:', typeof v, '| constructor:', v != null && v.constructor ? v.constructor.name : 'n/a');
          });
        }
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
      console.log('[CloudScat] obs attrs:', obsKey, attrs, attrs ? Object.keys(attrs) : 'null');
    } catch (e) {
      console.warn('[CloudScat] could not read obs attrs:', e);
    }
    obsAttrsTable.innerHTML = renderParamTable(attrs, OBS_FORMAT);
    obsAttrsSection.style.display = 'block';
  }

  // Render a parameter table from an attrs object, using a format spec array.
  // formatSpec: array of { key, label, convert, unit, decimals }
  function renderParamTable(attrs, formatSpec) {
    if (!attrs) return '<tr><td colspan="2" style="color:#8892a4;font-style:italic">No attributes</td></tr>';

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
      rows = '<tr><td colspan="2" style="color:#8892a4;font-style:italic">No attributes found</td></tr>';
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

    // Show relevant log toggle
    logImageWrap.style.display = tab === 'image' ? 'flex' : 'none';
    logLcWrap.style.display    = tab === 'lc'    ? 'flex' : 'none';

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

    // Build nested array (rows × cols), applying log scale if requested
    var z = [];
    var globalMin = Infinity;
    var globalMax = -Infinity;

    if (state.imageLog) {
      // Pre-compute log values
      var zLog = [];
      for (var r = 0; r < nrows; r++) {
        var row = [];
        for (var c = 0; c < ncols; c++) {
          var v = data[r * ncols + c];
          var lv = (v > 0) ? Math.log10(v) : NaN;
          row.push(lv);
          if (!isNaN(lv)) {
            if (lv < globalMin) globalMin = lv;
            if (lv > globalMax) globalMax = lv;
          }
        }
        zLog.push(row);
      }
      z = zLog;
    } else {
      for (var r = 0; r < nrows; r++) {
        var row = [];
        for (var c = 0; c < ncols; c++) {
          var v = data[r * ncols + c];
          row.push(v);
          if (v < globalMin) globalMin = v;
          if (v > globalMax) globalMax = v;
        }
        z.push(row);
      }
    }

    var trace = {
      type: 'heatmap',
      z: z,
      colorscale: 'Hot',
      zsmooth: false,
      showscale: true,
      colorbar: {
        thickness: 14,
        len: 0.8,
        tickfont: { color: '#e0e0e0', size: 11 },
        tickcolor: '#e0e0e0',
        outlinecolor: '#2a3a5e',
        exponentformat: 'power',
        title: {
          text: state.imageLog ? 'log₁₀(ph / sr / m² / src ph)' : 'ph / sr / m² / src ph',
          font: { color: '#8892a4', size: 11 },
        },
      },
    };

    if (state.imageLog && isFinite(globalMax)) {
      trace.zmax = globalMax;
      trace.zmin = globalMax - 4;
    }

    var savedX = null, savedY = null;
    if (preserveRange && plotImage._fullLayout) {
      savedX = plotImage._fullLayout.xaxis.range.slice();
      savedY = plotImage._fullLayout.yaxis.range.slice();
    }

    var layout = {
      paper_bgcolor: '#1a1a2e',
      plot_bgcolor:  '#1a1a2e',
      margin: { t: 20, r: 90, b: 60, l: 90 },
      xaxis: {
        title: { text: 'px', font: { color: '#8892a4' } },
        color: '#8892a4',
        gridcolor: '#2a3a5e',
        zerolinecolor: '#2a3a5e',
        exponentformat: 'power',
        scaleanchor: 'y',
        scaleratio: 1,
        range: savedX || undefined,
      },
      yaxis: {
        title: { text: 'py', font: { color: '#8892a4' } },
        color: '#8892a4',
        gridcolor: '#2a3a5e',
        zerolinecolor: '#2a3a5e',
        exponentformat: 'power',
        autorange: savedY ? false : 'reversed',
        range: savedY || undefined,
      },
    };

    var config = { responsive: true, displaylogo: false };

    Plotly.react(plotImage, [trace], layout, config);

    // Show log toggle
    logImageWrap.style.display = (state.currentTab === 'image') ? 'flex' : 'none';
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
      line: { color: '#e94560', width: 1.5 },
      name: obsKey,
    };

    var layout = {
      paper_bgcolor: '#1a1a2e',
      plot_bgcolor:  '#16213e',
      margin: { t: 20, r: 50, b: 80, l: 140 },
      xaxis: {
        title: { text: 'Time (ms)', font: { color: '#8892a4' }, standoff: 16 },
        color: '#8892a4',
        gridcolor: '#2a3a5e',
        zerolinecolor: '#2a3a5e',
        exponentformat: 'power',
        range: (preserveRange && plotLc._fullLayout) ? plotLc._fullLayout.xaxis.range.slice() : undefined,
      },
      yaxis: (function() {
        var ax = {
          title: { text: 'photons / m² / s / source photon', font: { color: '#8892a4' }, standoff: 28 },
          color: '#8892a4',
          gridcolor: '#2a3a5e',
          zerolinecolor: '#2a3a5e',
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

    // Show log toggle if on lc tab
    logLcWrap.style.display = (state.currentTab === 'lc') ? 'flex' : 'none';
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

  // Initialize tab controls visibility
  logImageWrap.style.display = 'none';
  logLcWrap.style.display    = 'none';

})();
