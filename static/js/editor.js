(function () {
  "use strict";
  var root = document.getElementById("gfv-result");
  if (!root) return;

  var cfg = {
    statusUrl: root.dataset.statusUrl,
    svgUrl: root.dataset.svgUrl,
    restyleUrl: root.dataset.restyleUrl,
    downloadBase: root.dataset.downloadBase,
    chartBase: root.dataset.chartBase,
    bundleUrl: root.dataset.bundleUrl,
    analyticsBase: root.dataset.analyticsBase
  };

  var el = {
    state: document.getElementById("gfv-state"),
    stateText: root.querySelector(".state-text"),
    error: document.getElementById("gfv-error"),
    warnings: document.getElementById("gfv-warnings"),
    workspace: document.getElementById("gfv-workspace"),
    figure: document.getElementById("gfv-figure"),
    famList: document.getElementById("gfv-fam-list"),
    famAll: document.getElementById("gfv-fam-all"),
    chromList: document.getElementById("gfv-chrom-list"),
    chromAll: document.getElementById("gfv-chrom-all"),
    showUnplaced: document.getElementById("gfv-show-unplaced"),
    genesChart: document.getElementById("gfv-genes-chart"),
    coloc: document.getElementById("gfv-coloc"),
    labelMode: document.getElementById("gfv-label-mode"),
    labelSize: document.getElementById("gfv-label-size"),
    labelSizeVal: document.getElementById("gfv-label-size-val"),
    legendSize: document.getElementById("gfv-legend-size"),
    legendSizeVal: document.getElementById("gfv-legend-size-val"),
    font: document.getElementById("gfv-font"),
    legendLoc: document.getElementById("gfv-legend-loc"),
    orientation: document.getElementById("gfv-orientation"),
    tickStyle: document.getElementById("gfv-tick-style"),
    centromeres: document.getElementById("gfv-centromeres"),
    showTitles: document.getElementById("gfv-show-titles"),
    apply: document.getElementById("gfv-apply"),
    reset: document.getElementById("gfv-reset"),
    applyMsg: document.getElementById("gfv-apply-msg"),
    analytics: document.getElementById("gfv-analytics"),
    analyticsSummary: document.getElementById("gfv-analytics-summary"),
    analyticsFigs: document.getElementById("gfv-analytics-figs"),
    genesSection: document.getElementById("gfv-genes-section"),
    dlBundle: document.getElementById("gfv-dl-bundle"),
    dlBundleAll: document.getElementById("gfv-dl-bundle-all"),
    zoomIn: document.getElementById("gfv-zoom-in"),
    zoomOut: document.getElementById("gfv-zoom-out"),
    zoomReset: document.getElementById("gfv-zoom-reset"),
    dlEdited: document.getElementById("gfv-dl-edited")
  };

  var status = null;
  var colorOverrides = {};
  var zoom = 1;

  function slug(s) {
    return String(s).replace(/[^A-Za-z0-9_.\-]/g, "_");
  }

  // ---------------------------------------------------------------- polling
  function poll() {
    fetch(cfg.statusUrl)
      .then(function (r) { return r.json(); })
      .then(function (s) {
        status = s;
        if (s.state === "error") {
          showError(s.error || "Rendering failed.");
        } else if (s.state === "done") {
          onDone(s);
        } else {
          el.stateText.textContent =
            s.state === "running" ? "Rendering your figure…" : "Queued…";
          setTimeout(poll, 1200);
        }
      })
      .catch(function () {
        setTimeout(poll, 2500);
      });
  }

  function showError(text) {
    el.state.hidden = true;
    el.error.hidden = false;
    el.error.textContent = text;
  }

  function onDone(s) {
    el.state.hidden = true;
    el.workspace.hidden = false;
    if (s.warnings && s.warnings.length) {
      el.warnings.hidden = false;
      el.warnings.innerHTML =
        "<b>Notes:</b><ul>" +
        s.warnings.map(function (w) { return "<li>" + escapeHtml(w) + "</li>"; }).join("") +
        "</ul>";
    }
    seedControls(s.style || {});
    buildFamilyList(s.color_map || {}, (s.style || {}).only_families || []);
    buildChromList(s.chromosomes || [], (s.style || {}).only_chromosomes || []);
    wireDownloads();
    loadSvg();
    if (s.has_analytics) loadAnalytics(s.color_map || {});
  }

  function escapeHtml(t) {
    var d = document.createElement("div");
    d.textContent = t;
    return d.innerHTML;
  }

  // ---------------------------------------------------------------- controls
  function seedControls(style) {
    if (style.label_mode) el.labelMode.value = style.label_mode;
    if (style.label_font_size) {
      el.labelSize.value = style.label_font_size;
      el.labelSizeVal.textContent = style.label_font_size;
    }
    if (style.legend_font_size) {
      el.legendSize.value = style.legend_font_size;
      el.legendSizeVal.textContent = style.legend_font_size;
    }
    if (style.legend_font_family) el.font.value = style.legend_font_family;
    if (style.legend_location) {
      el.legendLoc.value = style.legend_show === false ? "none" : style.legend_location;
    }
    if (style.tick_style) el.tickStyle.value = style.tick_style;
    if (style.orientation) el.orientation.value = style.orientation;
    el.centromeres.checked = !!style.show_centromeres;
    el.showUnplaced.checked = !!style.show_unplaced;
    if (el.showTitles) el.showTitles.checked = style.show_titles !== false;
  }

  el.famAll.addEventListener("change", function () {
    var on = el.famAll.checked;
    el.famList.querySelectorAll('input[type="checkbox"]').forEach(function (cb) {
      if (cb.checked !== on) {
        cb.checked = on;
        toggleFamily(cb.dataset.fam, on);
      }
    });
  });

  function syncFamAll() {
    var boxes = el.famList.querySelectorAll('input[type="checkbox"]');
    var checked = 0;
    boxes.forEach(function (cb) { if (cb.checked) checked++; });
    el.famAll.checked = checked === boxes.length;
    el.famAll.indeterminate = checked > 0 && checked < boxes.length;
  }

  // ------------------------------------------------------------ chromosomes
  function buildChromList(chroms, selected) {
    if (!el.chromList) return;
    el.chromList.innerHTML = "";
    var pick = selected && selected.length ? {} : null;
    (selected || []).forEach(function (c) { pick[c] = true; });
    chroms.forEach(function (name) {
      var row = document.createElement("div");
      row.className = "fam-row";
      var cb = document.createElement("input");
      cb.type = "checkbox";
      cb.checked = pick ? !!pick[name] : true;
      cb.dataset.chrom = name;
      cb.addEventListener("change", syncChromAll);
      var label = document.createElement("span");
      label.textContent = name;
      row.appendChild(cb);
      row.appendChild(label);
      el.chromList.appendChild(row);
    });
    syncChromAll();
  }

  function syncChromAll() {
    if (!el.chromAll) return;
    var boxes = el.chromList.querySelectorAll('input[type="checkbox"]');
    var checked = 0;
    boxes.forEach(function (cb) { if (cb.checked) checked++; });
    el.chromAll.checked = boxes.length > 0 && checked === boxes.length;
    el.chromAll.indeterminate = checked > 0 && checked < boxes.length;
  }

  if (el.chromAll) {
    el.chromAll.addEventListener("change", function () {
      var on = el.chromAll.checked;
      el.chromList.querySelectorAll('input[type="checkbox"]').forEach(function (cb) {
        cb.checked = on;
      });
      el.chromAll.indeterminate = false;
    });
  }

  // [] means "every chromosome" (also when the user leaves none ticked)
  function collectChroms() {
    if (!el.chromList) return [];
    var boxes = el.chromList.querySelectorAll('input[type="checkbox"]');
    var picked = [];
    boxes.forEach(function (cb) { if (cb.checked) picked.push(cb.dataset.chrom); });
    if (picked.length === 0 || picked.length === boxes.length) return [];
    return picked;
  }

  el.labelSize.addEventListener("input", function () {
    el.labelSizeVal.textContent = el.labelSize.value;
    applyLabelSize();
  });
  el.legendSize.addEventListener("input", function () {
    el.legendSizeVal.textContent = el.legendSize.value;
  });
  el.labelMode.addEventListener("change", function () {
    // live: hide labels if "none", otherwise a full regenerate is needed
    var hide = el.labelMode.value === "none";
    svgSelectAll('[id^="gfv-label-"]').forEach(function (n) { n.style.display = hide ? "none" : ""; });
  });

  function buildFamilyList(colorMap, selected) {
    el.famList.innerHTML = "";
    var pick = selected && selected.length ? {} : null;
    (selected || []).forEach(function (f) { pick[f] = true; });
    Object.keys(colorMap).forEach(function (fam) {
      var rgb = colorMap[fam];
      var hex = rgbToHex(rgb);
      var row = document.createElement("div");
      row.className = "fam-row";

      var cb = document.createElement("input");
      cb.type = "checkbox";
      cb.checked = pick ? !!pick[fam] : true;
      cb.dataset.fam = fam;
      cb.addEventListener("change", function () {
        toggleFamily(fam, cb.checked);
        syncFamAll();
      });

      var sw = document.createElement("input");
      sw.type = "color";
      sw.value = hex;
      sw.title = "Recolour " + fam;
      sw.addEventListener("input", function () { recolourFamily(fam, sw.value); });

      var name = document.createElement("span");
      name.textContent = fam;

      row.appendChild(cb);
      row.appendChild(sw);
      row.appendChild(name);
      el.famList.appendChild(row);
    });
    syncFamAll();
  }

  // [] means "every family" (also when the user leaves none ticked)
  function collectFamilies() {
    var boxes = el.famList.querySelectorAll('input[type="checkbox"]');
    var picked = [];
    boxes.forEach(function (cb) { if (cb.checked) picked.push(cb.dataset.fam); });
    if (picked.length === 0 || picked.length === boxes.length) return [];
    return picked;
  }

  function rgbToHex(rgb) {
    if (typeof rgb === "string") return rgb;
    function h(v) {
      var n = Math.max(0, Math.min(255, Math.round(v * 255)));
      return ("0" + n.toString(16)).slice(-2);
    }
    return "#" + h(rgb[0]) + h(rgb[1]) + h(rgb[2]);
  }

  // ---------------------------------------------------------------- svg
  function svg() { return el.figure.querySelector("svg"); }
  function svgSelectAll(sel) {
    var s = svg();
    return s ? Array.prototype.slice.call(s.querySelectorAll(sel)) : [];
  }

  function loadSvg() {
    fetch(cfg.svgUrl, { cache: "no-store" })
      .then(function (r) { return r.text(); })
      .then(function (txt) {
        el.figure.innerHTML = txt;
        var s = svg();
        if (!s) return;
        s.removeAttribute("width");
        s.removeAttribute("height");
        s.style.width = (100 * zoom) + "%";
        s.style.height = "auto";
        Object.keys(colorOverrides).forEach(function (fam) {
          recolourFamily(fam, colorOverrides[fam], true);
        });
        makeDraggable();
      });
  }

  function toggleFamily(fam, on) {
    var g = slug(fam);
    ["gfv-fam-" + g, "gfv-fam-" + g + "-tips"].forEach(function (id) {
      var node = svgById(id);
      if (node) node.style.display = on ? "" : "none";
    });
  }

  function svgById(id) {
    var s = svg();
    return s ? s.getElementById(id) : null;
  }

  function recolourFamily(fam, hex, quiet) {
    if (!quiet) colorOverrides[fam] = hex;
    var g = svgById("gfv-fam-" + slug(fam));
    if (g) {
      g.querySelectorAll("path, line").forEach(function (n) {
        n.setAttribute("stroke", hex);
      });
    }
    var tips = svgById("gfv-fam-" + slug(fam) + "-tips");
    if (tips) {
      tips.querySelectorAll("path").forEach(function (n) {
        n.setAttribute("fill", hex);
        n.setAttribute("stroke", hex);
      });
    }
  }

  // ---------------------------------------------------------------- drag
  function makeDraggable() {
    var targets = svgSelectAll('#gfv-legend, [id^="gfv-label-"], [id^="gfv-title"]');
    targets.forEach(function (node) {
      var offset = { x: 0, y: 0 };
      var start = null;
      node.style.cursor = "move";
      node.addEventListener("pointerdown", function (e) {
        start = { x: e.clientX, y: e.clientY, ox: offset.x, oy: offset.y };
        node.setPointerCapture(e.pointerId);
        e.preventDefault();
      });
      node.addEventListener("pointermove", function (e) {
        if (!start) return;
        var scale = svg().getBoundingClientRect().width / svg().viewBox.baseVal.width;
        offset.x = start.ox + (e.clientX - start.x) / scale;
        offset.y = start.oy + (e.clientY - start.y) / scale;
        node.setAttribute("transform", "translate(" + offset.x + "," + offset.y + ")");
      });
      node.addEventListener("pointerup", function () { start = null; });
    });
  }

  // ---------------------------------------------------------------- zoom
  function setZoom(z) {
    zoom = Math.max(0.4, Math.min(6, z));
    var s = svg();
    if (s) s.style.width = (100 * zoom) + "%";
  }
  el.zoomIn.addEventListener("click", function () { setZoom(zoom * 1.25); });
  el.zoomOut.addEventListener("click", function () { setZoom(zoom / 1.25); });
  el.zoomReset.addEventListener("click", function () { setZoom(1); });

  function applyLabelSize() {
    var px = parseFloat(el.labelSize.value) * (96 / 72);
    svgSelectAll('[id^="gfv-label-"] text, [id^="gfv-label-"]').forEach(function (n) {
      if (n.tagName === "text") n.style.fontSize = px + "px";
    });
  }

  // ---------------------------------------------------------------- regenerate
  el.apply.addEventListener("click", function () {
    el.apply.disabled = true;
    el.applyMsg.textContent = "Regenerating…";
    var payload = {
      label_mode: el.labelMode.value,
      label_font_size: parseFloat(el.labelSize.value),
      legend_font_size: parseFloat(el.legendSize.value),
      label_font_family: el.font.value,
      legend_font_family: el.font.value,
      chrom_label_font_family: el.font.value,
      title_font_family: el.font.value,
      tick_style: el.tickStyle.value,
      orientation: el.orientation.value,
      show_centromeres: el.centromeres.checked,
      show_unplaced: el.showUnplaced.checked,
      show_titles: el.showTitles ? el.showTitles.checked : true,
      only_chromosomes: collectChroms(),
      only_families: collectFamilies(),
      _colors: {}
    };
    if (el.legendLoc.value === "none") {
      payload.legend_show = false;
    } else {
      payload.legend_show = true;
      payload.legend_location = el.legendLoc.value;
    }
    Object.keys(colorOverrides).forEach(function (fam) {
      payload._colors[fam] = hexToRgb(colorOverrides[fam]);
    });

    fetch(cfg.restyleUrl, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload)
    })
      .then(function (r) { return r.json().then(function (b) { return { ok: r.ok, b: b }; }); })
      .then(function (res) {
        if (!res.ok) throw new Error(res.b.error || "Regenerate failed.");
        el.applyMsg.textContent = "Updated.";
        loadSvg();
        refreshChartImages();
      })
      .catch(function (e) { el.applyMsg.textContent = e.message; })
      .finally(function () { el.apply.disabled = false; });
  });

  el.reset.addEventListener("click", function () {
    colorOverrides = {};
    window.location.reload();
  });

  function hexToRgb(hex) {
    var m = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return m
      ? [parseInt(m[1], 16) / 255, parseInt(m[2], 16) / 255, parseInt(m[3], 16) / 255]
      : [0, 0, 0];
  }

  // ---------------------------------------------------------------- downloads
  function chartUrl(which, fmt) {
    var u = cfg.chartBase.replace("FMT", fmt);
    return which && which !== "genes_per_family" ? u + "?name=" + which : u;
  }

  function wireDownloads() {
    root.querySelectorAll(".dl-row [data-fmt]").forEach(function (a) {
      a.href = cfg.downloadBase.replace("FMT", a.dataset.fmt);
    });
    if (cfg.bundleUrl) {
      if (el.dlBundle) el.dlBundle.href = cfg.bundleUrl;
      if (el.dlBundleAll) el.dlBundleAll.href = cfg.bundleUrl + "?all=1";
    }
    el.dlEdited.addEventListener("click", function (ev) {
      ev.preventDefault();
      var s = svg();
      if (!s) return;
      var blob = new Blob(
        ['<?xml version="1.0" encoding="UTF-8"?>\n', new XMLSerializer().serializeToString(s)],
        { type: "image/svg+xml" }
      );
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url;
      a.download = "gfviewer-edited.svg";
      a.click();
      URL.revokeObjectURL(url);
    });
  }

  // ---------------------------------------------------------------- analytics
  function loadAnalytics(colorMap) {
    var url = cfg.analyticsBase.replace("NAME", "analytics_summary.json");
    fetch(url)
      .then(function (r) { return r.ok ? r.json() : null; })
      .then(function (data) {
        if (!data) return;
        el.analytics.hidden = false;
        el.analyticsSummary.innerHTML =
          renderSummary(data) + renderTelomereTable(data.telomere_bias_table);
        el.coloc.innerHTML = renderColocTable(data.colocalization_table);

        var gpf = data.genes_per_family || [];
        el.genesChart.innerHTML = renderGenesChart(data, colorMap || {});
        el.genesSection.hidden = gpf.length === 0;
        el.genesSection.querySelectorAll("[data-chart-fmt]").forEach(function (a) {
          a.href = cfg.chartBase.replace("FMT", a.dataset.chartFmt);
        });

        wireAnalyticsFigures(data);
        bindRowLimits(el.analytics);
        root.querySelectorAll("[data-name]").forEach(function (a) {
          a.href = cfg.analyticsBase.replace("NAME", a.dataset.name);
        });
      });
  }

  // after a Regenerate the server discards cached chart images so they pick up
  // the new style (e.g. titles on/off); repoint the previews with a cache-buster
  function refreshChartImages() {
    if (!el.analyticsFigs) return;
    var stamp = "_v=" + Date.now();
    el.analyticsFigs.querySelectorAll(".an-fig").forEach(function (block) {
      if (block.hidden) return;
      var img = block.querySelector("[data-chart-img]");
      if (!img) return;
      var u = chartUrl(block.dataset.chart, "png");
      img.src = u + (u.indexOf("?") < 0 ? "?" : "&") + stamp;
    });
  }

  // preview + download links for the stand-alone analytics figures
  function wireAnalyticsFigures(data) {
    if (!el.analyticsFigs) return;
    var have = {
      positional_profile: true,   // always produced when analytics run
      ripley: (data.ripley_clustering || []).length > 0,
      family_proximity: (data.family_proximity_order || []).length >= 2
    };
    var any = false;
    el.analyticsFigs.querySelectorAll(".an-fig").forEach(function (block) {
      var which = block.dataset.chart;
      if (!have[which]) { block.hidden = true; return; }
      block.hidden = false;
      any = true;
      var img = block.querySelector("[data-chart-img]");
      if (img && !img.getAttribute("src")) img.src = chartUrl(which, "png");
      block.querySelectorAll("[data-chart-fmt]").forEach(function (a) {
        a.href = chartUrl(which, a.dataset.chartFmt);
      });
    });
    el.analyticsFigs.hidden = !any;
  }

  // Cap each analytics table to N visible rows (default 10) with a scroll
  // region; the per-table <select> lets the user change N or show all.
  function bindRowLimits(root) {
    root.querySelectorAll(".table-block").forEach(function (block) {
      var scroll = block.querySelector(".table-scroll");
      var sel = block.querySelector(".rows-select");
      if (!scroll || !sel) return;
      function apply() {
        var n = parseInt(sel.value, 10);
        var table = scroll.querySelector("table");
        var rows = table ? table.querySelectorAll("tbody tr, tr") : [];
        if (!n || rows.length <= n) {
          scroll.style.maxHeight = "none";
          return;
        }
        var headH = table.tHead ? table.tHead.getBoundingClientRect().height
          : (rows[0] ? rows[0].getBoundingClientRect().height : 26);
        var rowH = 0, counted = 0;
        rows.forEach(function (r) {
          if (r.parentElement && r.parentElement.tagName === "THEAD") return;
          if (counted < 3) { rowH += r.getBoundingClientRect().height; counted++; }
        });
        rowH = counted ? rowH / counted : 26;
        scroll.style.maxHeight = Math.ceil(headH + n * rowH + 4) + "px";
      }
      sel.addEventListener("change", apply);
      apply();
    });
  }

  function rowLimitControl() {
    return '<label class="rows-ctl">Show ' +
      '<select class="rows-select">' +
      '<option value="10" selected>10 rows</option>' +
      '<option value="25">25 rows</option>' +
      '<option value="50">50 rows</option>' +
      '<option value="0">all rows</option>' +
      "</select></label>";
  }

  function renderSummary(d) {
    function stat(label, value) {
      return '<div class="stat"><span class="k">' + label + '</span><span class="v">' +
        value + "</span></div>";
    }
    var html = "";
    html += stat("Genes", d.n_genes);
    html += stat("Families", d.n_families);
    html += stat("Chromosomes", d.n_chromosomes);
    html += stat("Sub-telomeric", (100 * d.frac_subtelomeric).toFixed(1) + "%");
    if (d.largest_family) html += stat("Largest family", d.largest_family.gene_family);
    if (d.most_clustered_family) {
      html += stat("Most clustered",
        d.most_clustered_family.gene_family + " (" +
        (100 * d.most_clustered_family.frac_clustered).toFixed(0) + "%)");
    }
    if (d.n_genes_on_unplaced) {
      html += stat("On unplaced contigs",
        d.n_genes_on_unplaced + " gene" + (d.n_genes_on_unplaced === 1 ? "" : "s") +
        " / " + d.n_unplaced_contigs + " contig" +
        (d.n_unplaced_contigs === 1 ? "" : "s"));
    }
    return html;
  }

  function renderGenesChart(d, colorMap) {
    var rows = d.genes_per_family || [];
    if (!rows.length) return "";
    var max = rows.reduce(function (m, r) { return Math.max(m, r.total); }, 1);
    var anyUnplaced = rows.some(function (r) { return r.on_unplaced > 0; });
    var bars = rows.map(function (r) {
      var col = rgbToHex(colorMap[r.gene_family] || [0.45, 0.12, 0.17]);
      var wc = (100 * r.on_chromosomes / max).toFixed(2);
      var wu = (100 * r.on_unplaced / max).toFixed(2);
      var seg =
        '<span class="gpf-seg chrom" style="width:' + wc + '%;background:' + col + '"></span>' +
        (r.on_unplaced
          ? '<span class="gpf-seg unplaced" style="width:' + wu + '%"></span>'
          : "");
      var count = r.total + (r.on_unplaced
        ? " (" + r.on_chromosomes + " + " + r.on_unplaced + ")" : "");
      return '<div class="gpf-row"><span class="gpf-label" title="' +
        escapeHtml(r.gene_family) + '">' + escapeHtml(r.gene_family) + '</span>' +
        '<span class="gpf-bar">' + seg + '</span>' +
        '<span class="gpf-count">' + count + "</span></div>";
    }).join("");
    var legend = anyUnplaced
      ? '<p class="muted"><span class="gpf-key chrom"></span> on chromosomes ' +
        '&nbsp;&nbsp;<span class="gpf-key unplaced"></span> on unplaced / stray contigs</p>'
      : '<p class="muted">Bar length is the number of genes assigned to each family.</p>';
    return legend + '<div class="gpf">' + bars + "</div>";
  }

  function renderColocTable(rows) {
    if (!rows || !rows.length) {
      return '<h4>Pairwise co-localization</h4><p class="muted">Not computed for ' +
        "this run (enable it on the form &mdash; it is O(families&sup2;) and slow).</p>";
    }
    var head = "<tr><th>Family A</th><th>Family B</th><th>observed<br>close genes</th>" +
      "<th>null mean</th><th>p-value</th><th>sig.</th></tr>";
    var body = rows.map(function (r) {
      return "<tr" + (r.significant ? ' class="sig"' : "") + ">" +
        "<td>" + escapeHtml(r.family_a) + "</td>" +
        "<td>" + escapeHtml(r.family_b) + "</td>" +
        "<td>" + r.observed_close_pairs + "</td>" +
        "<td>" + r.null_mean.toFixed(2) + "</td>" +
        "<td>" + r.p_value.toFixed(3) + "</td>" +
        "<td>" + (r.significant ? "✓" : "") + "</td></tr>";
    }).join("");
    var ctl = rows.length > 10 ? rowLimitControl() : "";
    return '<div class="table-block">' +
      '<div class="table-head"><h4>Pairwise co-localization</h4>' + ctl + "</div>" +
      '<div class="table-scroll"><table class="analytics-table">' +
      "<thead>" + head + "</thead><tbody>" + body + "</tbody></table></div>" +
      '<p class="muted">Ordered by p-value. For each family pair, the number of ' +
      'genes of A within the co-localization window of a gene of B, compared ' +
      'with a permutation null (A re-placed at random on the same chromosomes). ' +
      'Small p &rArr; the two families sit together more than expected.</p></div>';
  }

  function renderTelomereTable(rows) {
    if (!rows || !rows.length) return "";
    var head =
      "<tr><th>Gene family</th><th>n genes</th><th>obs. mean<br>norm. dist.</th>" +
      "<th>null mean</th><th>direction</th><th>p (&rarr; telomere)</th>" +
      "<th>p (&rarr; interior)</th><th>sig.</th></tr>";
    var body = rows.map(function (r) {
      var p = Math.min(r.p_toward_telomere, r.p_toward_interior);
      return "<tr" + (r.significant ? ' class="sig"' : "") + ">" +
        "<td>" + escapeHtml(r.gene_family) + "</td>" +
        "<td>" + r.n_genes + "</td>" +
        "<td>" + r.observed_mean_norm_dist.toFixed(3) + "</td>" +
        "<td>" + r.null_mean_norm_dist.toFixed(3) + "</td>" +
        "<td>" + escapeHtml(r.direction) + "</td>" +
        "<td>" + r.p_toward_telomere.toFixed(3) + "</td>" +
        "<td>" + r.p_toward_interior.toFixed(3) + "</td>" +
        "<td>" + (r.significant ? "✓" : "") + "</td></tr>";
    }).join("");
    var ctl = rows.length > 10 ? rowLimitControl() : "";
    return '<div class="table-block">' +
      '<div class="table-head"><h4>Telomere-proximity test (per gene family)</h4>' +
      ctl + "</div>" +
      '<div class="table-scroll"><table class="analytics-table">' +
      "<thead>" + head + "</thead><tbody>" + body + "</tbody></table></div>" +
      '<p class="muted">Rows are ordered by significance (most telomere-proximal ' +
      'first). Permutation test vs. random placement on the same chromosomes; ' +
      'rows with min(p) &le; 0.05 are highlighted.</p></div>';
  }

  poll();
})();
