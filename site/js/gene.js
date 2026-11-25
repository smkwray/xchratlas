'use strict';

var CONFIG = { dataRoot: 'data' };

function fetchJSON(p) {
  return fetch(CONFIG.dataRoot + '/' + p).then(function (r) {
    if (!r.ok) throw new Error('Fetch failed: ' + p); return r.json();
  });
}
function fmt(n) { return Number(n).toLocaleString(); }
function fmtDomain(d) {
  if (d === null || d === undefined || d === '') return 'Unlabeled';
  return d.replace(/_/g, ' ').replace(/\b\w/g, function (c) { return c.toUpperCase(); });
}
function fmtPval(p) {
  if (p === null || p === undefined) return '\u2014';
  if (p < 1e-100) return '<1e-100';
  return p.toExponential(2);
}

/* ── State ───────────────────────────────────── */
var geneIndex = {};
var geneNames = [];
var eqtlOnly = false;

/* ── Init ────────────────────────────────────── */
function init() {
  fetchJSON('gene_index.json').then(function (data) {
    geneIndex = data;
    geneNames = Object.keys(data);

    var total = geneNames.length;
    var withHit = 0, withSupport = 0;
    geneNames.forEach(function (name) {
      var g = data[name];
      if (g.hits.some(function (h) { return h.eqtl_lookup_hit; })) withHit++;
      if (g.hits.some(function (h) { return h.eqtl_supported; })) withSupport++;
    });

    document.getElementById('gm-total').textContent = fmt(total);
    document.getElementById('gm-eqtl').textContent = fmt(withHit);
    document.getElementById('gm-supported').textContent = fmt(withSupport);

    var input = document.getElementById('gene-search');
    var toggle = document.getElementById('eqtl-only-toggle');

    input.addEventListener('input', function () { doSearch(input.value); });
    toggle.addEventListener('change', function () {
      eqtlOnly = this.checked;
      doSearch(input.value);
    });

    // Show initial state or pre-fill from URL
    var params = new URLSearchParams(location.search);
    var preGene = params.get('gene');
    if (preGene) {
      input.value = preGene;
    }
    doSearch(input.value);
  }).catch(function (err) {
    console.error('Gene index load failed:', err);
    document.getElementById('gene-results').innerHTML =
      '<div class="gene-empty">Failed to load gene index. Ensure gene_index.json exists in data/.</div>';
  });
}

/* ── Search ──────────────────────────────────── */
function doSearch(query) {
  var container = document.getElementById('gene-results');
  var counter = document.getElementById('gene-result-count');
  var q = query.trim().toLowerCase();

  if (q.length < 1) {
    container.innerHTML = '<div class="gene-empty">Start typing to search ' + fmt(geneNames.length) + ' chromosome X genes.</div>';
    counter.textContent = '';
    return;
  }

  var matches = geneNames.filter(function (name) {
    if (eqtlOnly && !geneIndex[name].hits.some(function (h) { return h.eqtl_lookup_hit; })) return false;
    return name.toLowerCase().includes(q);
  });

  // Sort: exact matches first, then prefix, then contains
  matches.sort(function (a, b) {
    var al = a.toLowerCase(), bl = b.toLowerCase();
    var aExact = al === q, bExact = bl === q;
    if (aExact !== bExact) return aExact ? -1 : 1;
    var aPrefix = al.indexOf(q) === 0, bPrefix = bl.indexOf(q) === 0;
    if (aPrefix !== bPrefix) return aPrefix ? -1 : 1;
    // Then by supported hit count
    var aSupp = geneIndex[a].hits.filter(function (h) { return h.eqtl_supported; }).length;
    var bSupp = geneIndex[b].hits.filter(function (h) { return h.eqtl_supported; }).length;
    if (aSupp !== bSupp) return bSupp - aSupp;
    return a.localeCompare(b);
  });

  var capped = matches.slice(0, 50);
  counter.textContent = matches.length + ' gene' + (matches.length !== 1 ? 's' : '') + ' found' +
    (matches.length > 50 ? ' (showing first 50)' : '');

  if (matches.length === 0) {
    container.innerHTML = '<div class="gene-empty">No genes match \u201c' + escHtml(query.trim()) + '\u201d' +
      (eqtlOnly ? ' with eQTL hits' : '') + '.</div>';
    return;
  }

  container.innerHTML = capped.map(function (name) { return renderGeneCard(name, q); }).join('');
}

function escHtml(s) {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}

function highlightName(name, query) {
  var idx = name.toLowerCase().indexOf(query);
  if (idx < 0) return escHtml(name);
  return escHtml(name.substring(0, idx)) +
    '<mark class="gene-highlight">' + escHtml(name.substring(idx, idx + query.length)) + '</mark>' +
    escHtml(name.substring(idx + query.length));
}

/* ── Render Gene Card ────────────────────────── */
function renderGeneCard(name, query) {
  var gene = geneIndex[name];
  var hits = gene.hits;
  var hitCount = hits.filter(function (h) { return h.eqtl_lookup_hit; }).length;
  var supportCount = hits.filter(function (h) { return h.eqtl_supported; }).length;
  var biotype = gene.biotype || 'unknown';
  var biotypeCls = biotype === 'protein_coding' ? 'biotype-protein' : 'biotype-other';

  // Group hits by trait (across panels)
  var traitGroups = {};
  hits.forEach(function (h) {
    var key = h.trait_id + '|' + h.panel_id;
    if (!traitGroups[key]) {
      traitGroups[key] = {
        trait_id: h.trait_id,
        trait_desc: h.trait_desc,
        panel_id: h.panel_id,
        panel_title: h.panel_title,
        domain: h.domain,
        best_supported: false,
        best_hit: false,
        best_pval: null,
        loci: [],
      };
    }
    var g = traitGroups[key];
    if (h.eqtl_supported) g.best_supported = true;
    if (h.eqtl_lookup_hit) g.best_hit = true;
    if (h.best_eqtl_pvalue !== null && (g.best_pval === null || h.best_eqtl_pvalue < g.best_pval)) {
      g.best_pval = h.best_eqtl_pvalue;
    }
    g.loci.push(h);
  });

  var groups = Object.keys(traitGroups).map(function (k) { return traitGroups[k]; });
  groups.sort(function (a, b) {
    if (a.best_supported !== b.best_supported) return a.best_supported ? -1 : 1;
    if (a.best_hit !== b.best_hit) return a.best_hit ? -1 : 1;
    if (a.best_pval !== b.best_pval) {
      if (a.best_pval === null) return 1;
      if (b.best_pval === null) return -1;
      return a.best_pval - b.best_pval;
    }
    return a.trait_desc.localeCompare(b.trait_desc);
  });

  var statusCls, statusLabel;
  if (supportCount > 0) {
    statusCls = 'is-supported'; statusLabel = supportCount + ' trait' + (supportCount > 1 ? 's' : '') + ' eQTL-supported';
  } else if (hitCount > 0) {
    statusCls = 'is-lookup-hit'; statusLabel = hitCount + ' trait' + (hitCount > 1 ? 's' : '') + ' with lookup hits';
  } else {
    statusCls = 'is-unsupported'; statusLabel = 'No eQTL hits';
  }

  var showCount = 5;
  var visible = groups.slice(0, showCount);
  var overflow = groups.length > showCount ? groups.length - showCount : 0;

  return '<div class="gene-card">' +
    '<div class="gene-card-header">' +
      '<div class="gene-card-title-row">' +
        '<span class="gene-card-name">' + highlightName(name, query) + '</span>' +
        '<span class="biotype-badge ' + biotypeCls + '">' + biotype.replace(/_/g, ' ') + '</span>' +
        (gene.gene_id ? '<span class="gene-card-id">' + gene.gene_id + '</span>' : '') +
      '</div>' +
      '<div class="gene-card-summary">' +
        '<span class="discovery-trait-status ' + statusCls + '">' + statusLabel + '</span>' +
        '<span class="gene-card-fact">' + hits.length + ' trait\u2013locus appearance' + (hits.length > 1 ? 's' : '') + '</span>' +
      '</div>' +
    '</div>' +
    '<div class="gene-card-traits">' +
      '<table class="gene-trait-table">' +
        '<thead><tr>' +
          '<th>Trait</th>' +
          '<th>Panel</th>' +
          '<th>Domain</th>' +
          '<th>eQTL Status</th>' +
          '<th>Best p-value</th>' +
          '<th>Loci</th>' +
        '</tr></thead>' +
        '<tbody>' +
          visible.map(function (g) { return traitRow(g); }).join('') +
          (overflow > 0 ? '<tr class="gene-overflow-row" data-gene="' + escHtml(name) + '"><td colspan="6" class="gene-overflow-cell">' +
            '<button class="gene-show-more" onclick="expandGene(this, \'' + escHtml(name).replace(/'/g, '\\\'') + '\')">Show ' + overflow + ' more trait' + (overflow > 1 ? 's' : '') + '</button>' +
          '</td></tr>' : '') +
        '</tbody>' +
      '</table>' +
    '</div>' +
  '</div>';
}

function traitRow(g) {
  var statusCls, statusLabel;
  if (g.best_supported) { statusCls = 'is-supported'; statusLabel = 'Supported'; }
  else if (g.best_hit) { statusCls = 'is-lookup-hit'; statusLabel = 'Lookup-hit'; }
  else { statusCls = 'is-unsupported'; statusLabel = 'None'; }

  var href = 'trait.html?panel=' + g.panel_id + '&trait=' + encodeURIComponent(g.trait_id);

  return '<tr class="gene-trait-row" onclick="location.href=\'' + href + '\'">' +
    '<td class="gene-trait-name">' + escHtml(g.trait_desc) + '</td>' +
    '<td class="gene-trait-panel">' + escHtml(g.panel_title) + '</td>' +
    '<td><span class="domain-badge">' + fmtDomain(g.domain) + '</span></td>' +
    '<td><span class="discovery-trait-status ' + statusCls + '">' + statusLabel + '</span></td>' +
    '<td class="num-cell">' + fmtPval(g.best_pval) + '</td>' +
    '<td class="num-cell">' + g.loci.length + '</td>' +
  '</tr>';
}

/* ── Expand overflow ─────────────────────────── */
function expandGene(btn, geneName) {
  var gene = geneIndex[geneName];
  if (!gene) return;
  var traitGroups = {};
  gene.hits.forEach(function (h) {
    var key = h.trait_id + '|' + h.panel_id;
    if (!traitGroups[key]) {
      traitGroups[key] = {
        trait_id: h.trait_id, trait_desc: h.trait_desc,
        panel_id: h.panel_id, panel_title: h.panel_title,
        domain: h.domain, best_supported: false, best_hit: false, best_pval: null, loci: [],
      };
    }
    var g = traitGroups[key];
    if (h.eqtl_supported) g.best_supported = true;
    if (h.eqtl_lookup_hit) g.best_hit = true;
    if (h.best_eqtl_pvalue !== null && (g.best_pval === null || h.best_eqtl_pvalue < g.best_pval)) g.best_pval = h.best_eqtl_pvalue;
    g.loci.push(h);
  });
  var groups = Object.keys(traitGroups).map(function (k) { return traitGroups[k]; });
  groups.sort(function (a, b) {
    if (a.best_supported !== b.best_supported) return a.best_supported ? -1 : 1;
    if (a.best_hit !== b.best_hit) return a.best_hit ? -1 : 1;
    if (a.best_pval !== b.best_pval) {
      if (a.best_pval === null) return 1;
      if (b.best_pval === null) return -1;
      return a.best_pval - b.best_pval;
    }
    return a.trait_desc.localeCompare(b.trait_desc);
  });
  var extra = groups.slice(5);
  var overflowRow = btn.closest('tr');
  overflowRow.outerHTML = extra.map(function (g) { return traitRow(g); }).join('');
}

document.addEventListener('DOMContentLoaded', init);
