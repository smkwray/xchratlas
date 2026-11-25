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
  if (!Number.isFinite(p)) return '\u2014';
  return p === 0
    ? '<abbr class="pval-bound" title="Reported as 0 at the precision limit">0 (reported)</abbr>'
    : p.toExponential(2);
}

/* ── State ───────────────────────────────────── */
var geneIndex = {};
var geneIds = [];
var eqtlOnly = false;

/* ── Init ────────────────────────────────────── */
function init() {
  fetchJSON('gene_index.json').then(function (data) {
    geneIndex = data;
    geneIds = Object.keys(data);

    eqtlOnly = document.getElementById('eqtl-only-toggle').checked;
    var total = geneIds.length;
    var withHit = 0, withSupport = 0;
    geneIds.forEach(function (name) {
      var g = data[name];
      if (g.returned_trait_count) withHit++;
      if (g.threshold_trait_count) withSupport++;
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
    if (preGene && !input.value) {
      input.value = preGene;
    }
    doSearch(input.value);
  }).catch(function (err) {
    console.error('Gene index load failed:', err);
    document.getElementById('gene-results').innerHTML =
      '<div class="gene-empty">Could not load gene search. <a href="index.html">Return to home</a>.</div>';
  });
}

/* ── Search ──────────────────────────────────── */
function doSearch(query) {
  var container = document.getElementById('gene-results');
  var counter = document.getElementById('gene-result-count');
  var q = query.trim().toLowerCase();

  if (q.length < 1) {
    container.innerHTML = '<div class="gene-empty">Search ' + fmt(geneIds.length) + ' indexed chromosome X genes by symbol or Ensembl ID.</div>';
    counter.textContent = '';
    return;
  }

  function matchRank(id) {
    var keys = [id.toLowerCase(), geneIndex[id].gene_name.toLowerCase()];
    if (keys.includes(q)) return 0;
    if (keys.some(function (key) { return key.startsWith(q); })) return 1;
    return 2;
  }
  var matches = geneIds.filter(function (id) {
    var gene = geneIndex[id];
    return (!eqtlOnly || gene.returned_trait_count) &&
      (id.toLowerCase().includes(q) || gene.gene_name.toLowerCase().includes(q));
  }).sort(function (a, b) {
    return matchRank(a) - matchRank(b) ||
      geneIndex[b].threshold_trait_count - geneIndex[a].threshold_trait_count ||
      geneIndex[a].gene_name.localeCompare(geneIndex[b].gene_name) || a.localeCompare(b);
  });

  var capped = matches.slice(0, 50);
  counter.textContent = matches.length + ' gene' + (matches.length !== 1 ? 's' : '') + ' found' +
    (matches.length > 50 ? ' (showing first 50)' : '');

  if (matches.length === 0) {
    container.innerHTML = '<div class="gene-empty">No indexed genes match \u201c' + escHtml(query.trim()) + '\u201d' +
      (eqtlOnly ? ' with eQTL records' : '') + '.</div>';
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
function groupAppearances(hits) {
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

  return groups;
}

function renderGeneCard(id, query) {
  var gene = geneIndex[id], hits = gene.hits;
  var groups = groupAppearances(hits);
  var biotype = gene.biotype || 'unknown';
  var biotypeCls = biotype === 'protein_coding' ? 'biotype-protein' : 'biotype-other';
  var threshold = gene.threshold_trait_count + ' ' + (gene.threshold_trait_count === 1 ? 'trait' : 'traits') + ' \u00b7 p \u2264 10\u207b\u2075';
  var returned = gene.returned_trait_count + ' ' + (gene.returned_trait_count === 1 ? 'trait' : 'traits') + ' \u00b7 eQTL records';
  var showCount = 5;
  var visible = groups.slice(0, showCount);
  var overflow = groups.length > showCount ? groups.length - showCount : 0;

  return '<div class="gene-card">' +
    '<div class="gene-card-header">' +
      '<div class="gene-card-title-row">' +
        '<span class="gene-card-name">' + highlightName(gene.gene_name, query) + '</span>' +
        '<span class="biotype-badge ' + biotypeCls + '">' + biotype.replace(/_/g, ' ') + '</span>' +
        (gene.gene_id ? '<span class="gene-card-id">' + gene.gene_id + '</span>' : '') +
      '</div>' +
      '<div class="gene-card-summary">' +
        '<span class="gene-card-fact">' + threshold + '</span>' +
        '<span class="gene-card-fact">' + returned + '</span>' +
        '<span class="gene-card-fact">' + hits.length + ' panel\u2013trait\u2013locus ' + (hits.length === 1 ? 'entry' : 'entries') + '</span>' +
      '</div>' +
    '</div>' +
    '<div class="gene-card-traits" tabindex="0" role="region" aria-label="Gene traits">' +
      '<table class="gene-trait-table">' +
        '<thead><tr>' +
          '<th>Trait</th>' +
          '<th>Panel</th>' +
          '<th>Category</th>' +
          '<th>eQTL</th>' +
          '<th>Best eQTL p</th>' +
          '<th>Loci</th>' +
        '</tr></thead>' +
        '<tbody>' +
          visible.map(function (g) { return traitRow(g); }).join('') +
          (overflow > 0 ? '<tr class="gene-overflow-row" data-gene="' + escHtml(id) + '"><td colspan="6" class="gene-overflow-cell">' +
            '<button class="gene-show-more" onclick="expandGene(this, \'' + escHtml(id).replace(/'/g, '\\\'') + '\')">Show ' + overflow + ' more panel ' + (overflow === 1 ? 'entry' : 'entries') + '</button>' +
          '</td></tr>' : '') +
        '</tbody>' +
      '</table>' +
    '</div>' +
  '</div>';
}

function traitRow(g) {
  var statusCls, statusLabel;
  if (g.best_supported) {
    statusCls = 'is-supported';
    var role = g.loci.some(function (h) { return h.eqtl_supported && h.eqtl_gene_role !== 'followup'; })
      ? 'Candidate' : 'Follow-up';
    statusLabel = role + ' \u00b7 p \u2264 10\u207b\u2075';
  } else if (g.best_hit) {
    statusCls = 'is-lookup-hit';
    statusLabel = 'Returned records';
  } else {
    statusCls = 'is-unsupported';
    statusLabel = g.loci.every(function (h) { return h.eqtl_lookup_status === 'skipped_no_rsid'; })
      ? 'Not assessed' : 'No records';
  }

  var roles = new Set(g.loci.map(function (h) {
    return h.eqtl_gene_role === 'followup' ? 'Follow-up' : 'Candidate';
  }));
  var href = 'trait.html?panel=' + g.panel_id + '&trait=' + encodeURIComponent(g.trait_id);

  return '<tr class="gene-trait-row">' +
    '<td class="gene-trait-name"><a href="' + href + '">' + escHtml(g.trait_desc) + '</a></td>' +
    '<td class="gene-trait-panel">' + escHtml(g.panel_title) + '</td>' +
    '<td><span class="domain-badge">' + fmtDomain(g.domain) + '</span></td>' +
    '<td><span class="discovery-trait-status ' + statusCls + '">' + statusLabel + '</span>' +
      (!g.best_supported ? '<div class="gene-role-sub">' + Array.from(roles).join(' / ') + '</div>' : '') + '</td>' +
    '<td class="num-cell">' + fmtPval(g.best_pval) + '</td>' +
    '<td class="num-cell">' + g.loci.length + '</td>' +
  '</tr>';
}

/* ── Expand overflow ─────────────────────────── */
function expandGene(btn, geneName) {
  var gene = geneIndex[geneName];
  if (!gene) return;
  var groups = groupAppearances(gene.hits);
  var extra = groups.slice(5);
  var overflowRow = btn.closest('tr');
  overflowRow.outerHTML = extra.map(function (g) { return traitRow(g); }).join('');
}

document.addEventListener('DOMContentLoaded', init);
