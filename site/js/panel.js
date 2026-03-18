'use strict';

var CONFIG = { dataRoot: 'data' };
var ACCENTS = {
  panel_expanded: 'teal', panel_mind_risk_core: 'amber',
  panel_blood_biochemistry_core: 'blue', panel_independent_set: 'slate'
};

function fetchJSON(p) {
  return fetch(CONFIG.dataRoot + '/' + p).then(function (r) {
    if (!r.ok) throw new Error('Fetch failed: ' + p); return r.json();
  });
}
function fmt(n) { return Number(n).toLocaleString(); }
function fmtDomain(d) { return d.replace(/_/g, ' ').replace(/\b\w/g, function (c) { return c.toUpperCase(); }); }
function pct(n, d) { return d ? Math.round(n / d * 100) : 0; }

var S = {
  panelId: null, panel: null, summary: null,
  traits: [], filtered: [],
  sortKey: 'x_evidence_score', sortDir: 'desc',
  searchTerm: '', domains: new Set(), eqtlOnly: false
};

function getParam(k) { return new URLSearchParams(location.search).get(k); }

/* ── Hero ──────────────────────────────────── */
function renderHero() {
  var p = S.panel;
  document.title = p.title + ' \u2014 chrXatlas';
  document.getElementById('panel-title').textContent = p.title;
  document.getElementById('panel-desc').textContent = p.description;
  document.getElementById('panel-breadcrumb-name').textContent = p.title;
  document.getElementById('ph-traits').textContent = p.trait_count;
  document.getElementById('ph-loci').textContent = fmt(p.locus_count);
  document.getElementById('ph-supported').textContent = p.supported_trait_count;
  document.getElementById('ph-rate').textContent = pct(p.supported_trait_count, p.trait_count) + '%';
}

/* ── Domain Filter Chips ───────────────────── */
function renderFilters() {
  var el = document.getElementById('domain-filters');
  if (!el || !S.summary) return;
  var ds = S.summary.domain_summary.slice().sort(function (a, b) { return b.locus_count - a.locus_count; });

  el.innerHTML = '<button class="filter-chip active" data-domain="all">All</button>' +
    ds.map(function (d) {
      return '<button class="filter-chip" data-domain="' + d.domain + '">' +
        fmtDomain(d.domain) + ' <span class="chip-count">' + d.trait_count + '</span></button>';
    }).join('');

  el.addEventListener('click', function (e) {
    var chip = e.target.closest('.filter-chip');
    if (!chip) return;
    var domain = chip.dataset.domain;
    if (domain === 'all') {
      S.domains.clear();
      el.querySelectorAll('.filter-chip').forEach(function (c) { c.classList.remove('active'); });
      chip.classList.add('active');
    } else {
      el.querySelector('[data-domain="all"]').classList.remove('active');
      chip.classList.toggle('active');
      if (S.domains.has(domain)) S.domains.delete(domain); else S.domains.add(domain);
      if (S.domains.size === 0) el.querySelector('[data-domain="all"]').classList.add('active');
    }
    applyFilters();
  });
}

/* ── Filters + Sort ────────────────────────── */
function applyFilters() {
  var q = S.searchTerm.toLowerCase();
  S.filtered = S.traits.filter(function (t) {
    if (q && !(
      t.description.toLowerCase().includes(q) ||
      t.domain.toLowerCase().includes(q) ||
      t.query_id.toLowerCase().includes(q) ||
      (t.top_candidate_genes || []).some(function (g) { return g.toLowerCase().includes(q); })
    )) return false;
    if (S.domains.size > 0 && !S.domains.has(t.domain)) return false;
    if (S.eqtlOnly && !t.eqtl_supported) return false;
    return true;
  });
  sortTraits();
  renderTable();
  var rc = document.getElementById('result-count');
  if (rc) rc.textContent = S.filtered.length + ' of ' + S.traits.length + ' traits';
}

function sortTraits() {
  var k = S.sortKey, d = S.sortDir === 'asc' ? 1 : -1;
  S.filtered.sort(function (a, b) {
    var av = a[k], bv = b[k];
    if (typeof av === 'string') return d * av.localeCompare(bv);
    return d * ((av || 0) - (bv || 0));
  });
}

/* ── Table ─────────────────────────────────── */
function renderTable() {
  var tbody = document.getElementById('traits-tbody');
  if (!tbody) return;
  if (S.filtered.length === 0) {
    tbody.innerHTML = '<tr><td colspan="7" class="table-empty">No traits match the current filters.</td></tr>';
    return;
  }
  tbody.innerHTML = S.filtered.map(function (t) {
    var sw = Math.min(t.x_evidence_score, 100);
    return '<tr class="trait-row" onclick="location.href=\'trait.html?panel=' + S.panelId + '&trait=' + encodeURIComponent(t.trait_id) + '\'">' +
      '<td class="trait-name-cell"><span class="trait-table-name">' + t.description + '</span></td>' +
      '<td><span class="domain-badge">' + fmtDomain(t.domain) + '</span></td>' +
      '<td class="score-cell"><div class="mini-score-bar"><div class="mini-score-fill" style="width:' + sw + '%"></div></div><span class="score-val">' + t.x_evidence_score.toFixed(0) + '</span></td>' +
      '<td class="num-cell">' + t.n_loci + '</td>' +
      '<td class="num-cell">' + fmt(t.eqtl_total_hit_count) + '</td>' +
      '<td><span class="trait-grade grade-' + t.coverage_grade.toLowerCase() + '">' + t.coverage_grade + '</span></td>' +
      '<td class="genes-cell">' + (t.top_candidate_genes || []).slice(0, 2).map(function (g) {
        return '<span class="gene-tag">' + g + '</span>';
      }).join('') + '</td></tr>';
  }).join('');
}

/* ── Sort Headers ──────────────────────────── */
function initSort() {
  document.querySelectorAll('.data-table th[data-sort]').forEach(function (th) {
    th.addEventListener('click', function () {
      var k = this.dataset.sort;
      if (S.sortKey === k) S.sortDir = S.sortDir === 'asc' ? 'desc' : 'asc';
      else { S.sortKey = k; S.sortDir = (k === 'description' || k === 'domain') ? 'asc' : 'desc'; }
      document.querySelectorAll('.data-table th[data-sort]').forEach(function (h) { h.classList.remove('sort-asc', 'sort-desc'); });
      this.classList.add('sort-' + S.sortDir);
      applyFilters();
    });
  });
}

/* ── Init ──────────────────────────────────── */
function init() {
  S.panelId = getParam('id');
  if (!S.panelId) { document.getElementById('panel-title').textContent = 'Panel not specified'; return; }

  fetchJSON('manifest.json').then(function (m) {
    S.panel = m.panels.find(function (p) { return p.panel_id === S.panelId; });
    if (!S.panel) throw new Error('Panel not found');
    renderHero();
    return Promise.all([fetchJSON(S.panel.files.summary), fetchJSON(S.panel.files.traits)]);
  }).then(function (r) {
    S.summary = r[0];
    S.traits = r[1];
    S.filtered = S.traits.slice();
    renderFilters();
    sortTraits();
    renderTable();
    var rc = document.getElementById('result-count');
    if (rc) rc.textContent = S.filtered.length + ' of ' + S.traits.length + ' traits';
    initSort();

    var si = document.getElementById('table-search');
    if (si) si.addEventListener('input', function () { S.searchTerm = this.value; applyFilters(); });
    var et = document.getElementById('eqtl-toggle');
    if (et) et.addEventListener('change', function () { S.eqtlOnly = this.checked; applyFilters(); });
  }).catch(function (e) {
    console.error(e);
    document.getElementById('panel-title').textContent = 'Error: ' + e.message;
  });
}

document.addEventListener('DOMContentLoaded', init);
