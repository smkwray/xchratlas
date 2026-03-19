'use strict';

var CONFIG = { dataRoot: 'data' };

function fetchJSON(p) {
  return fetch(CONFIG.dataRoot + '/' + p).then(function (r) {
    if (!r.ok) throw new Error('Fetch failed: ' + p); return r.json();
  });
}
function fmt(n) { return Number(n).toLocaleString(); }
function fmtDomain(d) { return d.replace(/_/g, ' ').replace(/\b\w/g, function (c) { return c.toUpperCase(); }); }
function fmtPval(v) { if (!v || v === null) return '\u2014'; if (v < 0.001) return v.toExponential(1); return v.toFixed(3); }
function fmtPos(s, e) { return 'X:' + fmt(s) + '\u2013' + fmt(e); }
function getParam(k) { return new URLSearchParams(location.search).get(k); }

/* ── Hero ──────────────────────────────────── */
function renderHero(data) {
  var t = data.trait, p = data.panel;
  document.title = t.description + ' \u2014 ' + p.title + ' \u2014 chrXatlas';
  var bpanel = document.getElementById('trait-breadcrumb-panel');
  bpanel.textContent = p.title;
  bpanel.href = 'panel.html?id=' + p.panel_id;
  document.getElementById('trait-breadcrumb-name').textContent = t.description;
  document.getElementById('trait-title').textContent = t.description;
  document.getElementById('trait-domain').textContent = fmtDomain(t.domain);
  var g = document.getElementById('trait-grade');
  g.textContent = t.coverage_grade;
  g.className = 'trait-grade grade-' + t.coverage_grade.toLowerCase();
  document.getElementById('trait-score-value').textContent = t.x_evidence_score.toFixed(0);
  document.getElementById('trait-score-fill').style.width = Math.min(t.x_evidence_score, 100) + '%';
  document.getElementById('ts-loci').textContent = t.n_loci;
  document.getElementById('ts-eqtl-loci').textContent = t.eqtl_supported_locus_count;
  document.getElementById('ts-eqtl-hits').textContent = fmt(t.eqtl_total_hit_count);
  document.getElementById('ts-pval').textContent = t.lead_neglog10_pvalue.toFixed(1);
}

/* ── Metadata ──────────────────────────────── */
function renderMeta(data) {
  var t = data.trait, el = document.getElementById('metadata-grid');
  if (!el) return;
  var items = [
    ['Phenocode', t.phenocode || '\u2014'],
    ['Trait type', t.trait_type || '\u2014'],
    ['Modifier', t.modifier || 'none'],
    ['Sample size', t.n_cases_hq_cohort_both_sexes ? fmt(t.n_cases_hq_cohort_both_sexes) : '\u2014'],
    ['Populations passing QC', t.num_pops_pass_qc || '\u2014'],
    ['Tier', t.tier || '\u2014'],
    ['Selection rank', t.selection_rank || '\u2014'],
    ['Selection score', t.selection_score || '\u2014'],
    ['Non-PAR loci', t.any_nonpar_locus ? 'Yes' : 'No'],
    ['PAR-only signal', t.par_only_signal ? 'Yes' : 'No'],
    ['Max independent set', t.in_max_independent_set ? 'Yes' : 'No'],
    ['Category', t.category || '\u2014']
  ];
  el.innerHTML = items.map(function (r) {
    return '<div class="metadata-item"><span class="metadata-label">' + r[0] +
      '</span><span class="metadata-value">' + r[1] + '</span></div>';
  }).join('');
}

/* ── Loci ──────────────────────────────────── */
function renderLoci(data) {
  var el = document.getElementById('loci-list');
  if (!el || !data.loci) return;
  document.getElementById('loci-heading').textContent = 'Chromosome X Loci (' + data.loci.length + ')';
  if (data.loci.length === 0) {
    el.innerHTML = '<p class="table-empty">No genome-wide significant loci on chromosome X for this trait.</p>';
    return;
  }
  el.innerHTML = data.loci.map(function (locus, i) {
    var genes = locus.candidate_genes || [];
    return '<details class="locus-card"' + (i < 3 ? ' open' : '') + '>' +
      '<summary class="locus-summary">' +
        '<div class="locus-summary-main">' +
          '<span class="locus-region region-' + locus.x_region.toLowerCase() + '">' + locus.x_region + '</span>' +
          '<span class="locus-position">' + fmtPos(locus.locus_start, locus.locus_end) + '</span>' +
          '<span class="locus-lead-variant">' + (locus.lead_rsid || locus.lead_varid) + '</span>' +
        '</div>' +
        '<div class="locus-summary-stats">' +
          '<span class="locus-pval">p: 10<sup>\u2212' + locus.lead_neglog10_pvalue.toFixed(1) + '</sup></span>' +
          '<span class="locus-top-gene">' + locus.top_gene_name + '</span>' +
          (locus.eqtl_supported
            ? '<span class="locus-eqtl-badge eqtl-yes">eQTL-supported</span>'
            : '<span class="locus-eqtl-badge eqtl-no">no eQTL</span>') +
          '<span class="locus-gene-count">' + genes.length + ' gene' + (genes.length !== 1 ? 's' : '') + '</span>' +
        '</div>' +
      '</summary>' +
      '<div class="locus-body">' +
        (genes.length > 0 ? geneTable(genes) : '<p style="color:var(--text-tertiary);font-size:0.88rem">No candidate genes mapped.</p>') +
      '</div>' +
    '</details>';
  }).join('');
}

function geneTable(genes) {
  return '<table class="candidate-table"><thead><tr>' +
    '<th>Rank</th><th>Gene</th><th>Biotype</th><th>Relation</th><th>Distance</th><th>eQTL</th><th>Studies</th><th>Best p</th>' +
    '</tr></thead><tbody>' +
    genes.map(function (g) {
      return '<tr class="' + (g.eqtl_supported ? 'gene-supported' : '') + '">' +
        '<td class="num-cell">' + g.candidate_rank + '</td>' +
        '<td><span class="gene-name">' + g.gene_name + '</span></td>' +
        '<td><span class="biotype-badge">' + g.gene_biotype.replace(/_/g, ' ') + '</span></td>' +
        '<td>' + g.mapping_relation.replace(/_/g, ' ') + '</td>' +
        '<td class="num-cell">' + fmt(g.distance_to_lead_bp) + ' bp</td>' +
        '<td>' + (g.eqtl_supported ? '<span class="eqtl-yes-sm">Yes</span>' : '\u2014') + '</td>' +
        '<td class="num-cell">' + g.eqtl_study_count + '</td>' +
        '<td class="num-cell">' + fmtPval(g.best_eqtl_pvalue) + '</td></tr>';
    }).join('') + '</tbody></table>';
}

/* ── Support Explanation ───────────────────── */
function fmtMb(bp) {
  if (bp >= 1e6) return (bp / 1e6).toFixed(1) + 'M';
  if (bp >= 1e3) return (bp / 1e3).toFixed(0) + 'K';
  return fmt(bp);
}

function renderSupport(data) {
  var el = document.getElementById('support-card');
  if (!el) return;
  var t = data.trait, loci = data.loci || [];
  var supported = t.eqtl_supported;
  var supLoci = loci.filter(function (l) { return l.eqtl_supported; });

  var statusClass = supported ? 'support-yes' : 'support-no';
  var statusLabel = supported ? 'eQTL-supported' : 'Not eQTL-supported';
  var statusIcon = supported ? '\u2713' : '\u2717';

  var html = '<div class="support-card ' + statusClass + '">' +
    '<div class="support-header">' +
      '<span class="support-status-badge ' + statusClass + '">' + statusIcon + ' ' + statusLabel + '</span>' +
      '<h3>' + (supported
        ? 'Why this trait is marked eQTL-supported'
        : 'Why this trait is not marked eQTL-supported') + '</h3>' +
    '</div>' +
    '<div class="support-body">' +
      '<div class="support-facts">' +
        '<div class="support-fact">' +
          '<span class="support-fact-label">eQTL-supported loci</span>' +
          '<span class="support-fact-value">' + t.eqtl_supported_locus_count + ' of ' + t.n_loci + '</span>' +
        '</div>' +
        '<div class="support-fact">' +
          '<span class="support-fact-label">Total eQTL hits</span>' +
          '<span class="support-fact-value">' + fmt(t.eqtl_total_hit_count) + '</span>' +
        '</div>' +
        '<div class="support-fact">' +
          '<span class="support-fact-label">Coverage grade</span>' +
          '<span class="support-fact-value"><span class="trait-grade grade-' + t.coverage_grade.toLowerCase() + '">' + t.coverage_grade + '</span></span>' +
        '</div>' +
      '</div>';

  if (supLoci.length > 0) {
    html += '<div class="support-loci">' +
      '<h4>Supporting ' + (supLoci.length === 1 ? 'locus' : 'loci') + '</h4>' +
      supLoci.map(function (l) {
        var topGene = (l.candidate_genes || []).find(function (g) { return g.eqtl_supported; });
        return '<div class="support-locus-row">' +
          '<span class="locus-region region-' + l.x_region.toLowerCase() + '">' + l.x_region + '</span>' +
          '<span class="support-locus-range">chrX ' + fmtMb(l.locus_start) + '\u2013' + fmtMb(l.locus_end) + '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Top gene:</strong> ' + l.top_gene_name +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Lead variant:</strong> ' + (l.lead_rsid || l.lead_varid) +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Lookup:</strong> ' + (l.eqtl_lookup_status || 'unknown') +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>eQTL hits:</strong> ' + l.eqtl_lookup_n_hits +
          '</span>' +
          (topGene ? '<span class="support-locus-detail"><strong>eQTL gene:</strong> <span class="gene-name">' + topGene.gene_name + '</span> (' + topGene.gene_biotype.replace(/_/g, ' ') + ', ' + topGene.eqtl_study_count + ' ' + (topGene.eqtl_study_count === 1 ? 'study' : 'studies') + ')</span>' : '') +
        '</div>';
      }).join('') +
    '</div>';
  } else {
    html += '<div class="support-loci">' +
      '<p class="support-none-detail">No chromosome X loci had eQTL-supported candidate-gene evidence in the current build.' +
      (t.eqtl_lookup_note ? ' ' + t.eqtl_lookup_note : '') + '</p>' +
    '</div>';
  }

  html += '<div class="support-rule">' +
      '<strong>Rule:</strong> trait-level eQTL support means at least one chrX locus had eQTL-supported candidate-gene evidence.' +
    '</div>' +
  '</div></div>';

  el.innerHTML = html;
}

/* ── Notes ─────────────────────────────────── */
function renderNotes(data) {
  var t = data.trait, el = document.getElementById('trait-notes-content');
  if (!el) return;
  var notes = [];
  if (t.confidence_notes) notes.push(t.confidence_notes);
  if (t.description_more) notes.push('<strong>Measurement:</strong> ' + t.description_more);
  if (t.eqtl_lookup_note) notes.push('<strong>eQTL lookup:</strong> ' + t.eqtl_lookup_note);
  el.innerHTML = notes.map(function (n) { return '<p>' + n + '</p>'; }).join('');
}

/* ── Init ──────────────────────────────────── */
function init() {
  var panelId = getParam('panel'), traitId = getParam('trait');
  if (!panelId || !traitId) {
    document.getElementById('trait-title').textContent = 'Trait not specified';
    return;
  }

  fetchJSON('manifest.json').then(function (m) {
    var panel = m.panels.find(function (p) { return p.panel_id === panelId; });
    if (!panel) throw new Error('Panel not found: ' + panelId);
    return fetchJSON(panel.files.trait_details_dir + '/' + traitId + '.json');
  }).then(function (data) {
    renderHero(data);
    renderSupport(data);
    renderMeta(data);
    renderLoci(data);
    renderNotes(data);
  }).catch(function (e) {
    console.error(e);
    document.getElementById('trait-title').textContent = 'Error: ' + e.message;
  });
}

document.addEventListener('DOMContentLoaded', init);
