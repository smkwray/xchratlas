'use strict';

var CONFIG = { dataRoot: 'data' };

function fetchJSON(p) {
  return fetch(CONFIG.dataRoot + '/' + p).then(function (r) {
    if (!r.ok) throw new Error('Fetch failed: ' + p); return r.json();
  });
}
function fmt(n) { return Number(n).toLocaleString(); }
function fmtDomain(d) { return d.replace(/_/g, ' ').replace(/\b\w/g, function (c) { return c.toUpperCase(); }); }
function fmtPval(p) {
  // copy: phase 2
  if (!Number.isFinite(p)) return '—';
  return p === 0 ? '0 (reported; precision limit)' : p.toExponential(2);
}
function fmtPos(s, e) { return 'X:' + fmt(s) + '\u2013' + fmt(e); }
function getParam(k) { return new URLSearchParams(location.search).get(k); }
function safeText(v, fallback) {
  if (v === null || v === undefined || v === '') return fallback || '';
  return String(v);
}
function displayTraitName(trait) {
  if (!trait) return '';
  if (trait.query_id === 'operation_code' && trait.coding_description) return trait.coding_description;
  return safeText(trait.description, safeText(trait.query_id, safeText(trait.trait_id, 'Unlabeled trait')));
}
function eqtlStatusMeta(item) {
  if (item.eqtl_supported) {
    return { key: 'supported', supportClass: 'support-yes', badgeClass: 'eqtl-yes', label: 'Candidate support', icon: '\u2713' };
  }
  if (item.eqtl_lookup_hit) {
    return { key: 'lookup', supportClass: 'support-mid', badgeClass: 'eqtl-mid', label: 'Returned records', icon: '\u2022' };
  }
  // copy: phase 2
  var label;
  if (item.n_loci === 0) label = 'No significant locus';
  else if (item.n_loci !== undefined) {
    var complete = item.eqtl_lookup_complete_locus_count;
    label = !complete ? 'Not assessed' : complete < item.n_loci ? 'Partly assessed' : 'Completed with no records';
  } else {
    label = ['complete', 'completed', 'done', 'found', 'ok', 'success', 'succeeded', 'not_found', 'no_hits', 'empty']
      .includes(item.eqtl_lookup_status) ? 'Completed with no records' : 'Not assessed';
  }
  return { key: 'none', supportClass: 'support-no', badgeClass: 'eqtl-no', label: label, icon: '\u2717' };
}
function pickSupportedGene(locus) {
  return (locus.candidate_genes || []).find(function (g) { return g.eqtl_supported && g.eqtl_gene_role !== 'followup'; });
}
function pickLookupGene(locus) {
  return (locus.candidate_genes || []).find(function (g) { return g.eqtl_lookup_hit && g.eqtl_gene_role !== 'followup'; }) ||
    (locus.candidate_genes || []).find(function (g) { return g.eqtl_lookup_hit; });
}

/* ── Hero ──────────────────────────────────── */
function renderHero(data) {
  var t = data.trait, p = data.panel;
  var displayName = displayTraitName(t);
  document.title = displayName + ' \u2014 ' + p.title + ' \u2014 chrXatlas';
  var bpanel = document.getElementById('trait-breadcrumb-panel');
  bpanel.textContent = p.title;
  bpanel.href = 'panel.html?id=' + p.panel_id;
  document.getElementById('trait-breadcrumb-name').textContent = displayName;
  document.getElementById('trait-title').textContent = displayName;
  document.getElementById('trait-domain').textContent = fmtDomain(t.domain);
  var g = document.getElementById('trait-grade');
  g.textContent = t.coverage_grade;
  g.className = 'trait-grade grade-' + t.coverage_grade.toLowerCase();
  document.getElementById('trait-score-value').textContent = t.x_evidence_score.toFixed(0);
  document.getElementById('trait-score-fill').style.width = Math.min(t.x_evidence_score, 100) + '%';
  document.getElementById('ts-loci').textContent = t.n_loci;
  document.getElementById('ts-eqtl-loci').textContent = t.eqtl_supported_locus_count;
  document.getElementById('ts-eqtl-hits').textContent = fmt(t.eqtl_total_hit_count);
  // copy: phase 2
  document.getElementById('ts-pval').textContent = Number.isFinite(t.lead_neglog10_pvalue)
    ? (t.lead_neglog10_is_lower_bound ? '≥' : '') + t.lead_neglog10_pvalue.toFixed(1) : '—';
}

/* ── Metadata ──────────────────────────────── */
function renderMeta(data) {
  var t = data.trait, el = document.getElementById('metadata-grid');
  if (!el) return;
  var sampleSize = t.effective_hq_case_count;
  if (sampleSize === null || sampleSize === undefined || sampleSize === '') {
    sampleSize = t.n_cases_hq_cohort_both_sexes;
  }
  var items = [
    ['Phenocode', t.phenocode || '\u2014'],
    ['Trait type', t.trait_type || '\u2014'],
    ['Modifier', t.modifier || 'none'],
    [(['continuous', 'biomarkers'].includes(t.trait_type) ? 'HQ measured participants' : 'HQ cases'), (sampleSize === null || sampleSize === undefined || sampleSize === '') ? '\u2014' : fmt(sampleSize)],
    ['Populations passing QC', t.num_pops_pass_qc || '\u2014'],
    ['Tier', t.tier || '\u2014'],
    ['Selection rank', t.selection_rank || '\u2014'],
    ['Selection score', t.selection_score || '\u2014'],
    ['Non-PAR loci', t.any_nonpar_locus ? 'Yes' : 'No'],
    ['PAR-only signal', t.par_only_signal ? 'Yes' : 'No'],
    ['Max independent set', t.in_max_independent_set ? 'Yes' : 'No'],
    ['Category', t.category || '\u2014']
  ];
  if (t.query_id === 'operation_code' && t.coding_description) {
    items.unshift(['Source label', t.description || '\u2014']);
    items.splice(items.length - 1, 0, ['Reported operation', t.coding_description]);
  }
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
    // copy: phase 2
    el.innerHTML = '<p class="table-empty">No significant chromosome X loci.</p>';
    return;
  }
  el.innerHTML = data.loci.map(function (locus, i) {
    var genes = locus.candidate_genes || [];
    var status = eqtlStatusMeta(locus);
    return '<details class="locus-card"' + (i < 3 ? ' open' : '') + '>' +
      '<summary class="locus-summary">' +
        '<div class="locus-summary-main">' +
          '<span class="locus-region region-' + locus.x_region.toLowerCase() + '">' + locus.x_region + '</span>' +
          '<span class="locus-position">' + fmtPos(locus.locus_start, locus.locus_end) + '</span>' +
          '<span class="locus-lead-variant">' + (locus.lead_rsid || locus.lead_varid) + '</span>' +
        '</div>' +
        '<div class="locus-summary-stats">' +
          '<span class="locus-pval">−log₁₀(p) ' + (locus.lead_neglog10_is_lower_bound ? '≥ ' : '') + locus.lead_neglog10_pvalue.toFixed(1) + '</span>' +
          '<span class="locus-top-gene">' + locus.top_gene_name + '</span>' +
          '<span class="locus-eqtl-badge ' + status.badgeClass + '">' + status.label + '</span>' +
          '<span class="locus-gene-count">' + genes.length + ' gene' + (genes.length !== 1 ? 's' : '') + '</span>' +
        '</div>' +
      '</summary>' +
      '<div class="locus-body" tabindex="0" role="region" aria-label="Locus details">' +
        // copy: phase 2
        (locus.lead_selection === 'clipped_unresolved' ? '<p>Representative after p-value clipping</p>' : '') +
        (locus.eqtl_query_rsid ? '<p>Lookup rsID: ' + locus.eqtl_query_rsid + '</p>' :
          (String(locus.eqtl_lookup_mode).startsWith('variant_recoder+') ? '<p>Recovered rsID unavailable</p>' : '')) +
        (genes.length > 0 ? geneTable(genes) : '<p style="color:var(--text-tertiary);font-size:0.88rem">No candidate genes mapped.</p>') +
      '</div>' +
    '</details>';
  }).join('');
}

function geneTable(genes) {
  return '<table class="candidate-table"><thead><tr>' +
    '<th>Rank</th><th>Gene</th><th>Biotype</th><th>Role</th><th>Relation</th><th>Distance</th><th>Lookup</th><th>p ≤ 10⁻⁵</th><th>Studies</th><th>Datasets</th><th>Best p</th>' +
    '</tr></thead><tbody>' +
    genes.map(function (g) {
      var role = g.eqtl_gene_role === 'followup' ? 'Follow-up' : 'Candidate';
      var lookup = g.eqtl_lookup_hit ? (g.eqtl_lookup_hit_count + ' hit' + (g.eqtl_lookup_hit_count === 1 ? '' : 's')) : '\u2014';
      var support = g.eqtl_supported
        ? '<span class="eqtl-yes-sm">Yes</span>'
        : (g.eqtl_lookup_hit ? '<span class="eqtl-mid-sm">Lookup-hit</span>' : '\u2014');
      return '<tr class="' + (g.eqtl_supported ? 'gene-supported' : '') + '">' +
        '<td class="num-cell">' + g.candidate_rank + '</td>' +
        '<td><span class="gene-name">' + g.gene_name + '</span></td>' +
        '<td><span class="biotype-badge">' + g.gene_biotype.replace(/_/g, ' ') + '</span></td>' +
        '<td><span class="gene-role-badge gene-role-' + (g.eqtl_gene_role === 'followup' ? 'followup' : 'candidate') + '">' + role + '</span></td>' +
        '<td>' + g.mapping_relation.replace(/_/g, ' ') + '</td>' +
        '<td class="num-cell">' + fmt(g.distance_to_lead_bp) + ' bp</td>' +
        '<td class="num-cell">' + lookup + '</td>' +
        '<td>' + support + '</td>' +
        '<td class="num-cell">' + g.eqtl_study_count + '</td>' +
        '<td class="num-cell">' + g.eqtl_dataset_count + '</td>' +
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
  var status = eqtlStatusMeta(t);
  var supLoci = loci.filter(function (l) { return l.eqtl_supported; });
  var hitLoci = loci.filter(function (l) { return l.eqtl_lookup_hit; });

  var html = '<div class="support-card ' + status.supportClass + '">' +
    '<div class="support-header">' +
      '<span class="support-status-badge ' + status.supportClass + '">' + status.icon + ' ' + status.label + '</span>' +
      '<h3>eQTL assessment</h3>' +
    '</div><div class="support-body">';
  // copy: phase 2
  if (t.n_loci) {
    var complete = t.eqtl_lookup_complete_locus_count;
    var partial = complete > 0 && complete < t.n_loci;
    html += '<p>' + (partial && status.label !== 'Partly assessed' ? 'Partly assessed · ' : '') +
      complete + ' of ' + t.n_loci + ' loci assessed</p>';
  }

  if (supLoci.length > 0) {
    html += '<div class="support-loci">' +
      '<h4>Supporting ' + (supLoci.length === 1 ? 'locus' : 'loci') + '</h4>' +
      supLoci.map(function (l) {
        var supportedGene = pickSupportedGene(l);
        return '<div class="support-locus-row">' +
          '<span class="locus-region region-' + l.x_region.toLowerCase() + '">' + l.x_region + '</span>' +
          '<span class="support-locus-range">chrX ' + fmtMb(l.locus_start) + '\u2013' + fmtMb(l.locus_end) + '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Positional top gene:</strong> ' + l.top_gene_name +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>' + (l.lead_selection === 'clipped_unresolved' ? 'Representative' : 'Lead variant') + ':</strong> ' + (l.lead_rsid || l.lead_varid) +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Lookup:</strong> ' + (l.eqtl_lookup_status || 'unknown') +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>eQTL hits:</strong> ' + l.eqtl_lookup_n_hits +
          '</span>' +
          (supportedGene ? '<span class="support-locus-detail"><strong>Supporting candidate:</strong> <span class="gene-name">' + supportedGene.gene_name + '</span> (' + supportedGene.eqtl_study_count + ' ' + (supportedGene.eqtl_study_count === 1 ? 'study' : 'studies') + ', best p ' + fmtPval(supportedGene.best_eqtl_pvalue) + ')</span>' : '') +
        '</div>';
      }).join('') +
    '</div>';
  } else if (hitLoci.length > 0) {
    html += '<div class="support-loci">' +
      '<h4>Lookup-hit loci</h4>' +
      hitLoci.map(function (l) {
        var lookupGene = pickLookupGene(l);
        return '<div class="support-locus-row">' +
          '<span class="locus-region region-' + l.x_region.toLowerCase() + '">' + l.x_region + '</span>' +
          '<span class="support-locus-range">chrX ' + fmtMb(l.locus_start) + '\u2013' + fmtMb(l.locus_end) + '</span>' +
          '<span class="support-locus-detail"><strong>Positional top gene:</strong> ' + l.top_gene_name + '</span>' +
          '<span class="support-locus-detail"><strong>Lead variant:</strong> ' + (l.lead_rsid || l.lead_varid) + '</span>' +
          '<span class="support-locus-detail"><strong>Lookup:</strong> ' + (l.eqtl_lookup_status || 'unknown') + '</span>' +
          '<span class="support-locus-detail"><strong>eQTL hits:</strong> ' + l.eqtl_lookup_n_hits + '</span>' +
          (lookupGene ? '<span class="support-locus-detail"><strong>Lookup gene:</strong> <span class="gene-name">' + lookupGene.gene_name + '</span> (best p ' + fmtPval(lookupGene.best_eqtl_pvalue) + ')</span>' : '') +
        '</div>';
      }).join('') +
    '</div>';
  }
  html += '</div></div>';

  el.innerHTML = html;
}

/* ── Notes ─────────────────────────────────── */
function renderNotes(data) {
  var t = data.trait, el = document.getElementById('trait-notes-content');
  if (!el) return;
  var notes = [];
  if (t.confidence_notes) notes.push(t.confidence_notes);
  if (t.description_more) notes.push('<strong>Measurement:</strong> ' + t.description_more);
  document.getElementById('trait-notes-section').hidden = notes.length === 0;
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
