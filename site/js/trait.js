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
  if (!Number.isFinite(p)) return '\u2014';
  return p === 0
    ? '<abbr class="pval-bound" title="Reported as 0 at the precision limit">0 (reported)</abbr>'
    : p.toExponential(2);
}
function fmtLookupStatus(s) {
  if (s === 'complete') return 'complete';
  if (s === 'no_hits') return 'no records';
  if (s === 'skipped_no_rsid') return 'not assessed';
  return s ? s.replace(/_/g, ' ') : 'unknown';
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
    return { key: 'supported', supportClass: 'support-yes', badgeClass: 'eqtl-yes', label: 'eQTL support', icon: '\u2713' };
  }
  if (item.eqtl_lookup_hit) {
    return { key: 'lookup', supportClass: 'support-mid', badgeClass: 'eqtl-mid', label: 'Returned records', icon: '\u2022' };
  }
  var label;
  if (item.n_loci === 0) label = 'No chrX loci';
  else if (item.n_loci !== undefined) {
    var complete = item.eqtl_lookup_complete_locus_count;
    label = !complete ? 'Not assessed' : complete < item.n_loci ? 'Partly assessed' : 'No records';
  } else {
    label = ['complete', 'completed', 'done', 'found', 'ok', 'success', 'succeeded', 'not_found', 'no_hits', 'empty']
      .includes(item.eqtl_lookup_status) ? 'No records' : 'Not assessed';
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
  var pvalEl = document.getElementById('ts-pval');
  if (Number.isFinite(t.lead_neglog10_pvalue)) {
    if (t.lead_neglog10_is_lower_bound) {
      pvalEl.innerHTML = '<abbr class="pval-bound" title="p-value reached the 1e-300 conversion floor">&ge; ' +
        t.lead_neglog10_pvalue.toFixed(1) + '</abbr>';
    } else {
      pvalEl.textContent = t.lead_neglog10_pvalue.toFixed(1);
    }
  } else {
    pvalEl.textContent = '\u2014';
  }
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
    ['Modifier', t.modifier || '\u2014'],
    ['Source file', t.pheno_sex || '\u2014'],
    [(['continuous', 'biomarkers'].includes(t.trait_type) ? 'HQ measured participants' : 'HQ cases'), (sampleSize === null || sampleSize === undefined || sampleSize === '') ? '\u2014' : fmt(sampleSize)],
    ['Populations passing QC', t.num_pops_pass_qc || '\u2014'],
    ['Non-PAR loci', t.any_nonpar_locus ? 'Yes' : 'No'],
    ['PAR-only signal', t.par_only_signal ? 'Yes' : 'No'],
    ['Max independent set', t.in_max_independent_set ? 'Yes' : 'No'],
    ['UK Biobank category', t.category || '\u2014']
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
  document.getElementById('loci-heading').textContent = 'Chromosome X loci (GRCh37; ' + data.loci.length + ')';
  if (data.loci.length === 0) {
    el.innerHTML = '<p class="table-empty">No tested chromosome X variant passed p &lt; 5 &times; 10<sup>&minus;8</sup>.</p>';
    return;
  }
  el.innerHTML = data.loci.map(function (locus, i) {
    var genes = locus.candidate_genes || [];
    var status = eqtlStatusMeta(locus);
    var pvalHtml = locus.lead_neglog10_is_lower_bound
      ? '<abbr class="pval-bound" title="p-value reached the 1e-300 conversion floor">&ge; ' + locus.lead_neglog10_pvalue.toFixed(1) + '</abbr>'
      : locus.lead_neglog10_pvalue.toFixed(1);
    return '<details class="locus-card"' + (i < 3 ? ' open' : '') + '>' +
      '<summary class="locus-summary">' +
        '<div class="locus-summary-main">' +
          '<span class="locus-region region-' + locus.x_region.toLowerCase() + '">' + locus.x_region + '</span>' +
          '<span class="locus-position">' + fmtPos(locus.locus_start, locus.locus_end) + '</span>' +
          '<span class="locus-lead-variant">' + (locus.lead_rsid || locus.lead_varid) + '</span>' +
        '</div>' +
        '<div class="locus-summary-stats">' +
          '<span class="locus-pval">&minus;log<sub>10</sub>(p) ' + pvalHtml + '</span>' +
          '<span class="locus-top-gene">' + locus.top_gene_name + '</span>' +
          '<span class="locus-eqtl-badge ' + status.badgeClass + '">' + status.label + '</span>' +
          '<span class="locus-gene-count">' + genes.length + ' gene' + (genes.length !== 1 ? 's' : '') + '</span>' +
        '</div>' +
      '</summary>' +
      '<div class="locus-body" tabindex="0" role="region" aria-label="Locus details">' +
        (locus.lead_selection === 'clipped_unresolved'
          ? '<p class="locus-meta-note">p-value reached the 10<sup>&minus;300</sup> floor; shown variant is the recorded representative, and gene/eQTL results refer to it.</p>'
          : '') +
        (locus.eqtl_query_rsid
          ? '<p class="locus-meta-note">eQTL lookup rsID: <code>' + locus.eqtl_query_rsid + '</code></p>'
          : (String(locus.eqtl_lookup_mode).startsWith('variant_recoder+')
            ? '<p class="locus-meta-note">eQTL lookup used a recovered rsID not retained in this export.</p>'
            : '')) +
        (genes.length > 0 ? geneTable(genes) : '<p class="locus-meta-note">No candidate genes mapped.</p>') +
      '</div>' +
    '</details>';
  }).join('');
}

function geneTable(genes) {
  return '<table class="candidate-table stack-table"><thead><tr>' +
    '<th>Rank</th><th>Gene</th><th>Biotype</th><th>Role</th><th>Mapping</th><th>Distance</th><th>Records</th><th>p &le; 10<sup>&minus;5</sup></th><th>Studies (&le;10<sup>&minus;5</sup>)</th><th>Datasets</th><th>Best eQTL p</th>' +
    '</tr></thead><tbody>' +
    genes.map(function (g) {
      var role = g.eqtl_gene_role === 'followup' ? 'Follow-up' : 'Candidate';
      var lookup = g.eqtl_lookup_hit ? fmt(g.eqtl_lookup_hit_count) : '\u2014';
      var support = g.eqtl_supported
        ? '<span class="eqtl-yes-sm">Yes</span>'
        : (g.eqtl_lookup_hit ? '<span class="eqtl-mid-sm">No</span>' : '\u2014');
      var geneHref = 'gene.html?gene=' + encodeURIComponent(g.gene_id || g.gene_name);
      return '<tr class="' + (g.eqtl_supported ? 'gene-supported' : '') + '">' +
        '<td class="num-cell" data-label="Rank">' + g.candidate_rank + '</td>' +
        '<td class="candidate-name"><a class="gene-name gene-name-link" href="' + geneHref + '">' + g.gene_name + '</a></td>' +
        '<td data-label="Biotype"><span class="biotype-badge">' + g.gene_biotype.replace(/_/g, ' ') + '</span></td>' +
        '<td data-label="Role"><span class="gene-role-badge gene-role-' + (g.eqtl_gene_role === 'followup' ? 'followup' : 'candidate') + '">' + role + '</span></td>' +
        '<td data-label="Mapping">' + g.mapping_relation.replace(/_/g, ' ') + '</td>' +
        '<td class="num-cell" data-label="Distance">' + fmt(g.distance_to_lead_bp) + ' bp</td>' +
        '<td class="num-cell" data-label="Records">' + lookup + '</td>' +
        '<td data-label="p ≤ 10⁻⁵">' + support + '</td>' +
        '<td class="num-cell" data-label="Studies (≤10⁻⁵)">' + g.eqtl_study_count + '</td>' +
        '<td class="num-cell" data-label="Datasets">' + g.eqtl_dataset_count + '</td>' +
        '<td class="num-cell" data-label="Best eQTL p">' + fmtPval(g.best_eqtl_pvalue) + '</td></tr>';
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
      '<h3>eQTL lookups</h3>' +
    '</div><div class="support-body">';

  if (t.n_loci) {
    var complete = t.eqtl_lookup_complete_locus_count;
    var summaryText = complete + ' of ' + t.n_loci + ' ' + (t.n_loci === 1 ? 'locus' : 'loci') + ' assessed';
    if (complete < t.n_loci) summaryText += ' (unassessed loci lacked a usable rsID)';
    if (supLoci.length === 0 && hitLoci.length === 0 && complete > 0) summaryText += '; no records returned';
    html += '<p class="support-assess">' + summaryText + '.</p>';
  } else {
    html += '<p class="support-assess">No significant chromosome X loci to query.</p>';
  }

  if (supLoci.length > 0) {
    html += '<div class="support-loci">' +
      '<h4>eQTL-supported ' + (supLoci.length === 1 ? 'locus' : 'loci') + '</h4>' +
      supLoci.map(function (l) {
        var supportedGene = pickSupportedGene(l);
        return '<div class="support-locus-row">' +
          '<span class="locus-region region-' + l.x_region.toLowerCase() + '">' + l.x_region + '</span>' +
          '<span class="support-locus-range">chrX ' + fmtMb(l.locus_start) + '\u2013' + fmtMb(l.locus_end) + '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Rank 1 gene:</strong> ' + l.top_gene_name +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>' + (l.lead_selection === 'clipped_unresolved' ? 'Representative' : 'Lead variant') + ':</strong> ' + (l.lead_rsid || l.lead_varid) +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Lookup:</strong> ' + fmtLookupStatus(l.eqtl_lookup_status) +
          '</span>' +
          '<span class="support-locus-detail">' +
            '<strong>Records:</strong> ' + fmt(l.eqtl_lookup_n_hits) +
          '</span>' +
          (supportedGene ? '<span class="support-locus-detail"><strong>Candidate (p &le; 10<sup>&minus;5</sup>):</strong> <span class="gene-name">' + supportedGene.gene_name + '</span> (' + supportedGene.eqtl_study_count + ' ' + (supportedGene.eqtl_study_count === 1 ? 'study' : 'studies') + ', best p ' + fmtPval(supportedGene.best_eqtl_pvalue) + ')</span>' : '') +
        '</div>';
      }).join('') +
    '</div>';
  } else if (hitLoci.length > 0) {
    html += '<div class="support-loci">' +
      '<h4>Loci with returned records</h4>' +
      hitLoci.map(function (l) {
        var lookupGene = pickLookupGene(l);
        return '<div class="support-locus-row">' +
          '<span class="locus-region region-' + l.x_region.toLowerCase() + '">' + l.x_region + '</span>' +
          '<span class="support-locus-range">chrX ' + fmtMb(l.locus_start) + '\u2013' + fmtMb(l.locus_end) + '</span>' +
          '<span class="support-locus-detail"><strong>Rank 1 gene:</strong> ' + l.top_gene_name + '</span>' +
          '<span class="support-locus-detail"><strong>' + (l.lead_selection === 'clipped_unresolved' ? 'Representative' : 'Lead variant') + ':</strong> ' + (l.lead_rsid || l.lead_varid) + '</span>' +
          '<span class="support-locus-detail"><strong>Lookup:</strong> ' + fmtLookupStatus(l.eqtl_lookup_status) + '</span>' +
          '<span class="support-locus-detail"><strong>Records:</strong> ' + fmt(l.eqtl_lookup_n_hits) + '</span>' +
          (lookupGene ? '<span class="support-locus-detail"><strong>Best lookup gene:</strong> <span class="gene-name">' + lookupGene.gene_name + '</span> (best p ' + fmtPval(lookupGene.best_eqtl_pvalue) + ')</span>' : '') +
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
  var section = document.getElementById('trait-notes-section');
  if (!el || !section) return;
  if (!t.description_more) {
    section.hidden = true;
    return;
  }
  section.hidden = false;
  el.innerHTML = '<p>' + t.description_more + '</p>';
}

/* ── Init ──────────────────────────────────── */
function init() {
  var panelId = getParam('panel'), traitId = getParam('trait');
  if (!panelId || !traitId) {
    document.getElementById('trait-title').textContent = 'No trait selected';
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
    document.getElementById('trait-title').textContent = 'Could not load trait';
  });
}

document.addEventListener('DOMContentLoaded', init);
