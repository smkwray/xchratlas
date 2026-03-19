'use strict';

/* ── Configuration ───────────────────────────── */
var CONFIG = {
  dataRoot: 'data',
  animDuration: 800,
  staggerMs: 30,
};

var PANEL_ACCENTS = {
  panel_expanded: 'teal',
  panel_mind_risk_core: 'amber',
  panel_blood_biochemistry_core: 'blue',
  panel_independent_set: 'slate',
};

/* ── Helpers ─────────────────────────────────── */
function fetchJSON(path) {
  return fetch(CONFIG.dataRoot + '/' + path)
    .then(function (res) {
      if (!res.ok) throw new Error('Fetch failed: ' + path + ' (' + res.status + ')');
      return res.json();
    });
}

function fmt(n) { return Number(n).toLocaleString(); }

function fmtDomain(d) {
  if (d === null || d === undefined || d === '') return 'Unlabeled';
  return d.replace(/_/g, ' ').replace(/\b\w/g, function (c) { return c.toUpperCase(); });
}

function safeText(value, fallback) {
  if (value === null || value === undefined || value === '') return fallback || '';
  return String(value);
}

function pct(n, d) { return d ? Math.round((n / d) * 100) : 0; }
function eqtlStatusMeta(item) {
  if (item.eqtl_supported) return { cls: 'is-supported', label: 'eQTL-supported' };
  if (item.eqtl_lookup_hit) return { cls: 'is-lookup-hit', label: 'Lookup-hit only' };
  return { cls: 'is-unsupported', label: 'No lookup-hit' };
}

function formatBuildDate(iso) {
  if (!iso) return null;
  var date = new Date(iso);
  if (Number.isNaN(date.getTime())) return iso;
  return date.getUTCFullYear() + '-' +
    String(date.getUTCMonth() + 1).padStart(2, '0') + '-' +
    String(date.getUTCDate()).padStart(2, '0');
}

function el(tag, cls, html) {
  var e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html) e.innerHTML = html;
  return e;
}

/* ── Counter Animation ───────────────────────── */
function animateValue(element, target, duration, suffix) {
  suffix = suffix || '';
  var startTime = performance.now();

  function tick(now) {
    var elapsed = now - startTime;
    var progress = Math.min(elapsed / duration, 1);
    var eased = 1 - Math.pow(1 - progress, 3);
    var current = Math.round(target * eased);
    element.textContent = fmt(current) + suffix;
    if (progress < 1) requestAnimationFrame(tick);
  }

  requestAnimationFrame(tick);
}

/* ── Scroll Observer ─────────────────────────── */
function initScrollAnimations() {
  var staggerIndex = 0;
  var observer = new IntersectionObserver(function (entries) {
    entries.forEach(function (entry) {
      if (!entry.isIntersecting) return;
      var delay = (staggerIndex++) * CONFIG.staggerMs;
      setTimeout(function () {
        entry.target.classList.add('is-visible');
        var counter = entry.target.querySelector('[data-count-to]');
        if (counter && !counter.dataset.counted) {
          var target = parseInt(counter.dataset.countTo, 10);
          var suffix = counter.dataset.countSuffix || '';
          animateValue(counter, target, CONFIG.animDuration, suffix);
          counter.dataset.counted = 'true';
        }
      }, delay);
      observer.unobserve(entry.target);
    });
  }, { threshold: 0.05, rootMargin: '0px 0px 0px 0px' });

  document.querySelectorAll('.animate-on-scroll').forEach(function (el) {
    observer.observe(el);
  });
}

/* ── Render: Hero Metrics ────────────────────── */
function renderHeroMetrics(manifest) {
  var featured = manifest.panels.filter(function (p) { return p.featured; });
  var totalTraits = 0, totalLoci = 0, totalSupported = 0;
  featured.forEach(function (p) {
    totalTraits += p.trait_count;
    totalLoci += p.locus_count;
    totalSupported += p.supported_trait_count;
  });
  var supportRate = pct(totalSupported, totalTraits);

  function setMetric(id, value, suffix) {
    var e = document.getElementById(id);
    if (e) {
      e.dataset.countTo = value;
      e.dataset.countSuffix = suffix || '';
      e.textContent = '0' + (suffix || '');
    }
  }

  setMetric('metric-traits', totalTraits, '');
  setMetric('metric-loci', totalLoci, '');
  setMetric('metric-support', supportRate, '%');

  var buildNote = document.getElementById('metric-build-note');
  if (buildNote) {
    var buildDate = formatBuildDate(manifest.generated_at);
    buildNote.textContent = buildDate ? 'Current data build: ' + buildDate : 'Current data build unavailable';
  }
}

/* ── Render: Panel Cards ─────────────────────── */
function renderPanelCards(manifest) {
  var grid = document.getElementById('panels-grid');
  if (!grid) return;

  manifest.panels.filter(function (p) { return p.featured; }).forEach(function (panel) {
    var accent = PANEL_ACCENTS[panel.panel_id] || 'teal';
    var supportRate = pct(panel.supported_trait_count, panel.trait_count);

    var card = el('div', 'panel-card animate-on-scroll');
    card.setAttribute('data-panel-accent', accent);
    card.setAttribute('data-panel-id', panel.panel_id);

    card.innerHTML =
      '<div class="panel-card-accent"></div>' +
      '<div class="panel-card-body">' +
        '<h3 class="panel-card-title">' + panel.title + '</h3>' +
        '<p class="panel-card-desc">' + panel.description + '</p>' +
        '<div class="panel-card-stats">' +
          '<div class="panel-stat"><span class="panel-stat-value">' + panel.trait_count + '</span><span class="panel-stat-label">traits</span></div>' +
          '<div class="panel-stat"><span class="panel-stat-value">' + fmt(panel.locus_count) + '</span><span class="panel-stat-label">loci</span></div>' +
          '<div class="panel-stat"><span class="panel-stat-value">' + supportRate + '%</span><span class="panel-stat-label">eQTL-supported</span></div>' +
        '</div>' +
        '<div class="panel-card-domains" id="domains-' + panel.panel_id + '"></div>' +
        '<div class="panel-card-top-traits" id="top-traits-' + panel.panel_id + '"></div>' +
        '<a href="panel.html?id=' + panel.panel_id + '" class="panel-card-cta">Explore panel &rarr;</a>' +
      '</div>';

    grid.appendChild(card);
  });
}

/* ── Render: Panel Details (post-summary load) ── */
function updatePanelDetails(summaries) {
  Object.keys(summaries).forEach(function (panelId) {
    var summary = summaries[panelId];

    var domainsEl = document.getElementById('domains-' + panelId);
    if (domainsEl && summary.domain_summary) {
      var domains = summary.domain_summary.slice(0, 6);
      var html = domains.map(function (d) {
        return '<span class="domain-badge">' + fmtDomain(d.domain) + '</span>';
      }).join('');
      if (summary.domain_summary.length > 6) {
        html += '<span class="domain-badge domain-badge-more">+' + (summary.domain_summary.length - 6) + ' more</span>';
      }
      domainsEl.innerHTML = html;
    }

    var traitsEl = document.getElementById('top-traits-' + panelId);
    if (traitsEl && summary.top_traits) {
      var top3 = summary.top_traits.slice(0, 3);
      traitsEl.innerHTML =
        '<div class="panel-top-traits-label">Top traits</div>' +
        top3.map(function (t) {
          return '<div class="panel-trait-row">' +
            '<span class="panel-trait-name">' + t.description + '</span>' +
            '<span class="panel-trait-loci">' + t.n_loci + ' loci</span>' +
          '</div>';
        }).join('');
    }
  });
}

/* ── Render: Domain Composition Bars ─────────── */
function renderDomainBars(summary) {
  var container = document.getElementById('domain-bars');
  if (!container || !summary.domain_summary) return;

  var domains = summary.domain_summary
    .filter(function (d) { return d.locus_count > 0; })
    .sort(function (a, b) { return b.locus_count - a.locus_count; });

  var maxLoci = domains[0] ? domains[0].locus_count : 1;

  domains.forEach(function (d) {
    var width = (d.locus_count / maxLoci) * 100;
    var bar = el('div', 'domain-bar');
    bar.innerHTML =
      '<span class="domain-name">' + fmtDomain(d.domain) + '</span>' +
      '<div class="domain-bar-track">' +
        '<div class="domain-bar-fill" style="width:' + width + '%">' +
          '<span class="domain-bar-value">' + d.locus_count + '</span>' +
        '</div>' +
      '</div>' +
      '<span class="domain-trait-count">' + d.trait_count + ' trait' + (d.trait_count > 1 ? 's' : '') + '</span>';
    container.appendChild(bar);
  });
}

/* ── Render: Top Traits Rail ─────────────────── */
function renderTopTraits(summary) {
  var rail = document.getElementById('traits-rail');
  if (!rail || !summary.top_traits) return;

  summary.top_traits.slice(0, 12).forEach(function (t) {
    var scoreWidth = Math.min(t.x_evidence_score, 100);
    var card = el('div', 'trait-card');
    card.innerHTML =
      '<div class="trait-card-header">' +
        '<span class="trait-domain-badge">' + fmtDomain(t.domain) + '</span>' +
        '<span class="trait-grade grade-' + t.coverage_grade.toLowerCase() + '">' + t.coverage_grade + '</span>' +
      '</div>' +
      '<h4 class="trait-name">' + t.description + '</h4>' +
      '<div class="trait-score">' +
        '<div class="trait-score-bar"><div class="trait-score-fill" style="width:' + scoreWidth + '%"></div></div>' +
        '<span class="trait-score-value">' + t.x_evidence_score.toFixed(0) + '</span>' +
      '</div>' +
      '<div class="trait-meta">' +
        '<div class="trait-meta-item"><span class="trait-meta-label">Loci</span><span class="trait-meta-value">' + t.n_loci + '</span></div>' +
        '<div class="trait-meta-item"><span class="trait-meta-label">eQTL hits</span><span class="trait-meta-value">' + fmt(t.eqtl_total_hit_count) + '</span></div>' +
      '</div>' +
      '<div class="trait-genes">' +
        t.top_candidate_genes.slice(0, 3).map(function (g) {
          return '<span class="gene-tag">' + g + '</span>';
        }).join('') +
      '</div>';
    rail.appendChild(card);
  });
}

/* ── Render: Panel Comparison ────────────────── */
function renderPanelComparison(manifest) {
  var grid = document.getElementById('comparison-grid');
  if (!grid) return;

  var panels = manifest.panels;
  var maxLoci = Math.max.apply(null, panels.map(function (p) { return p.locus_count; }));
  var maxTraits = Math.max.apply(null, panels.map(function (p) { return p.trait_count; }));

  panels.forEach(function (panel) {
    var accent = PANEL_ACCENTS[panel.panel_id] || 'slate';
    var supportRate = pct(panel.supported_trait_count, panel.trait_count);
    var lociW = (panel.locus_count / maxLoci) * 100;
    var traitW = (panel.trait_count / maxTraits) * 100;

    var item = el('div', 'comparison-item');
    item.setAttribute('data-panel-accent', accent);
    item.innerHTML =
      '<div class="comparison-header">' +
        '<span class="comparison-title">' + panel.title + '</span>' +
        (panel.featured ? '<span class="comparison-badge">Featured</span>' : '') +
      '</div>' +
      '<div class="comparison-bars">' +
        compRow('Traits', traitW, panel.trait_count) +
        compRow('Loci', lociW, fmt(panel.locus_count)) +
        compRow('eQTL support', supportRate, supportRate + '%') +
      '</div>';
    grid.appendChild(item);
  });
}

function compRow(label, width, value) {
  return '<div class="comparison-row">' +
    '<span class="comparison-label">' + label + '</span>' +
    '<div class="comparison-bar-track"><div class="comparison-bar-fill" style="width:' + width + '%"></div></div>' +
    '<span class="comparison-value">' + value + '</span>' +
  '</div>';
}

/* ── Render: Discovery Pool ──────────────────── */
function discoveryTraitRow(pool, trait) {
  var href = 'trait.html?panel=' + pool.panel_id + '&trait=' + encodeURIComponent(trait.trait_id);
  var status = eqtlStatusMeta(trait);
  var topGenes = (trait.top_candidate_genes || []).slice(0, 3);

  return '<a class="discovery-trait-row" href="' + href + '">' +
    '<div class="discovery-trait-main">' +
      '<span class="discovery-trait-name">' + trait.description + '</span>' +
      '<div class="discovery-trait-meta">' +
        '<span class="discovery-trait-domain">' + fmtDomain(trait.domain) + '</span>' +
        '<span class="discovery-trait-status ' + status.cls + '">' + status.label + '</span>' +
        '<span class="discovery-trait-fact">score ' + trait.x_evidence_score.toFixed(0) + '</span>' +
        '<span class="discovery-trait-fact">' + trait.n_loci + ' loci</span>' +
        '<span class="discovery-trait-fact">' + trait.eqtl_total_hit_count + ' eQTL hits</span>' +
      '</div>' +
    '</div>' +
    '<div class="discovery-trait-side">' +
      (topGenes.length
        ? '<span class="discovery-trait-genes">' + topGenes.join(', ') + '</span>'
        : '<span class="discovery-trait-genes discovery-trait-genes-empty">No highlighted genes</span>') +
    '</div>' +
  '</a>';
}

function discoveryDrawer(title, count, rows, subtitle, open) {
  return '<details class="discovery-drawer"' + (open ? ' open' : '') + '>' +
    '<summary class="discovery-drawer-summary">' +
      '<div>' +
        '<span class="discovery-drawer-title">' + title + '</span>' +
        '<span class="discovery-drawer-subtitle">' + subtitle + '</span>' +
      '</div>' +
      '<span class="discovery-drawer-count">' + fmt(count) + '</span>' +
    '</summary>' +
    '<div class="discovery-drawer-body">' + rows.join('') + '</div>' +
  '</details>';
}

function renderDiscoveryPool(manifest, summary, traits) {
  var card = document.getElementById('discovery-card');
  if (!card) return;
  var pool = manifest.panels.find(function (p) { return !p.featured; });
  if (!pool) return;

  var topDomains = summary.domain_summary
    .filter(function (d) { return d.locus_count > 0; })
    .sort(function (a, b) { return b.locus_count - a.locus_count; })
    .slice(0, 8);
  var allTraits = (traits || []).slice().sort(function (a, b) {
    return b.x_evidence_score - a.x_evidence_score ||
      b.n_loci - a.n_loci ||
      safeText(a.description).localeCompare(safeText(b.description));
  });
  var unsupportedTraits = allTraits.filter(function (t) { return !t.eqtl_supported; });

  card.innerHTML =
    '<div class="discovery-stats">' +
      dStat(pool.trait_count, 'traits scanned') +
      dStat(pool.lookup_hit_trait_count, 'eQTL lookup-hit') +
      dStat(pool.supported_trait_count, 'eQTL-supported') +
      dStat(fmt(pool.locus_count), 'loci') +
      dStat(pool.zero_locus_trait_count, 'zero-locus') +
    '</div>' +
    '<div class="discovery-domains">' +
      '<h4>Top domains by locus count</h4>' +
      '<div class="discovery-domain-list">' +
        topDomains.map(function (d) {
          return '<div class="discovery-domain-row">' +
            '<span class="discovery-domain-name">' + fmtDomain(d.domain) + '</span>' +
            '<span class="discovery-domain-value">' + d.locus_count + ' loci &middot; ' + d.trait_count + ' traits</span>' +
          '</div>';
        }).join('') +
      '</div>' +
    '</div>' +
    '<div class="discovery-traits">' +
      discoveryDrawer(
        'Traits without strict eQTL support',
        unsupportedTraits.length,
        unsupportedTraits.map(function (trait) { return discoveryTraitRow(pool, trait); }),
        'Includes zero-locus traits plus lookup-hit-only traits that did not meet the strict support rule.',
        true
      ) +
      discoveryDrawer(
        'All tested discovery traits',
        allTraits.length,
        allTraits.map(function (trait) { return discoveryTraitRow(pool, trait); }),
        'Full Pan-UKB max independent-set discovery pool used for curation.',
        false
      ) +
    '</div>' +
    '<p class="discovery-note">The discovery pool is a broad curation reservoir, not a featured panel. ' +
    'Traits with zero genome-wide significant chrX loci are retained for completeness but represent the ' +
    'absence of detectable signal in the current build.</p>';
}

function dStat(val, label) {
  return '<div class="discovery-stat">' +
    '<span class="discovery-stat-value">' + val + '</span>' +
    '<span class="discovery-stat-label">' + label + '</span>' +
  '</div>';
}

/* ── Unique Domains ──────────────────────────── */
function countUniqueDomains(summaries, panels) {
  var domains = {};
  panels.forEach(function (p) {
    if (p.featured && summaries[p.panel_id]) {
      summaries[p.panel_id].domain_summary.forEach(function (d) {
        domains[d.domain] = true;
      });
    }
  });
  return Object.keys(domains).length;
}

/* ── Smooth Scroll ───────────────────────────── */
function initSmoothScroll() {
  document.querySelectorAll('a[href^="#"]').forEach(function (link) {
    link.addEventListener('click', function (e) {
      var href = link.getAttribute('href');
      if (href === '#') return;
      var target = document.querySelector(href);
      if (target) {
        e.preventDefault();
        var h = document.querySelector('.header');
        var offset = h ? h.offsetHeight + 20 : 20;
        var top = target.getBoundingClientRect().top + window.scrollY - offset;
        window.scrollTo({ top: top, behavior: 'smooth' });
        var nav = document.getElementById('nav');
        if (nav) nav.classList.remove('open');
      }
    });
  });
}

/* ── Main ────────────────────────────────────── */
function init() {
  fetchJSON('manifest.json').then(function (manifest) {
    renderHeroMetrics(manifest);
    renderPanelCards(manifest);
    renderPanelComparison(manifest);

    var summaries = {};
    var panelTraits = {};
    var promises = manifest.panels.map(function (panel) {
      return fetchJSON(panel.files.summary).then(function (summary) {
        summaries[panel.panel_id] = summary;
      });
    });
    var discoveryPanel = manifest.panels.find(function (panel) { return !panel.featured; });
    if (discoveryPanel) {
      promises.push(
        fetchJSON(discoveryPanel.files.traits).then(function (traits) {
          panelTraits[discoveryPanel.panel_id] = traits;
        })
      );
    }

    return Promise.all(promises).then(function () {
      var domainCount = countUniqueDomains(summaries, manifest.panels);
      var domainEl = document.getElementById('metric-domains');
      if (domainEl) {
        domainEl.dataset.countTo = domainCount;
        domainEl.dataset.countSuffix = '';
        domainEl.textContent = '0';
      }

      updatePanelDetails(summaries);

      if (summaries.panel_expanded) {
        renderDomainBars(summaries.panel_expanded);
        renderTopTraits(summaries.panel_expanded);
      }

      if (summaries.panel_independent_set) {
        renderDiscoveryPool(manifest, summaries.panel_independent_set, panelTraits.panel_independent_set || []);
      }

      initScrollAnimations();
      initSmoothScroll();
      loadSearchIndex(manifest);
    });
  }).catch(function (err) {
    console.error('chrXatlas init failed:', err);
    var sub = document.querySelector('.hero-subtitle');
    if (sub) {
      sub.innerHTML += '<br><small style="color:#ef4444">Data loading error: ' + err.message +
        '. Ensure the data files are in site/data/.</small>';
    }
  });
}

/* ── Search ───────────────────────────────────── */
var searchIndex = [];

function loadSearchIndex(manifest) {
  var promises = manifest.panels.map(function (panel) {
    return fetchJSON(panel.files.traits).then(function (traits) {
      traits.forEach(function (t) {
        searchIndex.push({
          trait_id: t.trait_id, description: t.description,
          domain: t.domain, query_id: t.query_id,
          x_evidence_score: t.x_evidence_score,
          top_candidate_genes: t.top_candidate_genes || [],
          panel_id: panel.panel_id, panel_title: panel.title
        });
      });
    });
  });
  return Promise.all(promises).then(initSearch);
}

function initSearch() {
  var input = document.getElementById('search-input');
  var results = document.getElementById('search-results');
  if (!input || !results) return;

  input.addEventListener('input', function () {
    var q = input.value.toLowerCase().trim();
    if (q.length < 2) { results.innerHTML = ''; results.style.display = 'none'; return; }
    var matches = searchIndex.filter(function (t) {
      var d = (t.description || '').toLowerCase();
      var dm = (t.domain || '').toLowerCase();
      var qi = (t.query_id || '').toLowerCase();
      return d.includes(q) || dm.includes(q) || qi.includes(q) ||
        (t.top_candidate_genes || []).some(function (g) { return g && g.toLowerCase().includes(q); });
    }).slice(0, 12);

    if (matches.length === 0) {
      results.innerHTML = '<div class="search-empty">No traits found</div>';
      results.style.display = 'block'; return;
    }
    results.innerHTML = matches.map(function (t) {
      var score = t.x_evidence_score != null ? t.x_evidence_score.toFixed(0) : '—';
      return '<a href="trait.html?panel=' + t.panel_id + '&trait=' + encodeURIComponent(t.trait_id) + '" class="search-result-item">' +
        '<div class="search-result-name">' + (t.description || t.query_id) + '</div>' +
        '<div class="search-result-meta">' +
          '<span class="search-result-panel">' + (t.panel_title || '') + '</span>' +
          '<span class="search-result-domain">' + fmtDomain(t.domain || '') + '</span>' +
          '<span class="search-result-score">' + score + '</span>' +
        '</div></a>';
    }).join('');
    results.style.display = 'block';
  });

  document.addEventListener('click', function (e) {
    if (!e.target.closest('.search-container')) results.style.display = 'none';
  });
}

document.addEventListener('DOMContentLoaded', init);
