'use strict';

function listingStatus(trait) {
  var completed = trait.eqtl_lookup_complete_locus_count;
  var partial = completed > 0 && completed < trait.n_loci;
  if (trait.eqtl_supported) return { rank: 5, cls: 'is-supported', label: 'eQTL support' };
  if (trait.eqtl_lookup_hit) return { rank: 4, cls: 'is-lookup-hit', label: 'Returned records' };
  var rank = !trait.n_loci ? 0 : !completed ? 1 : partial ? 2 : 3;
  var labels = ['No chrX loci', 'Not assessed', 'Partly assessed', 'No records'];
  return { rank: rank, cls: 'is-unsupported', label: labels[rank] };
}

function listingGeneLink(symbol, className) {
  var label = symbol.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  return '<a class="gene-link ' + className + '" href="gene.html?gene=' + encodeURIComponent(symbol) + '">' + label + '</a>';
}
