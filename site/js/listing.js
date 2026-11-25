'use strict';

// copy: phase 2
function listingStatus(trait) {
  var completed = trait.eqtl_lookup_complete_locus_count;
  var partial = completed > 0 && completed < trait.n_loci;
  var suffix = partial ? ' · Partly assessed' : '';
  if (trait.eqtl_supported) return { rank: 5, cls: 'is-supported', label: 'Candidate support' + suffix };
  if (trait.eqtl_lookup_hit) return { rank: 4, cls: 'is-lookup-hit', label: 'Returned records' + suffix };
  var rank = !trait.n_loci ? 0 : !completed ? 1 : partial ? 2 : 3;
  var labels = ['No significant locus', 'Not assessed', 'Partly assessed', 'Completed with no records'];
  return { rank: rank, cls: 'is-unsupported', label: labels[rank] };
}
