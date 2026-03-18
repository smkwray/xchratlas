from __future__ import annotations

PANUKB_DOCS_URL = "https://pan-dev.ukbb.broadinstitute.org/docs/per-phenotype-files/index.html"
PANUKB_DOWNLOADS_URL = "https://pan.ukbb.broadinstitute.org/downloads/"
PANUKB_PHENOTYPE_MANIFEST_URL = (
    "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz"
)
PANUKB_VARIANT_MANIFEST_URL = (
    "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz"
)

ENSEMBL_GRCH37_REST = "https://grch37.rest.ensembl.org"
X_CHROM_LENGTH_GRCH37 = 155_270_560

PAR1_START = 60_001
PAR1_END = 2_699_520
PAR2_START = 154_931_044
PAR2_END = 155_260_560

DEFAULT_X_REGION_LABELS = {
    "PAR1": (PAR1_START, PAR1_END),
    "nonPAR": (PAR1_END + 1, PAR2_START - 1),
    "PAR2": (PAR2_START, PAR2_END),
}

DEFAULT_PANUKB_OUTPUT_COLUMNS = [
    "chr",
    "pos",
    "ref",
    "alt",
    "rsid",
    "varid",
    "af_meta_hq",
    "af_meta",
    "af_cases_meta_hq",
    "af_cases_meta",
    "af_controls_meta_hq",
    "af_controls_meta",
    "beta_meta_hq",
    "beta_meta",
    "se_meta_hq",
    "se_meta",
    "neglog10_pval_meta_hq",
    "neglog10_pval_meta",
    "pval_meta_hq",
    "pval_meta",
    "log_pval_meta_hq",
    "log_pval_meta",
    "neglog10_pval_heterogeneity_hq",
    "neglog10_pval_heterogeneity",
]

EQTL_DATA_ACCESS_URL = "https://www.ebi.ac.uk/eqtl/Data_access/"
EQTL_API_BASE_URL = "https://www.ebi.ac.uk/eqtl/api"
ENSEMBL_GRCH37_VARIANT_RECODER_URL = "https://grch37.rest.ensembl.org/variant_recoder/human"
ENSEMBL_GRCH37_REST_BASE_URL = "https://grch37.rest.ensembl.org"
ENSEMBL_GRCH38_SEQUENCE_BASE_URL = "https://rest.ensembl.org/sequence/region/human"
EQTL_CHRX_GENOTYPES_URL = (
    "https://raw.githubusercontent.com/eQTL-Catalogue/"
    "eQTL-Catalogue-resources/refs/heads/master/data_tables/chrX_genotypes.tsv"
)
EQTL_TABIX_PATHS_URL = (
    "https://raw.githubusercontent.com/eQTL-Catalogue/"
    "eQTL-Catalogue-resources/refs/heads/master/tabix/tabix_ftp_paths.tsv"
)
EQTL_TABIX_IMPORTED_PATHS_URL = (
    "https://raw.githubusercontent.com/eQTL-Catalogue/"
    "eQTL-Catalogue-resources/refs/heads/master/tabix/tabix_ftp_paths_imported.tsv"
)
