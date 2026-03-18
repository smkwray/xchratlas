import pandas as pd
import runpy
from pathlib import Path

from xatlas.eqtl import (
    build_locus_query_status,
    eqtl_variant_query_from_components,
    extract_variant_recoder_rsids,
    is_region_query_build_safe,
    normalize_rsid,
    normalize_build_name,
    normalize_rsid_from_locus,
    parse_eqtl_api_payload,
    variant_recoder_input_from_locus,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
PREP_EQTL_INDEX = runpy.run_path(str(REPO_ROOT / "scripts/04_prepare_eqtl_index.py"))
FETCH_EQTL = runpy.run_path(str(REPO_ROOT / "scripts/06_fetch_eqtl_region_hits.py"))


def test_normalize_rsid_accepts_common_forms():
    assert normalize_rsid("rs123") == "rs123"
    assert normalize_rsid(" RS456 ") == "rs456"
    assert normalize_rsid("RS789") == "rs789"


def test_normalize_rsid_rejects_invalid_values():
    assert normalize_rsid(None) is None
    assert normalize_rsid(pd.NA) is None
    assert normalize_rsid("123") is None
    assert normalize_rsid("rs0") is None
    assert normalize_rsid("rs12a3") is None


def test_normalize_rsid_from_locus_prefers_lead_column():
    row = pd.Series({"lead_rsid": "RS100", "rsid": "rs200"})
    assert normalize_rsid_from_locus(row) == "rs100"


def test_region_query_build_safety_requires_matching_known_builds():
    assert is_region_query_build_safe("GRCh37", "hg19")
    assert not is_region_query_build_safe("GRCh37", "GRCh38")
    assert not is_region_query_build_safe("GRCh37", None)
    assert is_region_query_build_safe("GRCh37", None, allow_unknown=True)


def test_normalize_build_name_exposes_canonical_form():
    assert normalize_build_name("hg19") == "GRCh37"
    assert normalize_build_name("b38") == "GRCh38"


def test_parse_eqtl_api_payload_handles_nested_shapes():
    payload = {
        "_embedded": {
            "associations": [
                {
                    "study": {"study_id": "BLUEPRINT"},
                    "molecular_trait": {"molecular_trait_id": "ENSG000001.1"},
                    "gene": {"gene_id": "ENSG000001"},
                    "variant": {"variant_id": "X_123_A_G_b38", "rsid": "RS123"},
                    "pvalue": 1.2e-8,
                }
            ]
        }
    }
    out = parse_eqtl_api_payload(payload)
    assert len(out) == 1
    assert {"study_id", "molecular_trait_id", "gene_id", "variant_id", "rsid", "pvalue"}.issubset(out.columns)
    assert out.loc[0, "study_id"] == "BLUEPRINT"
    assert out.loc[0, "molecular_trait_id"] == "ENSG000001.1"
    assert out.loc[0, "gene_id"] == "ENSG000001"
    assert out.loc[0, "variant_id"] == "X_123_A_G_b38"
    assert out.loc[0, "rsid"] == "rs123"
    assert float(out.loc[0, "pvalue"]) == 1.2e-8


def test_parse_eqtl_api_payload_uses_top_level_and_lead_fallback():
    payload = [{"study_id": "GTEx", "molecular_trait_id": "ENSG2", "variant_id": "X_9_A_C", "p_value": 0.02}]
    out = parse_eqtl_api_payload(payload, lead_rsid="RS999")
    assert len(out) == 1
    assert out.loc[0, "study_id"] == "GTEx"
    assert out.loc[0, "rsid"] == "rs999"
    assert float(out.loc[0, "pvalue"]) == 0.02


def test_build_locus_query_status_normalizes_fields():
    row = build_locus_query_status(
        locus_id="traitA:nonPAR:1-2",
        query_mode="rsid_api",
        query_status="ok",
        query_detail="1 hit",
        lead_rsid="RS77",
        source_build="hg19",
        target_build="GRCh38",
    )
    assert row["locus_id"] == "traitA:nonPAR:1-2"
    assert row["query_mode"] == "rsid_api"
    assert row["query_status"] == "ok"
    assert row["query_detail"] == "1 hit"
    assert row["lead_rsid"] == "rs77"
    assert row["source_build"] == "GRCh37"
    assert row["target_build"] == "GRCh38"


def test_variant_recoder_input_from_locus_uses_vcf_like_grch37_format():
    row = pd.Series({"lead_chr": "chrX", "lead_pos": 147366938, "lead_ref": "G", "lead_alt": "T"})
    assert variant_recoder_input_from_locus(row) == "X 147366938 . G T"


def test_extract_variant_recoder_rsids_prefers_first_normalized_rsid():
    payload = [
        {"T": {"input": "X 147366938 . G T", "id": ["RS889464307", "COSV123"]}},
        {"-": {"input": "X 70070428 . CAA C"}},
    ]
    assert extract_variant_recoder_rsids(payload) == {"X 147366938 . G T": "rs889464307"}


def test_eqtl_variant_query_from_components_formats_chr_variant_string():
    assert eqtl_variant_query_from_components("X", 70850578, "CAA", "C") == "chrX_70850578_CAA_C"


def test_canonicalize_study_id_maps_known_aliases():
    canonicalize_study_id = PREP_EQTL_INDEX["canonicalize_study_id"]
    assert canonicalize_study_id("GTEx_V8") == "GTEx"
    assert canonicalize_study_id("CommonMind_V2") == "CommonMind"


def test_choose_eqtl_study_specs_prefers_harmonized_study_and_group_columns():
    choose_eqtl_study_specs = FETCH_EQTL["choose_eqtl_study_specs"]
    eqtl_index = pd.DataFrame(
        [
            {
                "study_id": "GTEx",
                "study_accession": pd.NA,
                "dataset_id": pd.NA,
                "qtl_group": "Whole_Blood",
                "quant_method": "ge",
                "priority": 1,
            },
            {
                "study_id": "Alasoo_2018",
                "study_accession": "QTS000001",
                "dataset_id": "QTD000006",
                "qtl_group": "macrophage_IFNg",
                "quant_method": "ge",
                "priority": 2,
            },
        ]
    )
    specs = choose_eqtl_study_specs(eqtl_index, 3)
    assert specs == [
        {
            "dataset_id": "QTD000006",
            "study": "Alasoo_2018",
            "qtl_group": "macrophage_IFNg",
            "quant_method": "ge",
        },
    ]


def test_choose_eqtl_study_specs_diversifies_studies_before_extra_datasets():
    choose_eqtl_study_specs = FETCH_EQTL["choose_eqtl_study_specs"]
    eqtl_index = pd.DataFrame(
        [
            {"study_id": "GTEx", "dataset_id": "QTD000116", "qtl_group": "Whole_Blood", "quant_method": "ge", "priority": 1},
            {"study_id": "GTEx", "dataset_id": "QTD000121", "qtl_group": "Muscle", "quant_method": "ge", "priority": 1},
            {"study_id": "CommonMind", "dataset_id": "QTD000075", "qtl_group": "DLPFC", "quant_method": "ge", "priority": 2},
            {"study_id": "Braineac2", "dataset_id": "QTD000041", "qtl_group": "cortex", "quant_method": "ge", "priority": 3},
        ]
    )
    specs = choose_eqtl_study_specs(eqtl_index, 3)
    assert specs == [
        {
            "dataset_id": "QTD000116",
            "study": "GTEx",
            "qtl_group": "Whole_Blood",
            "quant_method": "ge",
        },
        {
            "dataset_id": "QTD000075",
            "study": "CommonMind",
            "qtl_group": "DLPFC",
            "quant_method": "ge",
        },
        {
            "dataset_id": "QTD000041",
            "study": "Braineac2",
            "qtl_group": "cortex",
            "quant_method": "ge",
        },
    ]


def test_resolve_missing_rsids_via_variant_recoder_maps_loci():
    resolve_missing_rsids_via_variant_recoder = FETCH_EQTL["resolve_missing_rsids_via_variant_recoder"]

    class FakeResponse:
        status_code = 200
        text = "ok"

        def json(self):
            return [
                {"T": {"input": "X 147366938 . G T", "id": ["rs889464307"]}},
                {"-": {"input": "X 70070428 . CAA C"}},
            ]

    class FakeSession:
        def post(self, url, json, timeout):
            assert url == "https://example.test/variant_recoder"
            assert json["ids"] == ["X 147366938 . G T", "X 70070428 . CAA C"]
            return FakeResponse()

    loci = pd.DataFrame(
        [
            {"locus_id": "l1", "lead_chr": "X", "lead_pos": 147366938, "lead_ref": "G", "lead_alt": "T", "lead_rsid": pd.NA},
            {"locus_id": "l2", "lead_chr": "X", "lead_pos": 70070428, "lead_ref": "CAA", "lead_alt": "C", "lead_rsid": pd.NA},
            {"locus_id": "l3", "lead_chr": "X", "lead_pos": 1, "lead_ref": "A", "lead_alt": "G", "lead_rsid": "rs1"},
        ]
    )
    resolved, error = resolve_missing_rsids_via_variant_recoder(
        FakeSession(),
        loci,
        source_build="GRCh37",
        variant_recoder_url="https://example.test/variant_recoder",
        timeout_seconds=10.0,
    )
    assert error is None
    assert resolved == {"l1": "rs889464307"}


def test_bridge_loci_to_eqtl_variants_maps_and_validates_exact_variant():
    bridge_loci_to_eqtl_variants = FETCH_EQTL["bridge_loci_to_eqtl_variants"]

    class FakeResponse:
        def __init__(self, *, status_code=200, text="", payload=None):
            self.status_code = status_code
            self.text = text
            self._payload = payload

        def json(self):
            if self._payload is None:
                raise ValueError("no json")
            return self._payload

    class FakeSession:
        def get(self, url, timeout):
            if "/map/human/GRCh37/" in url:
                return FakeResponse(
                    payload={
                        "mappings": [
                            {
                                "mapped": {
                                    "seq_region_name": "X",
                                    "start": 70850578,
                                    "end": 70850580,
                                    "strand": 1,
                                }
                            }
                        ]
                    }
                )
            if "/sequence/region/human/" in url:
                return FakeResponse(text="CAA")
            raise AssertionError(url)

    loci = pd.DataFrame(
        [
            {"locus_id": "l1", "lead_pos": 70070428, "lead_ref": "CAA", "lead_alt": "C"},
        ]
    )
    bridged, errors = bridge_loci_to_eqtl_variants(
        FakeSession(),
        loci,
        source_build="GRCh37",
        target_build="GRCh38",
        grch37_rest_base_url="https://grch37.example.test",
        grch38_sequence_base_url="https://grch38.example.test/sequence/region/human",
        timeout_seconds=10.0,
    )
    assert bridged == {"l1": "chrX_70850578_CAA_C"}
    assert errors == {}
