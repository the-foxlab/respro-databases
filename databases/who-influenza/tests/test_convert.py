"""Tests for WHO Influenza conversion logic (databases/who-influenza)."""

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import convert  # noqa: E402
import resolve_references  # noqa: E402
import who_parsers  # noqa: E402


class TestParseSubstitutionCell:
    """Substitution cell parsing: components, mixed alleles, deletions."""

    def test_single_substitution(self):
        comps, foot = convert.parse_substitution_cell("H275Y")
        assert comps == [{"kind": "sub", "ref": "H", "pos": 275, "alts": ["Y"]}]
        assert foot == ""

    def test_multi_allele_mixed(self):
        comps, foot = convert.parse_substitution_cell("D151E/N+H275Y")
        assert comps[0] == {"kind": "mixed", "ref": "D", "pos": 151,
                            "alts": ["E", "N"]}
        assert comps[1] == {"kind": "sub", "ref": "H", "pos": 275,
                            "alts": ["Y"]}
        assert foot == ""

    def test_combo_of_substitutions(self):
        comps, foot = convert.parse_substitution_cell("E119D+I222L")
        assert len(comps) == 2
        assert comps[0]["pos"] == 119 and comps[1]["pos"] == 222
        assert all(c["kind"] == "sub" for c in comps)

    def test_deletion_en_dash(self):
        comps, foot = convert.parse_substitution_cell("Del 245\u2013248")
        assert comps == [{"kind": "del", "start": 245, "end": 248}]
        assert foot == ""

    def test_deletion_hyphen_with_footnote(self):
        comps, foot = convert.parse_substitution_cell("Del 247-250g")
        assert comps == [{"kind": "del", "start": 247, "end": 250}]
        assert foot == "g"

    def test_footnote_suffix_on_substitution(self):
        comps, foot = convert.parse_substitution_cell("V116Aa")
        assert comps[0]["alts"] == ["A"]
        assert foot == "a"

    @pytest.mark.parametrize("cell", ["", "275", "H275", "Ins 245", "H275Y?"])
    def test_invalid_cells_raise(self, cell):
        with pytest.raises(ValueError):
            convert.parse_substitution_cell(cell)


class TestParseNaiCell:
    """NAI drug-cell parsing: phenotype tokens, fold ranges, footnotes."""

    def test_simple_phenotype_with_fold(self):
        phenotype, fold, note = convert.parse_nai_cell("NI (3)")
        assert phenotype == "normal inhibition"
        assert fold == "3"
        assert note == ""

    def test_greater_than_fold(self):
        phenotype, fold, note = convert.parse_nai_cell("HRI (>7692)")
        assert phenotype == "highly reduced inhibition"
        assert fold == "7692"
        assert ">7692" in note

    def test_less_than_fold(self):
        phenotype, fold, note = convert.parse_nai_cell("NI (<2)")
        assert phenotype == "normal inhibition"
        assert fold == "2"
        assert "<2" in note

    def test_range_takes_max_endpoint(self):
        phenotype, fold, note = convert.parse_nai_cell("RI/HRI (12\u2013138)")
        assert phenotype == "highly reduced inhibition"
        assert fold == "138"
        assert "12\u2013138" in note

    def test_slash_space_combo_picks_most_severe(self):
        phenotype, fold, note = convert.parse_nai_cell("NI /HRI (3\u20137310)")
        assert phenotype == "highly reduced inhibition"
        assert fold == "7310"

    def test_footnote_after_parenthesis(self):
        phenotype, fold, note = convert.parse_nai_cell("NI (<2)j")
        assert phenotype == "normal inhibition"
        assert fold == "2"

    def test_unknown_cell(self):
        assert convert.parse_nai_cell("?") == ("", "", "")

    def test_unknown_fold(self):
        phenotype, fold, note = convert.parse_nai_cell("NI (?)")
        assert phenotype == "normal inhibition"
        assert fold == ""
        assert note == ""

    def test_unparseable_raises(self):
        with pytest.raises(ValueError):
            convert.parse_nai_cell("275 (3)")

    def test_unknown_phenotype_token_fails_closed(self):
        with pytest.raises(KeyError):
            convert.parse_nai_cell("banana (3)")


class TestParsePaClassification:
    """PA classification vocabulary mapping."""

    @pytest.mark.parametrize("raw,expected", [
        ("Normal", "normal"),
        ("Suspected reduced", "suspected reduced"),
        ("Reduced", "reduced susceptibility"),
        ("?", ""),
        ("", ""),
    ])
    def test_mapping(self, raw, expected):
        assert convert.parse_pa_classification(raw) == expected

    def test_unknown_raises(self):
        with pytest.raises(ValueError):
            convert.parse_pa_classification("Extreme")


class TestParseFoldValue:
    """Fold-change cell parsing for the PA table."""

    def test_plain_value(self):
        fold, note = convert.parse_fold_value("3.2", metric="EC50")
        assert fold == "3.2"
        assert note == ""

    def test_greater_than(self):
        fold, note = convert.parse_fold_value(">7692", metric="EC50")
        assert fold == "7692"
        assert ">7692" in note

    def test_range_takes_max(self):
        fold, note = convert.parse_fold_value("12\u2013138", metric="EC50")
        assert fold == "138"
        assert "12\u2013138" in note

    def test_space_thousands_separator(self):
        fold, note = convert.parse_fold_value(">10 000", metric="EC50")
        assert fold == "10000"
        assert ">10 000" in note

    def test_comma_thousands_separator(self):
        fold, note = convert.parse_fold_value(">10,000", metric="EC50")
        assert fold == "10000"

    def test_unknown(self):
        assert convert.parse_fold_value("?", metric="EC50") == ("", "")


class TestCanonicalizeDeletion:
    """Deletion canonical form: anchored reference block."""

    SEQ = "MKTIIALSYIFCLVFAQ"*10 + "GSASGK" + "AAAA"*10

    def test_anchored_form(self):
        # Delete positions 172-175 (block "SASG"), anchor 171 = 'G'.
        canon = convert.canonicalize_deletion(172, 175, self.SEQ)
        anchor, reference, anchor_aa = canon
        assert anchor == 171
        assert reference == "GSASG"
        assert anchor_aa == "G"

    def test_out_of_range_returns_none(self):
        assert convert.canonicalize_deletion(len(self.SEQ) + 1,
                                             len(self.SEQ) + 4,
                                             self.SEQ) is None


class TestResolvePosition:
    """Position resolution: direct validation and H5N1 f-suffix offset."""

    # 30-residue synthetic reference; H at 21, N at 11.
    SEQ = "AAAAAAAAAAAAAAAAAAAANAAAAAAAAA"

    def test_direct_match(self):
        pos, note = convert.resolve_position(21, "N", self.SEQ, {}, "A(H3N2)")
        assert pos == 21
        assert note == ""

    def test_direct_mismatch_fails(self):
        pos, note = convert.resolve_position(21, "E", self.SEQ, {}, "A(H3N2)")
        assert pos is None
        assert "mismatch" in note

    def test_h5n1_footnote_offset(self):
        # 'H275Yf' with 20-residue stalk offset: full-length position 275
        # maps to reference position 255.
        seq = "A"*254 + "H" + "A"*20
        footnotes = {"f": "full-length N1 numbering"}
        pos, note = convert.resolve_position(275, "H", seq, footnotes, "H5N1")
        assert pos == 255
        assert "20" in note

    def test_offset_only_applies_to_h5n1(self):
        seq = "A"*254 + "H" + "A"*20
        footnotes = {"f": "full-length N1 numbering"}
        pos, _ = convert.resolve_position(275, "H", seq, footnotes, "A(H3N2)")
        assert pos is None


class TestMergePolicies:
    """Dedup/merge policies: phenotype severity, fold max, unique join."""

    def test_phenotype_most_severe_wins(self):
        phenotype, note = convert.merge_phenotype(
            ["normal inhibition", "highly reduced inhibition", "reduced inhibition"]
        )
        assert phenotype == "highly reduced inhibition"
        assert "Phenotypes in source: " in note
        assert "most severe used: highly reduced inhibition" in note

    def test_phenotype_single_no_note(self):
        phenotype, note = convert.merge_phenotype(["reduced inhibition"])
        assert phenotype == "reduced inhibition"
        assert note == ""

    def test_fold_max_endpoint_wins(self):
        fold, note = convert.merge_fold([14.0, 138.0])
        assert fold == "138"
        assert "14, 138" in note

    def test_fold_single_no_note(self):
        fold, note = convert.merge_fold([3.0])
        assert fold == "3"
        assert note == ""

    def test_join_unique_deduplicates_comma_joined_elements(self):
        assert convert.join_unique(
            ["doi:10.1/a,doi:10.1/b", "doi:10.1/b", "doi:10.1/c"]
        ) == "doi:10.1/a,doi:10.1/b,doi:10.1/c"


class TestExpandOrigin:
    """Origin-column abbreviation expansion per the WHO annotation footnotes."""

    def test_single_abbreviations(self):
        assert convert.expand_origin("Sur", "human-nai") == "surveillance studies"
        assert convert.expand_origin("RG", "avian-nai") == "reverse genetics"
        assert convert.expand_origin("RecNA", "human-nai") == "recombinant NA"
        assert convert.expand_origin("P-p", "human-nai") == "plaque purification"

    def test_semicolon_and_comma_segments(self):
        assert convert.expand_origin("Sur; RG", "human-nai") == (
            "surveillance studies; reverse genetics"
        )
        assert convert.expand_origin("RG, Sur, RecNA", "avian-nai") == (
            "reverse genetics, surveillance studies, recombinant NA"
        )
        assert convert.expand_origin("Sur;RG", "human-nai") == (
            "surveillance studies; reverse genetics"
        )

    def test_slash_and_plus_subsegments(self):
        assert convert.expand_origin("Clin/Ose", "human-nai") == (
            "clinical detection/oseltamivir used"
        )
        assert convert.expand_origin("Clin/Ose+Zan", "human-nai") == (
            "clinical detection/oseltamivir used+zanamivir used"
        )
        assert convert.expand_origin("RG; in vitro/Zan", "avian-nai") == (
            "reverse genetics; in vitro/zanamivir used"
        )

    def test_clin_and_no_differ_between_nai_and_pa(self):
        assert convert.expand_origin("Clin/No", "human-nai") == (
            "clinical detection/no NAI used"
        )
        assert convert.expand_origin("Sur/No", "pa") == (
            "surveillance studies/no baloxavir exposure"
        )
        assert convert.expand_origin("Cell/BXA", "pa") == (
            "cell culture/substitution selected under baloxavir selection pressure"
        )
        assert convert.expand_origin("Clin/BXA", "pa") == (
            "clinical trial/substitution selected under baloxavir selection pressure"
        )

    def test_trailing_and_doubled_separators_are_dropped(self):
        assert convert.expand_origin("in vitro; recNA; Sur;", "avian-nai") == (
            "in vitro; recombinant NA; surveillance studies"
        )
        assert convert.expand_origin("Sur;; ", "human-nai") == "surveillance studies"

    def test_case_is_normalized(self):
        assert convert.expand_origin("In vitro; recNA; sur;", "avian-nai") == (
            "in vitro; recombinant NA; surveillance studies"
        )
        assert convert.expand_origin("In vivo/Zan", "human-nai") == (
            "in vivo/zanamivir used"
        )

    def test_unknown_token_kept_verbatim(self, capsys):
        assert convert.expand_origin("Sur; NewTok", "human-nai") == (
            "surveillance studies; NewTok"
        )
        assert "NewTok" in capsys.readouterr().err

    def test_real_human_row_shape(self):
        assert convert.expand_origin(
            "Clin/Ose, Clin/Ose+Zan; in vitro/Zan; RG; Clin/Sur", "human-nai"
        ) == (
            "clinical detection/oseltamivir used, "
            "clinical detection/oseltamivir used+zanamivir used; "
            "in vitro/zanamivir used; reverse genetics; "
            "clinical detection/surveillance studies"
        )


class TestBuildComment:
    def test_parts_joined_with_pipe(self):
        assert convert.build_comment(["a", "b", "c"]) == "a | b | c"

    def test_empty_parts_dropped(self):
        assert convert.build_comment(["", "b", ""]) == "b"


class TestRuleRegistryDedup:
    """Atomic rule deduplication and merging across source tables."""

    def make_registry(self):
        return convert.RuleRegistry()

    def test_identical_rules_deduplicate(self):
        registry = self.make_registry()
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "doi:10.1/a", "WHO human NAI", [],
                            subtype="A(H5N1)")
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "ri", 14, "", "doi:10.1/b", "WHO avian NAI", [],
                            subtype="A(H5N1)")
        rules, formulas = convert.finalize_rules(registry)
        assert len(rules) == 1
        assert len(formulas) == 0
        row = rules[0]
        assert row["phenotype"] == "hri"
        assert row["fold_ic50"] == "138"
        assert row["publication"] == "doi:10.1/a,doi:10.1/b"
        assert row["source"] == "WHO human NAI,WHO avian NAI"
        assert "Phenotypes in source" in row["comment"]
        assert "Fold-change values in source" in row["comment"]
        assert "Fold-change values in source (WHO avian NAI): 14" in row["comment"]
        assert "Fold-change values in source (WHO human NAI): 138" in row["comment"]

    def test_process_table_merges_subtypes_and_qualifies_conflicting_evidence(self):
        registry = self.make_registry()
        sequences = {"NC_026434": "A" * 274 + "H" + "A" * 100}
        non_migrated = []
        for subtype, fold_range in [
            ("A(H1N1)", "221–2846"),
            ("A(H1N1)pdm09", "321–2597"),
        ]:
            parsed = who_parsers.ParsedTable(
                footer_date="2026-01-01",
                rows=[
                    who_parsers.MarkerRow(
                        page=1,
                        subtype=subtype,
                        substitution="H275Y",
                        oseltamivir=f"HRI ({fold_range})",
                        origin=(
                            "Sur; RG" if subtype == "A(H1N1)"
                            else "Clin/Ose, Clin/Sur; in vitro"
                        ),
                    )
                ],
            )
            convert.process_table(
                "human-nai", parsed, {}, sequences, registry, non_migrated
            )

        rules, _ = convert.finalize_rules(registry)

        assert non_migrated == []
        assert len(rules) == 1
        assert rules[0]["source"] == "WHO human NAI"
        assert rules[0]["fold_ic50"] == "2846"
        assert rules[0]["comment"] == (
            "Origin of virus (A(H1N1)): surveillance studies; reverse genetics | "
            "Origin of virus (A(H1N1)pdm09): clinical detection/oseltamivir used, "
            "clinical detection/surveillance studies; in vitro | "
            "IC50 fold-change in source (A(H1N1)): 221–2846 | "
            "IC50 fold-change in source (A(H1N1)pdm09): 321–2597"
        )

    def test_merged_origin_comments_carry_source_labels(self):
        registry = self.make_registry()
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "doi:10.1/a", "WHO human NAI",
                            ["Origin of virus: surveillance studies"],
                            subtype="H5N1")
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "doi:10.1/b", "WHO avian NAI",
                            ["Origin of virus: reverse genetics"],
                            subtype="H5N1")
        rules, _ = convert.finalize_rules(registry)
        assert len(rules) == 1
        assert rules[0]["comment"] == (
            "Origin of virus (WHO human NAI): surveillance studies | "
            "Origin of virus (WHO avian NAI): reverse genetics"
        )

    def test_single_source_origin_comment_has_no_source_tag(self):
        registry = self.make_registry()
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "doi:10.1/a", "WHO avian NAI",
                            ["Origin of virus: in vitro; "
                             "reverse genetics/zanamivir used"],
                            subtype="H5N1")
        rules, _ = convert.finalize_rules(registry)
        assert len(rules) == 1
        assert rules[0]["comment"] == (
            "Origin of virus: in vitro; reverse genetics/zanamivir used"
        )

    def test_identical_origin_from_same_source_deduplicates(self):
        registry = self.make_registry()
        comment = "Origin of virus: reverse genetics"
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "doi:10.1/a", "WHO avian NAI",
                            [comment], subtype="H5N1")
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "doi:10.1/b", "WHO avian NAI",
                            [comment], subtype="H5N1")
        rules, _ = convert.finalize_rules(registry)
        assert len(rules) == 1
        assert rules[0]["comment"] == "Origin of virus: reverse genetics"

    def test_identical_origin_across_sources_stays_unqualified(self):
        registry = self.make_registry()
        comment = "Origin of virus: reverse genetics"
        registry.add_atomic(
            "NA", "EF619973", 255, "H", "Y", "oseltamivir",
            "hri", 138, "", "doi:10.1/a", "WHO human NAI", [comment],
            subtype="H5N1",
        )
        registry.add_atomic(
            "NA", "EF619973", 255, "H", "Y", "oseltamivir",
            "hri", 138, "", "doi:10.1/b", "WHO avian NAI", [comment],
            subtype="H5N1",
        )

        rules, _ = convert.finalize_rules(registry)

        assert rules[0]["comment"] == comment

    def test_unresolved_reference_note_carries_source_label(self):
        registry = self.make_registry()
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "", "WHO human NAI",
                            ["reference 5 unresolved"], subtype="H5N1")
        rules, _ = convert.finalize_rules(registry)
        assert rules[0]["comment"] == "reference 5 unresolved"

    def test_different_antivirals_stay_separate(self):
        registry = self.make_registry()
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "oseltamivir",
                            "hri", 138, "", "doi:10.1/a", "WHO human NAI", [],
                            subtype="")
        registry.add_atomic("NA", "EF619973", 255, "H", "Y", "zanamivir",
                            "ni", 3, "", "doi:10.1/a", "WHO human NAI", [],
                            subtype="")
        rules, _ = convert.finalize_rules(registry)
        assert len(rules) == 2
        assert {r["antiviral"] for r in rules} == {"oseltamivir", "zanamivir"}


class TestFinalizeRulesFormulas:
    """Formula group materialization: members, support rows, OR groups."""

    def _registry_with_combo(self):
        registry = convert.RuleRegistry()
        components = [
            {"member_key": ("NA", "NC_007368", "E119D"),
             "canonical": ("NA", "NC_007368", 119, "E", "D"),
             "or_group": None, "position": 119, "reference": "E",
             "mutation": "D"},
            {"member_key": ("NA", "NC_007368", "I222L"),
             "canonical": ("NA", "NC_007368", 222, "I", "L"),
             "or_group": None, "position": 222, "reference": "I",
             "mutation": "L"},
        ]
        registry.add_atomic("NA", "NC_007368", 119, "E", "D", "oseltamivir",
                            "ni", 1, "", "doi:10.1/a", "WHO human NAI", [],
                            subtype="")
        registry.add_combo("NA", "NC_007368",
                           tuple(c["member_key"] for c in components),
                           "oseltamivir", "hri", 799, "", "doi:10.1/b",
                           "WHO human NAI", [], components, subtype="")
        return registry, components

    def test_combo_emits_formula_and_reuses_singleton_member(self):
        registry, _ = self._registry_with_combo()
        rules, formulas = convert.finalize_rules(registry)
        assert len(formulas) == 1
        formula = formulas[0]
        assert formula["group_id"] == "G00001"
        assert formula["expression"] == "(M00001 AND M00002)"
        assert formula["phenotype"] == "hri"
        # The singleton atomic rule (E119D) is reused as a member: it carries
        # the member_id and no blank-antiviral support row is created.
        singleton = [r for r in rules
                     if (r["position"], r["mutation"]) == (119, "D")]
        assert len(singleton) == 1
        assert singleton[0]["member_id"] == "M00001"
        assert singleton[0]["antiviral"] == "oseltamivir"
        support = [r for r in rules if r["antiviral"] == ""]
        assert len(support) == 1
        assert support[0]["member_id"] == "M00002"

    def test_mixed_alleles_form_or_group(self):
        registry = convert.RuleRegistry()
        components = [
            {"member_key": ("NA", "X", "D151E"),
             "canonical": ("NA", "X", 151, "D", "E"),
             "or_group": 0, "position": 151, "reference": "D",
             "mutation": "E"},
            {"member_key": ("NA", "X", "D151N"),
             "canonical": ("NA", "X", 151, "D", "N"),
             "or_group": 0, "position": 151, "reference": "D",
             "mutation": "N"},
            {"member_key": ("NA", "X", "H275Y"),
             "canonical": ("NA", "X", 275, "H", "Y"),
             "or_group": None, "position": 275, "reference": "H",
             "mutation": "Y"},
        ]
        registry.add_combo("NA", "X",
                           tuple(c["member_key"] for c in components),
                           "oseltamivir", "hri", 799, "", "doi:10.1/b",
                           "WHO human NAI", [], components, subtype="")
        _, formulas = convert.finalize_rules(registry)
        assert formulas[0]["expression"] == "((M00001 OR M00002) AND M00003)"

    def test_formula_groups_from_different_subtypes_merge_with_provenance(self):
        registry = convert.RuleRegistry()
        for subtype in ("A(H1N1)", "A(H1N1)pdm09"):
            components = [
                {
                    "member_key": ("NA", "NC_026434", "E119D"),
                    "canonical": ("NA", "NC_026434", 119, "E", "D"),
                    "or_group": None, "position": 119, "reference": "E",
                    "mutation": "D",
                },
                {
                    "member_key": ("NA", "NC_026434", "I222L"),
                    "canonical": ("NA", "NC_026434", 222, "I", "L"),
                    "or_group": None, "position": 222, "reference": "I",
                    "mutation": "L",
                },
            ]
            registry.add_combo(
                "NA", "NC_026434",
                tuple(component["member_key"] for component in components),
                "oseltamivir", "hri", 799, "", "doi:10.1/a",
                "WHO human NAI",
                [
                    "Origin of virus: surveillance studies"
                    if subtype == "A(H1N1)"
                    else "Origin of virus: clinical detection"
                ],
                components,
                subtype=subtype,
            )

        rules, formulas = convert.finalize_rules(registry)

        assert len(formulas) == 1
        assert len(rules) == 2
        assert formulas[0]["source"] == "WHO human NAI"
        assert formulas[0]["comment"] == (
            "Origin of virus (A(H1N1)): surveillance studies | "
            "Origin of virus (A(H1N1)pdm09): clinical detection"
        )


class TestCitationListParsing:
    """Citation list extraction from WHO PDF page text."""

    def test_sequential_citations(self):
        text = "\n".join([
            "1. Author A. Title one. Journal. 2020;12:1-10.",
            "2. Author B. Title two. Journal. 2021;13:2-20.",
            "3. Author C. Title three. Journal. 2022;14:3-30.",
        ])
        entries = resolve_references.parse_citation_list([text])
        assert sorted(entries) == [1, 2, 3]
        assert "Title two" in entries[2]

    def test_line_wrapped_citation(self):
        text = "\n".join([
            "1. Author A. A very long title that wraps across",
            "two lines here. Journal. 2020;12:1-10.",
            "2. Author B. Title two. Journal. 2021;13:2-20.",
        ])
        entries = resolve_references.parse_citation_list([text])
        assert sorted(entries) == [1, 2]
        assert "two lines here" in entries[1]

    def test_dash_page_number_not_treated_as_citation(self):
        text = "\n".join([
            "1. Author A. Title one. Journal. 2020;12:1019\u2013",
            "23. doi:10.1/aaa",
            "2. Author B. Title two. Journal. 2021;13:2-20.",
        ])
        entries = resolve_references.parse_citation_list([text])
        assert sorted(entries) == [1, 2]

    def test_glued_page_footer_is_stripped(self):
        text = "\n".join([
            "1. Author A. Title one. https://doi.org/10.1/aaa. "
            "Last updated on 17 July 2026 Page 5 of 7",
            "2. Author B. Title two. Journal. 2021;13:2-20.",
        ])
        entries = resolve_references.parse_citation_list([text])
        assert "Last updated" not in entries[1]
        assert "Page 5 of 7" not in entries[1]

    def test_doi_extraction_with_trailing_period(self):
        doi = resolve_references.citation_doi(
            "Some title. https://doi.org/10.1371/journal.ppat.1000933.")
        assert doi == "10.1371/journal.ppat.1000933"

    def test_doi_wrapped_at_hyphen_is_rejoined(self):
        # Euro Surveill. DOIs wrap after the journal-prefix hyphen; the
        # primary regex stops at the wrap whitespace and must not win.
        doi = resolve_references.citation_doi(
            "Takashita E, et al. Euro Surveill. 2019;24(12):1900170. "
            "https://doi.org/10.2807/1560- 7917.ES.2019.24.12.1900170.")
        assert doi == "10.2807/1560-7917.ES.2019.24.12.1900170"

    def test_doi_without_scheme_uses_repair_path(self):
        # Without an http(s) scheme the primary regex cannot match; the
        # repair path validates the text after "doi.org/" as a DOI.
        doi = resolve_references.citation_doi(
            "Some title. Available from doi.org/10.1016/j.virusres.2019.03.019.")
        assert doi == "10.1016/j.virusres.2019.03.019"

    def test_doi_repair_requires_prefix_continuation(self):
        # A truncated primary match must only be extended by a candidate that
        # continues it, never replaced by unrelated text.
        doi = resolve_references.citation_doi(
            "Title. https://doi.org/10.1016/j.antiviral.2020.104- (see also "
            "doi.org/10.9999/unrelated)")
        assert doi == "10.1016/j.antiviral.2020.104-"


class TestWhoParsersHelpers:
    """Parser helper functions."""

    def test_group_words_into_lines_by_y(self):
        words = [
            {"text": "H275Y", "x0": 10, "x1": 50, "top": 10, "bottom": 20},
            {"text": "NI", "x0": 100, "x1": 120, "top": 12, "bottom": 22},
            {"text": "V116A", "x0": 10, "x1": 50, "top": 40, "bottom": 50},
        ]
        lines = who_parsers._group_words_into_lines(words)
        assert len(lines) == 2
        assert [w["text"] for w in lines[0]] == ["H275Y", "NI"]

    def test_line_text_joins_with_spaces(self):
        line = [{"text": "Del", "x0": 0, "x1": 10},
                {"text": "245\u2013248", "x0": 12, "x1": 40}]
        assert who_parsers._line_text(line) == "Del 245\u2013248"
