import re
import sys
from io import StringIO
from pathlib import Path
from shutil import which

import pytest

from tola.assembly.gfa_stats import GfaStats, GfaStatsError, VgpStats


@pytest.fixture
def genome_file():
    data_dir = Path(__file__).parent.parent / "fasta"
    return data_dir / "test_other.fa"


@pytest.fixture
def gfa_after_file():
    return Path(__file__).parent / "gfa_stats_after.txt"


@pytest.fixture
def gfa_before_file():
    return Path(__file__).parent / "gfa_stats_after.txt"


def test_gfasta_stats(genome_file):
    if not which("gfastats"):
        pytest.skip("No gfastats executable")  # ty: ignore[too-many-positional-arguments]

    none_such = "no_such_dir/no_such_asm.fa.gz"
    gfa1 = GfaStats(Path(none_such))
    with pytest.raises(
        GfaStatsError,
        match=re.compile(
            rf'Error running.+"{none_such}".+file does not exist',
            flags=re.DOTALL,
        ),
    ):
        gfa1.run()

    gfa2 = GfaStats(genome_file).run()
    assert len(list(gfa2.stats())) == 92
    io = StringIO()
    for stat in gfa2.stats():
        io.write(f"{stat.as_dict()}")
    # ruff: disable[E501]
    assert io.getvalue() == (
        "{'gfa_name': '# scaffolds', 'col_name': 'scaffolds', 'category': 'scaffold', 'measure': 'count', 'value': 3}"
        "{'gfa_name': 'Total scaffold length', 'col_name': 'scaffold_length', 'category': 'scaffold', 'measure': 'total', 'value': 1061}"
        "{'gfa_name': 'Average scaffold length', 'col_name': 'scaffold_length_mean', 'category': 'scaffold', 'measure': 'mean_length', 'value': 353.67}"
        "{'gfa_name': 'Scaffold N50', 'col_name': 'scaffold_n50', 'category': 'scaffold', 'measure': 'n50', 'value': 609}"
        "{'gfa_name': 'Scaffold auN', 'col_name': 'scaffold_aun', 'category': 'scaffold', 'measure': 'aun', 'value': 454.05}"
        "{'gfa_name': 'Scaffold L50', 'col_name': 'scaffold_l50', 'category': 'scaffold', 'measure': 'l50', 'value': 1}"
        "{'gfa_name': 'Largest scaffold', 'col_name': 'scaffold_longest', 'category': 'scaffold', 'measure': 'longest_length', 'value': 609}"
        "{'gfa_name': 'Smallest scaffold', 'col_name': 'scaffold_shortest', 'category': 'scaffold', 'measure': 'shortest_length', 'value': 160}"
        "{'gfa_name': '# contigs', 'col_name': 'contigs', 'category': 'contig', 'measure': 'count', 'value': 13}"
        "{'gfa_name': 'Total contig length', 'col_name': 'contig_length', 'category': 'contig', 'measure': 'total', 'value': 1051}"
        "{'gfa_name': 'Average contig length', 'col_name': 'contig_length_mean', 'category': 'contig', 'measure': 'mean_length', 'value': 80.85}"
        "{'gfa_name': 'Contig N50', 'col_name': 'contig_n50', 'category': 'contig', 'measure': 'n50', 'value': 159}"
        "{'gfa_name': 'Contig auN', 'col_name': 'contig_aun', 'category': 'contig', 'measure': 'aun', 'value': 128.01}"
        "{'gfa_name': 'Contig L50', 'col_name': 'contig_l50', 'category': 'contig', 'measure': 'l50', 'value': 3}"
        "{'gfa_name': 'Largest contig', 'col_name': 'contig_longest', 'category': 'contig', 'measure': 'longest_length', 'value': 216}"
        "{'gfa_name': 'Smallest contig', 'col_name': 'contig_shortest', 'category': 'contig', 'measure': 'shortest_length', 'value': 15}"
        "{'gfa_name': '# gaps in scaffolds', 'col_name': 'gaps', 'category': 'gap', 'measure': 'count', 'value': 10}"
        "{'gfa_name': 'Total gap length in scaffolds', 'col_name': 'gap_length_in_scaffolds', 'category': 'gap', 'measure': 'total', 'value': 10}"
        "{'gfa_name': 'Average gap length in scaffolds', 'col_name': 'gap_length_in_scaffolds_mean', 'category': 'gap', 'measure': 'mean_length', 'value': 1.0}"
        "{'gfa_name': 'Gap N50 in scaffolds', 'col_name': 'gap_n50', 'category': 'gap', 'measure': 'n50', 'value': 1}"
        "{'gfa_name': 'Gap auN in scaffolds', 'col_name': 'gap_aun', 'category': 'gap', 'measure': 'aun', 'value': 1.0}"
        "{'gfa_name': 'Gap L50 in scaffolds', 'col_name': 'gap_l50', 'category': 'gap', 'measure': 'l50', 'value': 5}"
        "{'gfa_name': 'Largest gap in scaffolds', 'col_name': 'gap_in_scaffolds_longest', 'category': 'gap', 'measure': 'longest_length', 'value': 1}"
        "{'gfa_name': 'Smallest gap in scaffolds', 'col_name': 'gap_in_scaffolds_shortest', 'category': 'gap', 'measure': 'shortest_length', 'value': 1}"
        "{'gfa_name': 'Base composition (A:C:G:T)', 'col_name': 'base_counts_a_c_g_t', 'category': 'contig', 'measure': 'composition', 'value': [295, 207, 209, 325]}"
        "{'gfa_name': 'GC content %', 'col_name': 'gc_content_percent', 'category': 'contig', 'measure': 'gc_content', 'value': 40.15}"
        "{'gfa_name': '# soft-masked bases', 'col_name': 'soft-masked', 'category': 'soft-masked', 'measure': 'count', 'value': 1046}"
        "{'gfa_name': '# segments', 'col_name': 'segments', 'category': 'segment', 'measure': 'count', 'value': 13}"
        "{'gfa_name': 'Total segment length', 'col_name': 'segment_length', 'category': 'segment', 'measure': 'total', 'value': 1036}"
        "{'gfa_name': 'Average segment length', 'col_name': 'segment_length_mean', 'category': 'segment', 'measure': 'mean_length', 'value': 79.69}"
        "{'gfa_name': '# gaps', 'col_name': 'gaps', 'category': 'gap', 'measure': 'count', 'value': 10}"
        "{'gfa_name': '# paths', 'col_name': 'paths', 'category': 'path', 'measure': 'count', 'value': 3}"
        "{'gfa_name': 'Scaffold N10', 'col_name': 'scaffold_n10', 'category': 'scaffold', 'measure': 'n10', 'value': 609}"
        "{'gfa_name': 'Scaffold N20', 'col_name': 'scaffold_n20', 'category': 'scaffold', 'measure': 'n20', 'value': 609}"
        "{'gfa_name': 'Scaffold N30', 'col_name': 'scaffold_n30', 'category': 'scaffold', 'measure': 'n30', 'value': 609}"
        "{'gfa_name': 'Scaffold N40', 'col_name': 'scaffold_n40', 'category': 'scaffold', 'measure': 'n40', 'value': 609}"
        "{'gfa_name': 'Scaffold N50', 'col_name': 'scaffold_n50', 'category': 'scaffold', 'measure': 'n50', 'value': 609}"
        "{'gfa_name': 'Scaffold N60', 'col_name': 'scaffold_n60', 'category': 'scaffold', 'measure': 'n60', 'value': 292}"
        "{'gfa_name': 'Scaffold N70', 'col_name': 'scaffold_n70', 'category': 'scaffold', 'measure': 'n70', 'value': 292}"
        "{'gfa_name': 'Scaffold N80', 'col_name': 'scaffold_n80', 'category': 'scaffold', 'measure': 'n80', 'value': 292}"
        "{'gfa_name': 'Scaffold N90', 'col_name': 'scaffold_n90', 'category': 'scaffold', 'measure': 'n90', 'value': 160}"
        "{'gfa_name': 'Scaffold N100', 'col_name': 'scaffold_n100', 'category': 'scaffold', 'measure': 'n100', 'value': 160}"
        "{'gfa_name': 'Scaffold L10', 'col_name': 'scaffold_l10', 'category': 'scaffold', 'measure': 'l10', 'value': 1}"
        "{'gfa_name': 'Scaffold L20', 'col_name': 'scaffold_l20', 'category': 'scaffold', 'measure': 'l20', 'value': 1}"
        "{'gfa_name': 'Scaffold L30', 'col_name': 'scaffold_l30', 'category': 'scaffold', 'measure': 'l30', 'value': 1}"
        "{'gfa_name': 'Scaffold L40', 'col_name': 'scaffold_l40', 'category': 'scaffold', 'measure': 'l40', 'value': 1}"
        "{'gfa_name': 'Scaffold L50', 'col_name': 'scaffold_l50', 'category': 'scaffold', 'measure': 'l50', 'value': 1}"
        "{'gfa_name': 'Scaffold L60', 'col_name': 'scaffold_l60', 'category': 'scaffold', 'measure': 'l60', 'value': 2}"
        "{'gfa_name': 'Scaffold L70', 'col_name': 'scaffold_l70', 'category': 'scaffold', 'measure': 'l70', 'value': 2}"
        "{'gfa_name': 'Scaffold L80', 'col_name': 'scaffold_l80', 'category': 'scaffold', 'measure': 'l80', 'value': 2}"
        "{'gfa_name': 'Scaffold L90', 'col_name': 'scaffold_l90', 'category': 'scaffold', 'measure': 'l90', 'value': 3}"
        "{'gfa_name': 'Scaffold L100', 'col_name': 'scaffold_l100', 'category': 'scaffold', 'measure': 'l100', 'value': 3}"
        "{'gfa_name': 'Contig N10', 'col_name': 'contig_n10', 'category': 'contig', 'measure': 'n10', 'value': 216}"
        "{'gfa_name': 'Contig N20', 'col_name': 'contig_n20', 'category': 'contig', 'measure': 'n20', 'value': 216}"
        "{'gfa_name': 'Contig N30', 'col_name': 'contig_n30', 'category': 'contig', 'measure': 'n30', 'value': 160}"
        "{'gfa_name': 'Contig N40', 'col_name': 'contig_n40', 'category': 'contig', 'measure': 'n40', 'value': 159}"
        "{'gfa_name': 'Contig N50', 'col_name': 'contig_n50', 'category': 'contig', 'measure': 'n50', 'value': 159}"
        "{'gfa_name': 'Contig N60', 'col_name': 'contig_n60', 'category': 'contig', 'measure': 'n60', 'value': 121}"
        "{'gfa_name': 'Contig N70', 'col_name': 'contig_n70', 'category': 'contig', 'measure': 'n70', 'value': 102}"
        "{'gfa_name': 'Contig N80', 'col_name': 'contig_n80', 'category': 'contig', 'measure': 'n80', 'value': 50}"
        "{'gfa_name': 'Contig N90', 'col_name': 'contig_n90', 'category': 'contig', 'measure': 'n90', 'value': 43}"
        "{'gfa_name': 'Contig N100', 'col_name': 'contig_n100', 'category': 'contig', 'measure': 'n100', 'value': 15}"
        "{'gfa_name': 'Contig L10', 'col_name': 'contig_l10', 'category': 'contig', 'measure': 'l10', 'value': 1}"
        "{'gfa_name': 'Contig L20', 'col_name': 'contig_l20', 'category': 'contig', 'measure': 'l20', 'value': 1}"
        "{'gfa_name': 'Contig L30', 'col_name': 'contig_l30', 'category': 'contig', 'measure': 'l30', 'value': 2}"
        "{'gfa_name': 'Contig L40', 'col_name': 'contig_l40', 'category': 'contig', 'measure': 'l40', 'value': 3}"
        "{'gfa_name': 'Contig L50', 'col_name': 'contig_l50', 'category': 'contig', 'measure': 'l50', 'value': 3}"
        "{'gfa_name': 'Contig L60', 'col_name': 'contig_l60', 'category': 'contig', 'measure': 'l60', 'value': 4}"
        "{'gfa_name': 'Contig L70', 'col_name': 'contig_l70', 'category': 'contig', 'measure': 'l70', 'value': 5}"
        "{'gfa_name': 'Contig L80', 'col_name': 'contig_l80', 'category': 'contig', 'measure': 'l80', 'value': 7}"
        "{'gfa_name': 'Contig L90', 'col_name': 'contig_l90', 'category': 'contig', 'measure': 'l90', 'value': 9}"
        "{'gfa_name': 'Contig L100', 'col_name': 'contig_l100', 'category': 'contig', 'measure': 'l100', 'value': 13}"
        "{'gfa_name': 'Gap N10', 'col_name': 'gap_n10', 'category': 'gap', 'measure': 'n10', 'value': 1}"
        "{'gfa_name': 'Gap N20', 'col_name': 'gap_n20', 'category': 'gap', 'measure': 'n20', 'value': 1}"
        "{'gfa_name': 'Gap N30', 'col_name': 'gap_n30', 'category': 'gap', 'measure': 'n30', 'value': 1}"
        "{'gfa_name': 'Gap N40', 'col_name': 'gap_n40', 'category': 'gap', 'measure': 'n40', 'value': 1}"
        "{'gfa_name': 'Gap N50', 'col_name': 'gap_n50', 'category': 'gap', 'measure': 'n50', 'value': 1}"
        "{'gfa_name': 'Gap N60', 'col_name': 'gap_n60', 'category': 'gap', 'measure': 'n60', 'value': 1}"
        "{'gfa_name': 'Gap N70', 'col_name': 'gap_n70', 'category': 'gap', 'measure': 'n70', 'value': 1}"
        "{'gfa_name': 'Gap N80', 'col_name': 'gap_n80', 'category': 'gap', 'measure': 'n80', 'value': 1}"
        "{'gfa_name': 'Gap N90', 'col_name': 'gap_n90', 'category': 'gap', 'measure': 'n90', 'value': 1}"
        "{'gfa_name': 'Gap N100', 'col_name': 'gap_n100', 'category': 'gap', 'measure': 'n100', 'value': 1}"
        "{'gfa_name': 'Gap L10', 'col_name': 'gap_l10', 'category': 'gap', 'measure': 'l10', 'value': 1}"
        "{'gfa_name': 'Gap L20', 'col_name': 'gap_l20', 'category': 'gap', 'measure': 'l20', 'value': 2}"
        "{'gfa_name': 'Gap L30', 'col_name': 'gap_l30', 'category': 'gap', 'measure': 'l30', 'value': 3}"
        "{'gfa_name': 'Gap L40', 'col_name': 'gap_l40', 'category': 'gap', 'measure': 'l40', 'value': 4}"
        "{'gfa_name': 'Gap L50', 'col_name': 'gap_l50', 'category': 'gap', 'measure': 'l50', 'value': 5}"
        "{'gfa_name': 'Gap L60', 'col_name': 'gap_l60', 'category': 'gap', 'measure': 'l60', 'value': 6}"
        "{'gfa_name': 'Gap L70', 'col_name': 'gap_l70', 'category': 'gap', 'measure': 'l70', 'value': 7}"
        "{'gfa_name': 'Gap L80', 'col_name': 'gap_l80', 'category': 'gap', 'measure': 'l80', 'value': 8}"
        "{'gfa_name': 'Gap L90', 'col_name': 'gap_l90', 'category': 'gap', 'measure': 'l90', 'value': 9}"
        "{'gfa_name': 'Gap L100', 'col_name': 'gap_l100', 'category': 'gap', 'measure': 'l100', 'value': 10}"
    )
    # ruff: enable[E501]


def test_all_stat_fields_are_filled(gfa_after_file):
    gfa = GfaStats(Path("test_all_stats_populated"))
    gfa.parse(gfa_after_file.open())
    for stat in gfa.stats():
        for attr in (
            "gfa_name",
            "col_name",
            "category",
            "measure",
            "value",
        ):
            assert hasattr(stat, attr) and getattr(stat, attr) is not None


def test_parse_gfa_output(gfa_after_file):
    gfa = GfaStats(Path("working/curation/nxChoDudi2.2.primary.curated.fa.gz"))
    gfa.parse(gfa_after_file.open())

    assert gfa.short_stats_string() == (
        "# scaffolds\t296\n"
        "Total scaffold length\t139,596,165\n"
        "Average scaffold length\t471,608.67\n"
        "Scaffold N50\t16,442,996\n"
        "Scaffold auN\t16,337,482.47\n"
        "Scaffold L50\t4\n"
        "Largest scaffold\t23,458,917\n"
        "Smallest scaffold\t1,000\n"
        "# contigs\t322\n"
        "Total contig length\t139,593,065\n"
        "Average contig length\t433,518.84\n"
        "Contig N50\t4,166,186\n"
        "Contig auN\t6,365,325.74\n"
        "Contig L50\t9\n"
        "Largest contig\t16,442,996\n"
        "Smallest contig\t1,000\n"
        "# gaps in scaffolds\t26\n"
        "Total gap length in scaffolds\t3,100\n"
        "Average gap length in scaffolds\t119.23\n"
        "Gap N50 in scaffolds\t100\n"
        "Gap auN in scaffolds\t132.26\n"
        "Gap L50 in scaffolds\t11\n"
        "Largest gap in scaffolds\t200\n"
        "Smallest gap in scaffolds\t100\n"
        "Base composition (A:C:G:T)\t44548621:25208559:25300801:44535084\n"
        "GC content %\t36.18\n"
        "# soft-masked bases\t0\n"
        "# segments\t322\n"
        "Total segment length\t139,593,065\n"
        "Average segment length\t433,518.84\n"
        "# gaps\t26\n"
        "# paths\t296\n"
    )

    assert gfa.stats_string() == (
        "# scaffolds\t296\n"
        "Total scaffold length\t139,596,165\n"
        "Average scaffold length\t471,608.67\n"
        "Scaffold N50\t16,442,996\n"
        "Scaffold auN\t16,337,482.47\n"
        "Scaffold L50\t4\n"
        "Largest scaffold\t23,458,917\n"
        "Smallest scaffold\t1,000\n"
        "# contigs\t322\n"
        "Total contig length\t139,593,065\n"
        "Average contig length\t433,518.84\n"
        "Contig N50\t4,166,186\n"
        "Contig auN\t6,365,325.74\n"
        "Contig L50\t9\n"
        "Largest contig\t16,442,996\n"
        "Smallest contig\t1,000\n"
        "# gaps in scaffolds\t26\n"
        "Total gap length in scaffolds\t3,100\n"
        "Average gap length in scaffolds\t119.23\n"
        "Gap N50 in scaffolds\t100\n"
        "Gap auN in scaffolds\t132.26\n"
        "Gap L50 in scaffolds\t11\n"
        "Largest gap in scaffolds\t200\n"
        "Smallest gap in scaffolds\t100\n"
        "Base composition (A:C:G:T)\t44548621:25208559:25300801:44535084\n"
        "GC content %\t36.18\n"
        "# soft-masked bases\t0\n"
        "# segments\t322\n"
        "Total segment length\t139,593,065\n"
        "Average segment length\t433,518.84\n"
        "# gaps\t26\n"
        "# paths\t296\n"
        "Scaffold N10\t23,458,917\n"
        "Scaffold N20\t18,591,642\n"
        "Scaffold N30\t18,591,642\n"
        "Scaffold N40\t18,554,823\n"
        "Scaffold N50\t16,442,996\n"
        "Scaffold N60\t16,131,304\n"
        "Scaffold N70\t15,994,921\n"
        "Scaffold N80\t15,839,936\n"
        "Scaffold N90\t1,375,489\n"
        "Scaffold N100\t1,000\n"
        "Scaffold L10\t1\n"
        "Scaffold L20\t2\n"
        "Scaffold L30\t2\n"
        "Scaffold L40\t3\n"
        "Scaffold L50\t4\n"
        "Scaffold L60\t5\n"
        "Scaffold L70\t6\n"
        "Scaffold L80\t7\n"
        "Scaffold L90\t8\n"
        "Scaffold L100\t296\n"
        "Contig N10\t16,442,996\n"
        "Contig N20\t12,188,499\n"
        "Contig N30\t6,827,089\n"
        "Contig N40\t5,877,668\n"
        "Contig N50\t4,166,186\n"
        "Contig N60\t3,689,304\n"
        "Contig N70\t3,054,886\n"
        "Contig N80\t1,886,264\n"
        "Contig N90\t716,691\n"
        "Contig N100\t1,000\n"
        "Contig L10\t1\n"
        "Contig L20\t2\n"
        "Contig L30\t4\n"
        "Contig L40\t6\n"
        "Contig L50\t9\n"
        "Contig L60\t12\n"
        "Contig L70\t16\n"
        "Contig L80\t22\n"
        "Contig L90\t31\n"
        "Contig L100\t322\n"
        "Gap N10\t200\n"
        "Gap N20\t200\n"
        "Gap N30\t200\n"
        "Gap N40\t100\n"
        "Gap N50\t100\n"
        "Gap N60\t100\n"
        "Gap N70\t100\n"
        "Gap N80\t100\n"
        "Gap N90\t100\n"
        "Gap N100\t100\n"
        "Gap L10\t2\n"
        "Gap L20\t4\n"
        "Gap L30\t5\n"
        "Gap L40\t8\n"
        "Gap L50\t11\n"
        "Gap L60\t14\n"
        "Gap L70\t17\n"
        "Gap L80\t20\n"
        "Gap L90\t23\n"
        "Gap L100\t26\n"
    )


def test_vgp_stats(gfa_before_file, gfa_after_file):
    before = GfaStats(Path("draft"))
    after = GfaStats(Path("curated"))

    before.parse(gfa_before_file.open())
    after.parse(gfa_after_file.open())

    vgp_stats = VgpStats()
    vgp_stats.load_before_stats(before)
    vgp_stats.load_after_stats(after)

    assert str(vgp_stats) == (
        "scaffolds\n"
        "total 139596165 139596165\n"
        "count 296 296\n"
        "N50 16442996 16442996\n"
        "L50 4 4\n"
        "N90 1375489 1375489\n"
        "L90 8 8\n"
        "\n"
        "contigs\n"
        "total 139593065 139593065\n"
        "count 322 322\n"
        "N50 4166186 4166186\n"
        "L50 9 9\n"
        "N90 716691 716691\n"
        "L90 31 31\n"
    )


if __name__ == "__main__":
    file = Path(__file__).parent / "gfa_stats.txt"
    gfa = GfaStats(Path("none"))
    gfa.parse(file.open())
    sys.stdout.write(gfa.as_ndjson())
