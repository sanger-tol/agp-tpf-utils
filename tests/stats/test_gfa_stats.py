from io import StringIO
from pathlib import Path

import pytest

from tola.assembly.gfa_stats import GfaStats, GfaStatsError


@pytest.fixture
def genome_file():
    data_dir = Path(__file__).parent.parent / "fasta"
    return data_dir / "test_other.fa"


def test_gfasta_stats(genome_file):
    pass


def test_parse_gfa_output(genome_file):
    gfa_file = Path(__file__).parent / "gfa_stats.txt"
    io = StringIO(gfa_file.read_text())

    gfa = GfaStats(Path("none"))
    gfa.parse(io)

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


if __name__ == "__main__":
    gfa_file = Path(__file__).parent / "gfa_stats.txt"
    io = StringIO(gfa_file.read_text())

    gfa = GfaStats(Path("none"))
    gfa.parse(io)
    print(f"{gfa!r}")
