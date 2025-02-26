from textwrap import dedent

from tola.assembly.terminal_table import (
    Table,
    TableCell,
    TableHeader,
    TableLine,
    TableRow,
)


def test_simple_table():
    tbl = Table(
        header=TableHeader(
            [
                TableCell([TableLine("Hap1")]),
                TableCell([TableLine("Hap2")]),
                TableCell([TableLine("Hap3")]),
            ]
        ),
        rows=[
            TableRow(
                [
                    TableCell(
                        [
                            TableLine("Scaffold1"),
                            TableLine("25_111_000 bp"),
                        ]
                    ),
                    TableCell(),
                    TableCell(
                        [
                            TableLine("Scaffold2"),
                            TableLine("24_000_000 bp"),
                            TableLine("Singleton"),
                        ]
                    ),
                ]
            ),
            TableRow(
                [
                    TableCell(
                        [
                            TableLine("Scaffold3"),
                            TableLine("19_330_000 bp"),
                        ]
                    ),
                    TableCell(
                        [
                            TableLine("Scaffold4"),
                            TableLine("18_980_000 bp"),
                        ]
                    ),
                ]
            ),
        ],
    )
    expected = dedent("""
        ┌───────────────┬───────────────┬───────────────┐
        │     Hap1      │     Hap2      │     Hap3      │
        ├───────────────┼───────────────┼───────────────┤
        │   Scaffold1   │               │   Scaffold2   │
        │ 25_111_000 bp │               │ 24_000_000 bp │
        │               │               │   Singleton   │
        ├───────────────┼───────────────┼───────────────┤
        │   Scaffold3   │   Scaffold4   │               │
        │ 19_330_000 bp │ 18_980_000 bp │               │
        └───────────────┴───────────────┴───────────────┘
         """).lstrip()
    assert tbl.render() == expected


def test_no_header_table():
    tbl = Table(
        rows=[
            TableRow(
                [
                    TableCell(
                        [
                            TableLine("Scaffold1"),
                            TableLine("25_111_000 bp"),
                        ]
                    ),
                    TableCell(
                        [
                            TableLine("Scaffold2"),
                            TableLine("24_000_000 bp"),
                            TableLine("Singleton"),
                        ]
                    ),
                ]
            ),
        ],
    )
    expected = dedent("""
        ┌───────────────┬───────────────┐
        │   Scaffold1   │   Scaffold2   │
        │ 25_111_000 bp │ 24_000_000 bp │
        │               │   Singleton   │
        └───────────────┴───────────────┘
        """).lstrip()
    assert tbl.render() == expected
