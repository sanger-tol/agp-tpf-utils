from textwrap import dedent

from tola.assembly.terminal_table import (
    CellLine,
    Table,
    TableCell,
    TableHeader,
    TableRow,
)


def test_simple_table():
    tbl = Table(
        header=TableHeader(
            [
                TableCell([CellLine("Hap1")]),
                TableCell([CellLine("Hap2")]),
                TableCell([CellLine("Hap3")]),
            ]
        ),
        rows=[
            TableRow(
                [
                    TableCell(
                        [
                            CellLine("Scaffold1"),
                            CellLine("25_111_000 bp"),
                        ]
                    ),
                    TableCell(),
                    TableCell(
                        [
                            CellLine("Scaffold2"),
                            CellLine("24_000_000 bp"),
                            CellLine("Singleton"),
                        ]
                    ),
                ]
            ),
            TableRow(
                [
                    TableCell(
                        [
                            CellLine("Scaffold3"),
                            CellLine("19_330_000 bp"),
                        ]
                    ),
                    TableCell(
                        [
                            CellLine("Scaffold4"),
                            CellLine("18_980_000 bp"),
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
                            CellLine("Scaffold1"),
                            CellLine("25_111_000 bp"),
                        ]
                    ),
                    TableCell(
                        [
                            CellLine("Scaffold2"),
                            CellLine("24_000_000 bp"),
                            CellLine("Singleton"),
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
