from textwrap import dedent

from tola.assembly.terminal_table import (
    CellLine,
    TableCell,
    TableHeader,
    TableRow,
    TerminalTable,
)


def test_simple_table():
    tbl = TerminalTable(
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
    tbl = TerminalTable(
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


def test_contiguous_ranges():
    # 012 56
    assert list(
        TerminalTable.contiguous_ranges(
            [0, 1, 6],
            7,
        )
    ) == [
        (0, 2),
        (5, 6),
    ]

    # 012345
    assert list(
        TerminalTable.contiguous_ranges(
            [0, 1, 2, 4],
            7,
        )
    ) == [
        (0, 5),
    ]

    # 0123 5678
    assert list(
        TerminalTable.contiguous_ranges(
            [0, 8],
            9,
            3,
        )
    ) == [
        (0, 3),
        (5, 8),
    ]
