import io

import click


class TableLine:
    __slots__ = "text", "render"

    def __init__(self, text: tuple[str], render=None):
        self.text = text if isinstance(text, tuple) else (text,)
        self.render = render

    def length(self):
        return sum(len(x) for x in self.text)

    def format(self, width):
        txt = rdr(self.text) if (rdr := self.render) else "".join(self.text)
        pad = width - self.length()
        l_pad = " " * (pad // 2)
        r_pad = " " * (pad // 2 + pad % 2)
        return l_pad + txt + r_pad


class TableCell:
    __slots__ = "lines"

    def __init__(self, lines: list[TableLine] = None):
        self.lines = lines if lines else []

    def add_line(self, line: TableLine):
        self.lines.append(line)

    def max_line_length(self):
        if lns := self.lines:
            return max(x.length() for x in lns)
        else:
            return 0

    def line_count(self):
        return len(self.lines)

    def formatted_line(self, width, index):
        if index < len(self.lines):
            return self.lines[index].format(width)
        else:
            return " " * width


class TableRow:
    __slots__ = "cells"

    def __init__(self, cells: list[TableCell] = None):
        self.cells = cells if cells else []

    def add_cell(self, cell: TableCell):
        self.cells.append(cell)

    def get_cell(self, index):
        if index < len(self.cells):
            return self.cells[index]
        else:
            return None

    def max_line_length(self):
        return max(c.max_line_length() for c in self.cells)

    def column_count(self):
        return len(self.cells)

    def line_count(self):
        return max(c.line_count() for c in self.cells)


class TableHeader(TableRow):
    pass


class Table:
    __slots__ = "header", "rows", "buffer"

    def __init__(self, header: TableHeader = None, rows: list[TableRow] = None):
        self.header = header
        self.rows = rows if rows else []
        self.buffer = io.StringIO()

    def add_row(self, row: TableRow):
        self.rows.append(row)

    def render(self):
        # Empty the buffer
        self.buffer.seek(0)
        self.buffer.truncate(0)

        n_cols = self.column_count()
        width = self.max_line_length()
        pad_width = width + 2

        self.buffer.write("┌" + "┬".join(n_cols * ["─" * pad_width]) + "┐\n")

        ruler = "├" + "┼".join(n_cols * ["─" * pad_width]) + "┤\n"

        if hdr := self.header:
            self.render_row(hdr, n_cols, width)

        for row in self.rows:
            self.buffer.write(ruler)
            self.render_row(row, n_cols, width)

        self.buffer.write("└" + "┴".join(n_cols * ["─" * pad_width]) + "┘\n")

        return self.buffer.getvalue()

    def render_row(self, row, n_cols, width):
        for line_i in range(row.line_count()):
            out = []
            for col_i in range(n_cols):
                if cell := row.get_cell(col_i):
                    out.append(cell.formatted_line(width, line_i))
                else:
                    out.append(" " * width)
            self.buffer.write("│ " + " │ ".join(out) + " │\n")

    def column_count(self):
        col_count = max(x.column_count() for x in self.rows)
        if hdr := self.header:
            hdr_cols = hdr.column_count()
            if hdr_cols > col_count:
                col_count = hdr_cols
        return col_count

    def max_line_length(self):
        max_ll = 0
        if hdr := self.header:
            max_ll = hdr.max_line_length()
        for r in self.rows:
            n = r.max_line_length()
            if n > max_ll:
                max_ll = n
        return max_ll
