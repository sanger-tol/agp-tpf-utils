import io

import click


def bold_red(txt):
    return click.style("".join(txt), bold=True, fg="red")


def bold(txt):
    return click.style("".join(txt), bold=True)


class CellLine:
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

    def __init__(self, lines: list[CellLine] = None):
        self.lines = lines if lines else []

    def add_line(self, line: CellLine):
        self.lines.append(line)

    def new_line(self, *args, **kwargs):
        line = CellLine(*args, **kwargs)
        self.lines.append(line)
        return line

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
        if cells and not isinstance(cells, list):
            lst = list(cells)
            msg = f"Not a list: {cells = !r} {lst = !r}"
            raise ValueError(msg)
        self.cells = cells if cells else []

    def add_cell(self, cell: TableCell):
        self.cells.append(cell)

    def new_cell(self, *args, **kwargs):
        cell = TableCell(*args, **kwargs)
        self.cells.append(cell)
        return cell

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


class TerminalTableError(Exception):
    """Unexpected usage of Table and associated Classes"""


class TerminalTable:
    """A class for rendering a table in the terminal"""

    __slots__ = "header", "rows", "buffer", "errors"

    def __init__(self, header: TableHeader = None, rows: list[TableRow] = None):
        self.header = header
        self.rows = rows if rows else []
        self.buffer = io.StringIO()
        self.errors = {}

    def add_row(self, row: TableRow):
        self.rows.append(row)

    def new_header(self, *args, **kwargs):
        if self.header:
            msg = "TerminalTable header is already set"
            raise TerminalTableError(msg)
        hdr = TableHeader(*args, **kwargs)
        self.header = hdr
        return hdr

    def new_row(self, *args, **kwargs):
        row = TableRow(*args, **kwargs)
        self.rows.append(row)
        return row

    def current_row_index(self):
        if n := len(self.rows):
            return n - 1
        return None

    def mark_error(self):
        i = self.current_row_index()
        err = self.errors
        err[i] = err.get(i, 0) + 1

    def error_render(self, context=1):
        if not self.errors:
            return None

        for start, end in self.contiguous_ranges(
            self.errors.keys(),
            len(self.rows),
            context,
        ):
            yield self.render(range(start, end + 1))

    @staticmethod
    def contiguous_ranges(indices: list[int], length: int, context=1):
        """
        Given a list of indices within a list of the supplied `length`,
        returns a list of contiguous ranges padded by `context`.
        """
        idxs = sorted(indices)
        max_i = length - 1
        while idxs:
            x = idxs.pop(0)
            start = max(x - context, 0)
            end = min(x + context, max_i)
            while idxs:
                y = idxs[0]
                y_start = max(y - context, 0)
                y_end = min(y + context, max_i)
                if y_start > end + 1:
                    # There's a gap
                    break
                idxs.pop(0)
                end = y_end
            yield start, end

    def render(self, row_indices=None):
        all_rows = self.rows
        if not row_indices:
            row_indices = range(len(all_rows))

        # Empty the buffer
        self.buffer.seek(0)
        self.buffer.truncate(0)

        n_cols = self.column_count()
        width = self.max_line_length()
        pad_width = width + 2

        self.buffer.write("┌" + "┬".join(n_cols * ["─" * pad_width]) + "┐\n")

        if hdr := self.header:
            self.render_row(hdr, n_cols, width)

        ruler = "├" + "┼".join(n_cols * ["─" * pad_width]) + "┤\n"

        for i in row_indices:
            row = all_rows[i]
            if hdr:
                self.buffer.write(ruler)
            else:
                hdr = True
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
