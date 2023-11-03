""" Functions shared by test files. """

from textwrap import dedent


def strip_leading_spaces(txt):
    """
    Removes leading blank lines and de-indents text so that test data can be
    indented to the code, making it more readable.
    """
    return dedent(txt).lstrip()
