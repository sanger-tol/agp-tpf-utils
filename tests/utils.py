"""
For functions shared by test files.
"""

from textwrap import dedent


def strip_leading_spaces(txt):
    """
    Removes leading blank lines and de-indents text, so that the test data in
    this file can be indented to the code.
    """
    return dedent(txt).lstrip()
