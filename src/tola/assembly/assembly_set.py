"""
Stores a set of assemblies for an organism.
"""

from collections.abc import ItemsView, KeysView, MutableMapping, ValuesView

from tola.assembly.assembly import Assembly
from tola.assembly.naming_utils import natural_key


class AssemblySetError(Exception):
    """
    Error from `AssemblySet` use.
    """


class AssemblySet(MutableMapping):
    def __init__(self, asm_dict: dict[str | None, Assembly] | None = None):
        self.__asm_dict: dict[str | None, Assembly] = {**asm_dict} if asm_dict else {}

    def __getitem__(self, name: str | None) -> Assembly:
        return self.__asm_dict[name]

    def __setitem__(self, name: str | None, asm: Assembly):
        self.__asm_dict[name] = asm

    def __delitem__(self, name: str | None):
        del self.__asm_dict[name]

    def __len__(self):
        return len(self.__asm_dict)

    def __iter__(self):
        yield from self.__asm_dict

    def items(self) -> ItemsView[str | None, Assembly]:
        return self.__asm_dict.items()

    def keys(self) -> KeysView[str | None]:
        return self.__asm_dict.keys()

    def values(self) -> ValuesView[Assembly]:
        return self.__asm_dict.values()

    def curated(self) -> dict[str | None, Assembly]:
        """
        Returns a list of all the curated assemblies
        """

        return {n: asm for n, asm in self.__asm_dict.items() if asm.curated}

    def main_assembly(self) -> tuple[str | None, Assembly]:
        """
        Returns the primary curated assembly in the set, which is either
        stored under the keys "Primary" or `None`, or is the first given by
        the list of sorted, cuated assemblies, sorted naturally, *i.e.*
        `hap1` from `["hap1", "hap2", ...]`.
        """

        curated = self.curated()
        if not curated:
            msg = (
                "Cannot choose a main Assembly from AssemblySet"
                " containing zero curated assemblies"
            )
            raise AssertionError(msg)

        for asm_key in "Primary", None:
            if asm := curated.get(asm_key):
                return asm_key, asm

        by_name = sorted(curated, key=natural_key)  # ty: ignore[no-matching-overload]
        main = by_name[0]
        return main, curated[main]
