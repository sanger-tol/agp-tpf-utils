class Gap:
    __slots__ = "_length", "_gap_type", "_meta"

    def __init__(self, length, gap_type, meta=None):
        self._length = int(length)
        self._gap_type = str(gap_type)
        self._meta = meta

    @property
    def length(self):
        return self._length

    @property
    def gap_type(self):
        return self._gap_type

    @property
    def meta(self):
        return self._meta

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(length={self.length}, gap_type='{self.gap_type}'"
            + (f", meta='{self.meta}')" if self.meta else ")")
        )

    def __str__(self):
        return f"Gap:{self.length} {self.gap_type}" + (
            f" {self.meta}" if self.meta else ""
        )
