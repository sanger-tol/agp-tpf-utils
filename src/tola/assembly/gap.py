class Gap:
    __slots__ = "_length", "_gap_type"

    def __init__(self, length, gap_type, tags=()):
        self._length = int(length)
        self._gap_type = str(gap_type)

    @property
    def length(self):
        return self._length

    @property
    def gap_type(self):
        return self._gap_type

    def __repr__(self):
        return f"{self.__class__.__name__}(length={self.length}, gap_type='{self.gap_type}')"

    def __str__(self):
        return f"Gap:{self.length} {self.gap_type}"
