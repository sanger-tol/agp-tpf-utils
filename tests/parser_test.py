import io

import pytest
from tola.assembly.parser import parse_agp, parse_tpf

from .utils import strip_leading_spaces


def test_parse_agp():
    agp = strip_leading_spaces(
        """
        ##agp-version 2.1
        #
        # DESCRIPTION: Generated by PretextView Version 0.2.5
        # HiC MAP RESOLUTION: 8666.611572 bp/texel

        Scaffold_1	1	21337197	1	W	scaffold_1	1	21337197	+	Painted
        Scaffold_1	21337198	21337297	2	U	100	scaffold	yes	proximity_ligation
        Scaffold_1	21337298	21917959	3	W	scaffold_21	1	580662	+
        Scaffold_1	21917960	21918059	4	U	100	scaffold	yes	proximity_ligation
        Scaffold_1	21918060	24379376	5	W	scaffold_1	21770529	24231845	-	Painted
        Scaffold_2	1	3206646	1	W	scaffold_2	1	3206646	+	Painted
        Scaffold_2	3206647	3206746	2	U	100	scaffold	yes	proximity_ligation
        Scaffold_2	3206747	3267412	3	W	scaffold_67	1	60666	+	Painted	X
        Scaffold_2	3267413	3267512	4	U	100	scaffold	yes	proximity_ligation
        Scaffold_2	3267513	28348686	5	W	scaffold_2	3206647	28287820	?	Painted
        """,
    )

    fh = io.StringIO(agp)
    a1 = parse_agp(fh, "aaBbbCccc1")
    assert str(a1) == strip_leading_spaces(
        """
        Assembly: aaBbbCccc1
          # DESCRIPTION: Generated by PretextView Version 0.2.5
          # HiC MAP RESOLUTION: 8666.611572 bp/texel

          Scaffold_1
                21_337_197  scaffold_1:1-21337197(+) Painted
                       100  Gap:100 scaffold
                   580_662  scaffold_21:1-580662(+)
                       100  Gap:100 scaffold
                 2_461_317  scaffold_1:21770529-24231845(-) Painted

          Scaffold_2
                 3_206_646  scaffold_2:1-3206646(+) Painted
                       100  Gap:100 scaffold
                    60_666  scaffold_67:1-60666(+) Painted X
                       100  Gap:100 scaffold
                25_081_174  scaffold_2:3206647-28287820(.) Painted
        """,
    )


def test_parse_tpf():
    with pytest.raises(ValueError, match=r"Gap line before first sequence fragment"):
        parse_tpf(io.StringIO("GAP	TYPE-2	200"), "gap_first")

    with pytest.raises(ValueError, match=r"Unexpected name format"):
        parse_tpf(
            io.StringIO("?	frag	scaffold_1	PLUS"),
            "bad_fragment_name",
        )

    with pytest.raises(ValueError, match=r"Wrong field count"):
        parse_tpf(
            io.StringIO("?	scaffold_2:166926-629099"),
            "too_few_fields",
        )

    tpf = strip_leading_spaces(
        """
        ?	scaffold_1:1-93024	scaffold_1	PLUS
        GAP	TYPE-2	200
        ?	scaffold_1:93225-232397	scaffold_1	PLUS
        GAP	TYPE-2	200
        ?	scaffold_1:232598-261916	scaffold_1	PLUS
        GAP	TYPE-2	200
        ?	scaffold_1:262117-906261	scaffold_1	PLUS
        ?	scaffold_2:1-166725	scaffold_2	PLUS
        GAP	TYPE-2	200
        ?	scaffold_2:166926-629099	scaffold_2	MINUS
        GAP	TYPE-2	200
        ?	scaffold_2:629300-719848	scaffold_2	MINUS
        GAP	TYPE-2	200
        ?	scaffold_2:720049-3207246	scaffold_2	PLUS
        GAP	SHORT-ARM	200
        ?	scaffold_2:3207447-3240707	scaffold_2	PLUS
        """,
    )
    fh = io.StringIO(tpf)
    a1 = parse_tpf(fh, "aaBbbCccc1")
    assert str(a1) == strip_leading_spaces(
        """
        Assembly: aaBbbCccc1

          scaffold_1
                    93_024  scaffold_1:1-93024(+)
                       200  Gap:200 scaffold
                   139_173  scaffold_1:93225-232397(+)
                       200  Gap:200 scaffold
                    29_319  scaffold_1:232598-261916(+)
                       200  Gap:200 scaffold
                   644_145  scaffold_1:262117-906261(+)

          scaffold_2
                   166_725  scaffold_2:1-166725(+)
                       200  Gap:200 scaffold
                   462_174  scaffold_2:166926-629099(-)
                       200  Gap:200 scaffold
                    90_549  scaffold_2:629300-719848(-)
                       200  Gap:200 scaffold
                 2_487_198  scaffold_2:720049-3207246(+)
                       200  Gap:200 short_arm
                    33_261  scaffold_2:3207447-3240707(+)
        """,
    )


if __name__ == "__main__":
    test_parse_agp()
    test_parse_tpf()
