import os
import unittest

import pytest
from support import CACHE

import hgvs.assemblymapper
import hgvs.dataproviders.uta


@pytest.mark.issues
class Test_Issues(unittest.TestCase):
    """
    HGVS-730 fixes an issue where AARefAlt.format() sometimes returned the three letter amino acid
    when the single letter format was requested.
    """

    @classmethod
    def setUpClass(self):
        self.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
        )
        self.am37 = hgvs.assemblymapper.AssemblyMapper(
            self.hdp, replace_reference=True, assembly_name="GRCh37", alt_aln_method="splign"
        )
        self.hp = hgvs.parser.Parser()

    def test_730_format_startloss_as_configured(self):
        """
        Parse a start loss and make sure the one letter and three letter protein genotypes come back as expected.
        """
        var_g = self.hp.parse_hgvs_variant("NC_000016.9:g.89985662_89985667del")
        var_c = self.am37.g_to_c(var_g, str("NM_002386.3"))
        var_p = self.am37.c_to_p(var_c)
        var_p_one_letter = var_p.format(conf={"p_3_letter": False})
        var_p_three_letter = var_p.format(conf={"p_3_letter": True})

        self.assertEqual(var_p_one_letter, "NP_002377.4:p.M1?")
        self.assertEqual(var_p_three_letter, "NP_002377.4:p.Met1?")
