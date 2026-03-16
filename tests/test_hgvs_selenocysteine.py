"""Tests for correct HGVS translations of selenocysteine variants.

Selenocysteine (Sec) is encoded by the UGA codon, which is normally a stop codon.
This test verifies that variants affecting selenocysteine positions are correctly
translated to protein variants.
"""

import os
import unittest

import pytest

import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
from support import CACHE


@pytest.mark.quick
class TestSelenocysteineTranslations(unittest.TestCase):
    """Test correct HGVS translations for selenocysteine variants."""

    @classmethod
    def setUpClass(cls):
        cls.hdp = hgvs.dataproviders.uta.connect(
            mode=os.environ.get("HGVS_CACHE_MODE", "run"), cache=CACHE
        )
        cls.vm = hgvs.variantmapper.VariantMapper(cls.hdp)
        cls.hp = hgvs.parser.Parser()
        cls.am38 = hgvs.assemblymapper.AssemblyMapper(
            cls.hdp, assembly_name="GRCh38", alt_aln_method="splign"
        )

    def test_selenocysteine_coding_deletion(self):
        """Test c_to_p conversion for deletion affecting selenocysteine codon.

        NM_020451.3:c.380_382del deletes the UGA codon that encodes selenocysteine.
        This should result in a protein variant that correctly represents the deletion.
        """
        hgvs_c = "NM_020451.3:c.380_382del"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.vm.c_to_p(var_c)

        # Verify that the protein variant is created successfully
        assert var_p is not None
        # The variant should be a deletion at the protein level
        # The exact position and format may vary, but it should be a valid protein variant
        assert var_p.type == "p"
        assert str(var_p).startswith("NP_")

        assert str(var_p) == "NP_065184.2:p.(Sec127del)"

    def test_selenocysteine_genomic_to_protein(self):
        """Test g_to_c and c_to_p conversion for genomic variant affecting selenocysteine.

        NC_000001.10:g.26139280T>G is a genomic variant that should map to
        a coding variant and then to a protein variant affecting selenocysteine.
        """
        hgvs_g = "NC_000001.11:g.25812789T>C"
        var_g = self.hp.parse_hgvs_variant(hgvs_g)

        # Map genomic variant to coding variant for NM_020451.3
        # First, we need to find which transcript this maps to
        # For this test, we'll try to map it to NM_020451.3 if it overlaps

        var_c = self.am38.g_to_c(var_g, "NM_020451.3")
        var_p = self.am38.c_to_p(var_c)

        # Verify that the protein variant is created successfully
        assert var_p is not None
        assert var_p.type == "p"
        assert str(var_p).startswith("NP_")
        assert str(var_p) == "NP_065184.2:p.(Sec462Arg)"

    def test_selenocysteine_protein_variant_format(self):
        """Test that protein variants involving selenocysteine use correct format.

        Selenocysteine should be represented as 'Sec' in three-letter format
        or 'U' in single-letter format when present in protein variants.
        """
        hgvs_c = "NM_020451.3:c.380_382del"
        var_c = self.hp.parse_hgvs_variant(hgvs_c)
        var_p = self.vm.c_to_p(var_c)

        # Get the protein variant string
        var_p_str = str(var_p)

        # Check that it's a valid protein variant format
        assert ":p." in var_p_str

        # Test both one-letter and three-letter formats
        var_p_one_letter = var_p.format(conf={"p_3_letter": False})
        var_p_three_letter = var_p.format(conf={"p_3_letter": True})

        # Both should be valid protein variant strings
        assert ":p." in var_p_one_letter
        assert ":p." in var_p_three_letter

        assert str(var_p) == "NP_065184.2:p.(Sec127del)"
        assert str(var_p_one_letter) == "NP_065184.2:p.(U127del)"
        assert str(var_p_three_letter) == "NP_065184.2:p.(Sec127del)"


if __name__ == "__main__":
    unittest.main()

# <LICENSE>
# Copyright 2024 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
