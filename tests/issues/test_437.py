"""https://github.com/biocommons/hgvs/issues/437"""


def test_437(parser, am):
    """https://github.com/biocommons/hgvs/issues/437"""

    hgvs = "NR_003051.3:n.-19_-18insACT"
    hgvs = "NM_033486.1:c.119A>G"
    am38._validator = None
    var = parser.parse(hgvs)
    var_g = am.t_to_g(var)
    print(var_g)


if __name__ == "__main__":
    from hgvs.easy import *
    test_437(parser=parser, am=am37)