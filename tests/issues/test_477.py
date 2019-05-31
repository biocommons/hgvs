import pytest

from hgvs.exceptions import HGVSInvalidIntervalError

tests = (
    # {"c": "", "g": "", "rs": "" },

    # GPHB5, GRCh37 https://www.ncbi.nlm.nih.gov/gene/122876
    {
        "c": "NM_145171.3:c.-63A>G",
        "g": "NC_000014.8:g.63785599T>C",
        "rs": "GPHB5/GRCh37/rs1299953722",
        "ex": HGVSInvalidIntervalError
    },
    {
        "c": "NM_145171.3:c.-56G>A",
        "g": "NC_000014.8:g.63785592C>T",
        "rs": "GPHB5/GRCh37/rs982881702"
    },
    {
        "c": "NM_145171.3:c.2T>C",
        "g": "NC_000014.8:g.63784562A>G",
        "rs": "GPHB5/GRCh37/rs1221379530"
    },
    {
        "c": "NM_145171.3:c.388A>G",
        "g": "NC_000014.8:g.63779647T>C",
        "rs": "GPHB5/GRCh37/rs1380832691"
    },
    {
        "c": "NM_145171.3:c.*4C>T",
        "g": "NC_000014.8:g.63779638G>A",
        "rs": "GPHB5/GRCh37/rs753041439"
    },
    {
        "c": "NM_145171.3:c.*84A>G",
        "g": "NC_000014.8:g.63779558T>C",
        "rs": "GPHB5/GRCh37/rs1204774077"
    },
    {
        "c": "NM_145171.3:c.*99G>A",
        "g": "NC_000014.8:g.63779543C>T",
        "rs": "GPHB5/GRCh37/rs144659601",
        "ex": HGVSInvalidIntervalError
    },

    # GPHB5, GRCh37 https://www.ncbi.nlm.nih.gov/gene/122876
    {
        "c": "NM_145171.3:c.-63A>G",
        "g": "NC_000014.9:g.63318885T>C",
        "rs": "GPHB5/GRCh38/rs1299953722",
        "ex": HGVSInvalidIntervalError
    },
    {
        "c": "NM_145171.3:c.-56G>A",
        "g": "NC_000014.9:g.63318878C>T",
        "rs": "GPHB5/GRCh38/rs982881702"
    },
    {
        "c": "NM_145171.3:c.2T>C",
        "g": "NC_000014.9:g.63317848A>G",
        "rs": "GPHB5/GRCh38/rs1221379530"
    },
    {
        "c": "NM_145171.3:c.388A>G",
        "g": "NC_000014.9:g.63312933T>C",
        "rs": "GPHB5/GRCh38/rs1380832691"
    },
    {
        "c": "NM_145171.3:c.*4C>T",
        "g": "NC_000014.9:g.63312924G>A",
        "rs": "GPHB5/GRCh38/rs753041439"
    },
    {
        "c": "NM_145171.3:c.*84A>G",
        "g": "NC_000014.9:g.63312844T>C",
        "rs": "GPHB5/GRCh38/rs1204774077"
    },
    {
        "c": "NM_145171.3:c.*99G>A",
        "g": "NC_000014.9:g.63312829C>T",
        "rs": "GPHB5/GRCh38/rs144659601",
        "ex": HGVSInvalidIntervalError
    },

    # COX6A2 https://www.ncbi.nlm.nih.gov/gene/1339
    {
        "c": "NM_005205.3:c.-106G>A",
        "g": "NC_000016.10:g.31428431C>T",
        "rs": "COX6A2/GRCh38/rs1033792906",
        "ex": HGVSInvalidIntervalError
    },
    {
        "c": "NM_005205.3:c.-96C>T",
        "g": "NC_000016.10:g.31428421G>A",
        "rs": "COX6A2/GRCh38/rs755670336"
    },
    {
        "c": "NM_005205.3:c.2T>C",
        "g": "NC_000016.10:g.31428324A>G",
        "rs": "COX6A2/GRCh38/rs200780049"
    },
    {
        "c": "NM_005205.3:c.293G>A",
        "g": "NC_000016.10:g.31427775C>T",
        "rs": "COX6A2/GRCh38/rs764753905"
    },
    {
        "c": "NM_005205.3:c.*3C>T",
        "g": "NC_000016.10:g.31427771G>A",
        "rs": "COX6A2/GRCh38/rs909673485"
    },
    {
        "c": "NM_005205.3:c.*42G>C",
        "g": "NC_000016.10:g.31427732C>G",
        "rs": "COX6A2/GRCh38/rs375688325"
    },
    {
        "c": "NM_005205.3:c.*43A>G",
        "g": "NC_000016.10:g.31427731T>C",
        "rs": "COX6A2/GRCh38/rs961248971"
    },
    {
        "c": "NM_005205.3:c.*44G>A",
        "g": "NC_000016.10:g.31427730C>T",
        "rs": "COX6A2/GRCh38/rs756406653",
        "ex": HGVSInvalidIntervalError
    },
)


@pytest.mark.parametrize("pair", tests, ids=[p["rs"] for p in tests])
def test_pair(parser, am38, pair):
    var_c = parser.parse(pair["c"])
    var_g = parser.parse(pair["g"])
    if "ex" in pair:
        with pytest.raises(pair["ex"]):
            var_gtoc = am38.g_to_c(var_g, var_c.ac)
    else:
        var_gtoc = am38.g_to_c(var_g, var_c.ac)
        assert pair["c"] == str(var_gtoc)
