import re

from hgvs.enums import ValidationLevel


g_ac  = "|".join([r"^(?:NC_|NG_|NT_|NW_|NZ_|CM)\d+\.\d+$", r"^LRG_\d+$"])
cr_ac = "|".join([r"^(?:NM_|XM_|ENST)\d+\.\d+$",           r"^LRG_\d+t\d+$"])
n_ac  = "|".join([r"^(?:NR_|XR_|ENST)\d+\.\d+$",           r"^LRG_\d+t\d+$"])
p_ac  = "|".join([r"^(?:NP_|XP_|ENSP)\d+\.\d+$",           r"^LRG_\d+p\d+$"])
m_ac  = "|".join([r"^(?:NC_)\d+\.\d+$"])

valid_pairs = {
    "c": re.compile(cr_ac),
    "g": re.compile(g_ac),
    "m": re.compile(m_ac),
    "n": re.compile(n_ac),
    "p": re.compile(p_ac),
    "r": re.compile(cr_ac),
}

invalid_pairs = {
    "c": re.compile("|".join((g_ac, p_ac))),
    "g": re.compile("|".join((cr_ac, n_ac, p_ac))),
    "m": re.compile("|".join((cr_ac, n_ac, p_ac))),
    "n": re.compile("|".join((g_ac, p_ac))),
    "p": re.compile("|".join((g_ac, cr_ac, n_ac))),
    "r": re.compile("|".join((g_ac, p_ac))),
}


def validate_type_ac_pair(type, ac):
    """validate that accession is correct for variant type AND that
    accession is fully specified.

    """

    assert type in valid_pairs, "Unknown variant type " + type
    if valid_pairs[type].match(ac):
        return (ValidationLevel.VALID,
                "Accession ({ac}) is compatible with variant type {type}".format(ac=ac, type=type))
    elif invalid_pairs[type].match(ac):
        return (ValidationLevel.ERROR,
                "Accession ({ac}) is not compatible with variant type {type}".format(ac=ac, type=type))
    else:
        return (ValidationLevel.WARNING,
                "Accession ({ac}) is not known to be compatible with variant type {type}".format(ac=ac, type=type))
