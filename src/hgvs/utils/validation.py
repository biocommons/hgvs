import re

from hgvs.enums import ValidationLevel

g_ac = r"|".join([r"^(?:NC_|NG_|NT_|NW_|NZ_|CM)\d+\.\d+$", r"^LRG_\d+$"])
cr_ac = r"|".join([r"^(?:NM_|XM_|ENST)\d+\.\d+$", r"^LRG_\d+t\d+$"])
n_ac = r"|".join([r"^(?:NR_|XR_|ENST)\d+\.\d+$", r"^LRG_\d+t\d+$"])
p_ac = r"|".join([r"^(?:NP_|XP_|ENSP)\d+\.\d+$", r"^LRG_\d+p\d+$"])
m_ac = r"|".join([r"^(?:NC_)\d+\.\d+$"])

valid_pairs = {
    "c": re.compile(cr_ac),
    "g": re.compile(g_ac),
    "m": re.compile(m_ac),
    "n": re.compile(n_ac),
    "p": re.compile(p_ac),
    "r": re.compile(cr_ac),
}

invalid_pairs = {
    "c": re.compile(r"|".join((g_ac, p_ac))),
    "g": re.compile(r"|".join((cr_ac, n_ac, p_ac))),
    "m": re.compile(r"|".join((cr_ac, n_ac, p_ac))),
    "n": re.compile(r"|".join((g_ac, p_ac))),
    "p": re.compile(r"|".join((g_ac, cr_ac, n_ac))),
    "r": re.compile(r"|".join((g_ac, p_ac))),
}


def validate_type_ac_pair(variant_type, ac) -> tuple[ValidationLevel, str]:
    """validate that accession is correct for variant type AND that
    accession is fully specified.

    """
    if variant_type not in valid_pairs:
        raise ValueError(f"Unknown variant type {variant_type}")
    if valid_pairs[variant_type].match(ac):
        return (
            ValidationLevel.VALID,
            f"Accession ({ac}) is compatible with variant type {variant_type}",
        )
    if invalid_pairs[variant_type].match(ac):
        return (
            ValidationLevel.ERROR,
            f"Accession ({ac}) is not compatible with variant type {variant_type}",
        )
    return (
        ValidationLevel.WARNING,
        f"Accession ({ac}) is not known to be compatible with variant type {variant_type}",
    )
