from hgvs.utils.orderedenum import OrderedEnum

Datum = OrderedEnum("Datum", "SEQ_START CDS_START CDS_END")

ValidationLevel = OrderedEnum("ValidationLevel", "VALID WARNING ERROR")

PrevalidationLevel = OrderedEnum("PrevalidationLevel", "NONE INTRINSIC EXTRINSIC")
