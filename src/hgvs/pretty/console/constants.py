### Color constants for rendering text on the console

TGREEN = "\033[32m"  # Green Text
TGREENBG = "\033[30;42m"
TRED = "\033[31m"  # Red Text
TREDBG = "\033[30;41m"
TBLUE = "\033[34m"  # Blue Text
TBLUEBG = "\033[30;44m"
TPURPLE = "\033[35m"  # Purple Text
TPURPLEBG = "\033[30;45m"
TYELLOW = "\033[33m"  # Yellow Text
TYELLOWBG = "\033[30;43m"

ENDC = "\033[m"  # reset to the defaults

# standard colors for nucleotides and special codons
COLOR_MAP = {
    "A": TGREEN,
    "T": TRED,
    "C": TYELLOW,
    "G": TBLUE,
    "N": TPURPLE,
    "init_met": TGREEN,
    "stop_codon": TRED,
    "codon1": TBLUE,
    "codon2": TYELLOW,
    "tx_ref_disagree": TRED,
    "del": TRED,
}
