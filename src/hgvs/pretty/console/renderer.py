from abc import ABC, abstractmethod

from hgvs.pretty.console.constants import ENDC, TPURPLE, TYELLOW


def colorize_hgvs(hgvs_str: str) -> str:
    """ Takes a string representation of a hgvs Sequence Variant and renders it with console colors.
    """

    spl = hgvs_str.split(":")
    var_str = TPURPLE + spl[0] + ENDC
    var_str += ":"

    sec = spl[1].split(".")
    var_str += TYELLOW + sec[0] + ENDC
    var_str += "."
    var_str += sec[1]

    return var_str

class BasicRenderer(ABC):
    def __init__(self, config, orientation: int):
        self.config = config
        self.orientation = orientation

    @abstractmethod
    def legend(self):
        pass

    @abstractmethod
    def display(self):
        pass
