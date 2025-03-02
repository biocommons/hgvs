from abc import ABC, abstractmethod

from hgvs.pretty.console.constants import ENDC, TPURPLE, TYELLOW
from hgvs.pretty.models import PrettyConfig


def colorize_hgvs(hgvs_str: str) -> str:
    """Takes a string representation of a hgvs Sequence Variant and renders it with console colors."""

    spl = hgvs_str.split(":")
    var_str = TPURPLE + spl[0] + ENDC
    var_str += ":"

    sec = spl[1].split(".")
    var_str += TYELLOW + sec[0] + ENDC
    var_str += "."
    var_str += sec[1]

    return var_str


class BasicRenderer(ABC):
    """
    BasicRenderer is an abstract base class that provides a template for rendering objects with a specific configuration and orientation.
    Attributes:
        config: Configuration settings for the renderer.
        orientation (int): Orientation setting for the renderer.
    Methods:
        legend():
            Abstract method to generate a legend for the rendered object.
        display():
            Abstract method to display the rendered object.
    """

    def __init__(self, config: PrettyConfig, orientation: int):
        self.config = config
        self.orientation = orientation

    @abstractmethod
    def legend(self):
        pass

    @abstractmethod
    def display(self):
        pass
