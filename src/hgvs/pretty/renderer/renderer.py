from abc import ABC, abstractmethod


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
