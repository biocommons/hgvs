#
# Base class for inputs
#
# TODO - work out what the actual inputs would be

import abc

class BaseInputSource(object):

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_sequence(self, accession):
        """Get sequence
        :param tag: accession accession number
        :return sequence as string
        """
        # TODO - decide if this is better just labeled as generic "id" as opposed to accession
        return None

    # TODO - actually define this one in some way
    @abc.abstractmethod
    def get_chromosome(self):
        return None




def main():
    pass


if __name__ == "__main__":
    main()

