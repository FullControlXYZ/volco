from abc import ABC, abstractmethod, abstractproperty


class Instruction(ABC):
    @abstractmethod
    def read(self):
        pass

    @abstractproperty
    def movements(self):
        pass

    @abstractproperty
    def coordinate_limits(self):
        pass

    @abstractproperty
    def number_printed_filaments(self):
        pass

    @abstractproperty
    def filaments_coordinates(self):
        pass

    @abstractproperty
    def default_nozzle_speed(self):
        pass
