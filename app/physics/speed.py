from abc import ABC, abstractmethod, abstractproperty


class Speed(ABC):
    @abstractmethod
    def calculate_displacements(self):
        pass

    @abstractproperty
    def travel_length(self):
        pass

    @abstractproperty
    def target_speed(self):
        pass

    @abstractproperty
    def threshold_speed(self):
        pass

    @abstractproperty
    def acceleration(self):
        pass

    @abstractproperty
    def total_time(self):
        pass

    @abstractproperty
    def speed_profile(self):
        pass
