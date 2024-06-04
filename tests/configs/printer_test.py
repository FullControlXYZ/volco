import pytest

from app.configs.printer import Printer


class TestPrinter:
    def test_should_load_printer_settings(self):
        printer = Printer(config_path="tests/fixtures/printer_settings.json")

        assert printer.nozzle_jerk_speed == 40.0
        assert printer.extruder_jerk_speed == 5.0
        assert printer.nozzle_acceleration == 1200.0
        assert printer.extruder_acceleration == 1200.0
        assert printer.bulk_filament_diameter == 1.75
        assert printer.nozzle_diameter == 0.4
