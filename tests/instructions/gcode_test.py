import pytest

from app.instructions.gcode import Gcode


class TestGcode:
    def test_should_read_the_gcode(self):
        gcode = Gcode(
            gcode_path="tests/fixtures/gcode_example.gcode", default_nozzle_speed=40.0
        )

        gcode.read()

        assert gcode.number_printed_filaments == 3

        assert len(gcode.movements) == 8

        assert len(gcode.filaments_coordinates) == 3

        assert gcode.coordinate_limits["x"] == [10.0, 14.0]
        assert gcode.coordinate_limits["y"] == [8.0, 12.0]
        assert gcode.coordinate_limits["z"] == [0.0, 0.7]
