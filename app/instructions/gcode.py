import logging
import math
import re

from app.instructions.instruction import Instruction

logger = logging.getLogger(__name__)


class Gcode(Instruction):
    def __init__(self, gcode_path, default_nozzle_speed, printer):
        self.gcode_path = gcode_path
        self._movements = list()
        self._coordinate_limits = {}
        self._number_printed_filaments = 0
        self._filaments_coordinates = list()
        self._default_nozzle_speed = default_nozzle_speed
        self._printer = printer  # Needed for E to volume conversion

    @property
    def movements(self):
        return self._movements

    @property
    def coordinate_limits(self):
        return self._coordinate_limits

    @property
    def number_printed_filaments(self):
        return self._number_printed_filaments

    @property
    def filaments_coordinates(self):
        return self._filaments_coordinates

    @property
    def default_nozzle_speed(self):
        return self._default_nozzle_speed

    def read(self):
        gcode_file = open(self.gcode_path, "r")

        # Default values
        flag_relative = 0  # G90 -> absolute printing
        e_relative = 0  # M82 -> absolute extrusion
        unit_mode = 'mm'  # default units ('mm' or 'inches')

        # Coord list = [Xabs, Yabs, Zabs, Erelative, Vprint]
        movements = list()
        movements.append([0.0, 0.0, 0.0, 0.0, self.default_nozzle_speed])

        vprint = self.default_nozzle_speed

        extrusion_old = 0.0

        logger.info("Processing .gcode ...")

        for line_raw in gcode_file:
            # Strip UTF-8 BOM if present
            line_raw = line_raw.lstrip('\ufeff')
            # 1. Strip comments
            line = line_raw.split(";", 1)[0]
            # 2. Convert to uppercase
            line = line.upper()
            # 3. Strip leading/trailing whitespace
            line = line.strip()

            if not line:  # Skip empty lines or lines that were only comments
                continue

            parts = line.split()

            # 4. Handle line numbers (N codes)
            if parts[0].startswith('N') and len(parts[0]) > 1 and parts[0][1:].isdigit():
                parts.pop(0)
                if not parts:  # Line might have only contained N number
                    continue

            command = parts[0]
            params = {}
            for part in parts[1:]:
                # 5. Missing numeric value
                if len(part) == 1 and part[0] in "XYZEF":
                    logger.warning(
                        f"Missing numeric value for parameter '{part}' from line: {line_raw.strip()}")
                    raise ValueError("Please correct gcode format and retry")
                if len(part) > 1 and part[0] in "GMXYZEF":
                    val_str = part[1:]
                    # 6. Malformed number
                    if not re.match(r'^[-+]?(?:\d+\.?\d*|\.\d+)$', val_str):
                        logger.warning(
                            f"Malformed numeric value '{val_str}' in param '{part}' from line: {line_raw.strip()}")
                        raise ValueError(
                            "Please correct gcode format and retry")
                    try:
                        params[part[0]] = float(val_str)
                    except ValueError:
                        logger.warning(
                            f"Could not parse parameter value in '{part}' from line: {line_raw.strip()}")
                        raise ValueError(
                            "Please correct gcode format and retry")
                # Silently ignore other parts/parameters we don't understand

            # Convert coordinates from inches to mm if needed
            if unit_mode == 'inches':
                for axis in ("X", "Y", "Z"):
                    if axis in params:
                        params[axis] *= 25.4

            # Update feedrate if F parameter exists
            if "F" in params:
                vprint = params["F"] / 60.0  # transforming to mm/s

            # Update positioning modes
            if command == "G90":
                flag_relative = 0  # G90 -> absolute positioning
            elif command == "G91":
                flag_relative = 1  # G91 -> relative positioning
            # Units: G20 (inches), G21 (mm)
            elif command == "G20":
                unit_mode = 'inches'
            elif command == "G21":
                unit_mode = 'mm'

            # Update extrusion modes
            if command == "M82":
                e_relative = 0  # M82 -> absolute extrusion
            elif command == "M83":
                e_relative = 1  # M83 -> relative extrusion
            # Warn on unsupported M-codes
            elif command.startswith("M"):
                logger.warning(
                    f"Unsupported M-code '{command}' on line: {line_raw.strip()}")

            # Handle position reset (critical for absolute extrusion)
            if command == "G92":
                # Handle G92 resets
                if "E" in params:
                    if params["E"] == 0.0:
                        extrusion_old = 0.0
                        logger.debug("Resetting extrusion reference (G92 E0)")
                    else:
                        logger.warning(
                            f"G92 sets E to non-zero ({params['E']}) on line: {line_raw.strip()} – unexpected extrusion reset")
                        raise ValueError(
                            "Please correct gcode format and retry")
                if any(axis in params for axis in ("X", "Y", "Z")):
                    axes = [axis for axis in params if axis in "XYZ"]
                    logger.warning(
                        f"G92 resets position axes {axes} on line: {line_raw.strip()} – unsupported coordinate reset")
                    raise ValueError("Please correct gcode format and retry")

            # Handle movement commands
            elif command == "G1" or command == "G0":
                # Check if there's actual movement data
                if any(axis in params for axis in "XYZE"):
                    coord_new, extrusion_old = self._define_movement(
                        params,  # Pass parsed parameters
                        movements[-1],
                        flag_relative,
                        e_relative,
                        vprint,
                        extrusion_old,
                    )
                    movements.append(coord_new)

        gcode_file.close()

        # Ensure coordinate limits calculation happens *after* the loop
        xlim, ylim, zlim, nfil, coord_fil = self._max_min_extru_coordinates(
            movements)

        logger.info("Done processing .gcode!")

        self._number_printed_filaments = nfil
        self._movements = movements
        self._coordinate_limits = {"x": xlim, "y": ylim, "z": zlim}
        self._filaments_coordinates = coord_fil

    def _define_movement(
        self, params, coord_old, flag_relative, e_relative, vprint, extrusion_old
    ):
        # Initialize new coordinates with old ones
        xnew, ynew, znew = coord_old[0], coord_old[1], coord_old[2]
        enew = 0.0  # Default to no extrusion for this move

        # Update coordinates based on parameters and relative/absolute mode
        if "X" in params:
            xnew = params["X"] + flag_relative * coord_old[0]
        if "Y" in params:
            ynew = params["Y"] + flag_relative * coord_old[1]
        if "Z" in params:
            znew = params["Z"] + flag_relative * coord_old[2]

        # Handle extrusion
        if "E" in params:
            enow = params["E"]
            if e_relative == 0:  # Absolute extrusion
                enew = enow - extrusion_old  # Calculate relative extrusion for this move
                extrusion_old = enow  # Update the absolute reference
            else:  # Relative extrusion
                enew = enow
                # extrusion_old doesn't need updating in relative mode, but G92 E0 resets it

        # Store the calculated movement details
        # Note: enew here represents the *relative* extrusion for this specific segment
        coord_new = [xnew, ynew, znew, enew, vprint]

        return coord_new, extrusion_old  # Return updated absolute reference for E

    """
    This function returns the minimum and maximum printing coordinates;
    """

    def _max_min_extru_coordinates(self, coord_list):
        extru_now = coord_list[0][3]

        xlist = list()
        ylist = list()
        zlist = list()

        nfil = 0  # total number of filaments

        coord_fil = (
            list()
        )  # list of initial and final coordinates of each printed filament
        # coord_fil[0] = [coord_i, coord_f] -> of filament 1

        for index in range(1, len(coord_list)):

            coord_now = coord_list[index]

            # sum the relative extrusion coordinates.
            extru_now += coord_now[3]
            # this is done in order to account for retractiong movements.

            zlist.append(coord_now[2])

            if extru_now > 0:
                nfil += 1

                # initial coordinates of the filament
                coord_old = coord_list[index - 1]

                # Convert E to volume
                volume = extru_now * math.pi * \
                    (self._printer.feedstock_filament_diameter/2)**2

                coord_fil.append([coord_old, coord_now, volume])

                extru_now = 0.0  # reference is set to zero after printing a filament

                xlist.append(coord_old[0])
                xlist.append(coord_now[0])

                ylist.append(coord_old[1])
                ylist.append(coord_now[1])

                zlist.append(coord_old[2])
                zlist.append(coord_now[2])

        xlim = [min(xlist), max(xlist)]
        ylim = [min(ylist), max(ylist)]
        zlim = [0.0, max(zlist)]

        return xlim, ylim, zlim, nfil, coord_fil
