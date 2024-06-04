import logging

from app.instructions.instruction import Instruction

logger = logging.getLogger(__name__)


class Gcode(Instruction):
    def __init__(self, gcode_path, default_nozzle_speed):
        self.gcode_path = gcode_path
        self._movements = list()
        self._coordinate_limits = {}
        self._number_printed_filaments = 0
        self._filaments_coordinates = list()
        self._default_nozzle_speed = default_nozzle_speed

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

        # Coord list = [Xabs, Yabs, Zabs, Erelative, Vprint]
        movements = list()
        movements.append([0.0, 0.0, 0.0, 0.0, self.default_nozzle_speed])

        vprint = self.default_nozzle_speed

        extrusion_old = 0.0

        logger.info("Processing .gcode ...")

        for line in gcode_file:

            # check if the printing speed is changed
            if line.find("F") >= 0:
                str_now = line.split("F")
                str_now = str_now[-1]

                flag_speed = -1
                cont_speed = -1
                while flag_speed == -1:
                    cont_speed += 1

                    try:
                        str_i = str_now[cont_speed]
                        if (
                            str_i == " "
                            or str_i == "\n"
                            or str_i == ";"
                            or str_i == r"\\"
                        ):
                            flag_speed = 1
                    except:
                        str_i = "\n"
                        flag_speed = 1
                        cont_speed = len(str_now)

                vprint = float(str_now[0:cont_speed]) / 60.0  # transforming to mm/s

            if line.find("G90") >= 0:
                flag_relative = 0  # G90 -> absolute printing
            elif line.find("G91") >= 0:
                flag_relative = 1  # G91 -> relative printing

            if line.find("M82") >= 0:
                e_relative = 0  # M82 -> absolute extrusion
            elif line.find("M83") >= 0:
                e_relative = 1  # M83 -> relative extrusion

            elif len(line) >= 2:
                aux_line = line[0:2]  # reads the first 2 characters

                if aux_line == "G1" or aux_line == "G0":  # it means there is movement
                    coord_new, extrusion_old = self._define_movement(
                        line,
                        movements[-1],
                        flag_relative,
                        e_relative,
                        vprint,
                        extrusion_old,
                    )

                    movements.append(coord_new)

        gcode_file.close()

        xlim, ylim, zlim, nfil, coord_fil = self._max_min_extru_coordinates(movements)

        logger.info("Done processing .gcode!")

        self._number_printed_filaments = nfil
        self._movements = movements
        self._coordinate_limits = {"x": xlim, "y": ylim, "z": zlim}
        self._filaments_coordinates = coord_fil

    def _define_movement(
        self, str_move, coord_old, flag_relative, e_relative, vprint, extrusion_old
    ):
        str_move = str_move.split(" ")

        str_move[-1] = str_move[-1].replace("\n", "")

        xnew, ynew, znew, enew = coord_old[0], coord_old[1], coord_old[2], 0.0

        for str_i in str_move:
            if str_i[0] == "X":  # nozzle moves in the x direction
                xnew = float(str_i[1::]) + flag_relative * coord_old[0]

            elif str_i[0] == "Y":  # nozzle moves in the y direction
                ynew = float(str_i[1::]) + flag_relative * coord_old[1]

            elif str_i[0] == "Z":  # nozzle moves in the z direction
                znew = float(str_i[1::]) + flag_relative * coord_old[2]

            elif str_i[0] == "E":  # transform E in relative coordinates
                enow = float(str_i[1::])

                if e_relative == 0:  # if it is absolute
                    enew = enow - extrusion_old
                    extrusion_old = enow
                else:
                    enew = enow
                    extrusion_old = 0.0

        coord_new = [xnew, ynew, znew, enew, vprint]

        return coord_new, extrusion_old

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

            extru_now += coord_now[3]  # sum the relative extrusion coordinates.
            # this is done in order to account for retractiong movements.

            zlist.append(coord_now[2])

            if extru_now > 0:
                nfil += 1

                coord_old = coord_list[index - 1]  # initial coordinates of the filament

                coord_fil.append([coord_old, coord_now, extru_now])

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
