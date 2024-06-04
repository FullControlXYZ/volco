import logging


logger = logging.getLogger(__name__)


def default_increase_solver_tolerance(point_a, point_b):
    return False


class BisectionMethod:
    def execute(
        self,
        fun,
        initial_point,
        tolerance,
        increment,
        args,
        fun_increase_tolerance=default_increase_solver_tolerance,
    ):
        if initial_point == 0.0:
            point_b = increment
        else:
            point_b = initial_point

        _, _, point_a, point_b = self._loop_fb(fun, point_b, increment, args)

        return self._loop_fc(
            fun, point_a, point_b, tolerance, fun_increase_tolerance, args
        )

    def _loop_fb(self, fun, point_b, inc, args):
        fb = -1.0

        point_a = 0.0

        while fb < 0.0:

            fb, out = fun(point_b, *args)

            if fb < 0:
                point_a = point_b
                point_b += inc

        return fb, out, point_a, point_b

    def _loop_fc(self, fun, point_a, point_b, tolerance, fun_increase_tolerance, args):
        fc = 2.0 * tolerance

        while abs(fc) > tolerance:
            point_c = (point_a + point_b) * 0.5

            fc, out = fun(point_c, *args)

            if fun_increase_tolerance(point_a, point_b):
                logger.debug(
                    "[BisectionMethod]: increasing tolerance because point_b and point_a are too close"
                )
                tolerance = tolerance * 10

            if fc < 0.0:
                point_a = point_c
            else:
                point_b = point_c

        return fc, out
