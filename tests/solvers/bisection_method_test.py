import pytest

from app.solvers.bisection_method import BisectionMethod


def second_degree_fun(x, arg_a, arg_b):
    return x**2 + x * arg_a + arg_b, x


class TestBisectionMethod:
    def test_should_find_roots_of_a_second_degree_function(self):
        # fun = x^2 + 2x - 1

        args = 2, -1
        fn, x_result = BisectionMethod().execute(
            fun=second_degree_fun,
            initial_point=0,
            tolerance=1e-6,
            increment=0.1,
            args=args,
        )

        assert abs(fn) < 1e-6

        expected_x = -1.0 + 2**0.5
        assert abs(x_result - expected_x) < 1e-6
