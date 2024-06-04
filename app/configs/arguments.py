import argparse


class Arguments:
    @staticmethod
    def get_options():
        parser = argparse.ArgumentParser()
        parser.add_argument("--gcode", type=str)
        parser.add_argument("--sim", type=str)
        parser.add_argument("--printer", type=str)

        return parser.parse_args()
