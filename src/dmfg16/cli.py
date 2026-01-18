import argparse
from ._workflow import main

class CLICommand:
    """Run DMF-based double-ended TS search with Gaussian.

    This command converts a Gaussian QST input into a TS input with
    a single TS guess structure using DMF and performs the TS optimization.
    """

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument

        add(
            "-n", "--npoints",
            type=int,
            default=3,
            help=(
                "Number of energy evaluation points along the reaction path "
                "(default: 3)"
            )
        )

        add(
            "-e", "--equidistant",
            action="store_true",
            help=(
                "Keep the distnace between energy evaluation points equal"
            )
        )

        add(
            "-t", "--tol",
            default="middle",
            help=(
                'Condition for stop DMF optimization: "loose", "middle", or "tight"'
                '(default: "middle")'
            )
        )

        add(
            "-d", "--dir",
            default="dmf",
            help='DMF working directory (default: "dmf")'
        )

        add(
            "--exe",
            default="g16",
            help='Gaussian executable to use (default: "g16")'
        )

        add(
            "-P", "--parallel",
            action="store_true",
            help="Enable DMF-level parallel execution"
        )

        add(
            "-C", "--only-ts-guess",
            action="store_true",
            help=(
                "Stop after TS-guess generation phase "
                "(skip subsequent opt=ts calculation)"
            )
        )

        add(
            "-R", "--restart",
            action="store_true",
            help=(
                "Restart DMF optimization  "
            )
        )

    @staticmethod
    def run(args=None):
        main(args)


def g16(prog="dmf-g16", args=None):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="DMF Gaussian wrapper"
    )

    CLICommand.add_arguments(parser)

    ns = parser.parse_args(args)
    CLICommand.run(ns)
