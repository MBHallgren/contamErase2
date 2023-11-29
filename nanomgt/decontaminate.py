import os
import sys
import logging

from nanomgt.decontaminate_nanopore import nanopore_decontamination


def decontaminate(arguments):
    """
    Perform decontamination based on the provided arguments.

    Args:
        arguments: Parsed command-line arguments.
    """
    set_up_output_and_check_input(arguments)

    if arguments.nanopore is not None:
        nanopore_decontamination(arguments)


def set_up_output_and_check_input(arguments):
    """
    Set up the output directory and check the existence of input files.

    Args:
        arguments: Parsed command-line arguments.
    """
    # Create the output directory if it doesn't exist
    if not os.path.exists(arguments.output):
        os.makedirs(arguments.output)

    # Check if a nanopore input file is provided and if it exists
    if arguments.nanopore is not None:
        if not os.path.isfile(arguments.nanopore):
            print('Input file does not exist')
            sys.exit()