import os
import sys
import logging

from nanodecon.decontaminate_nanopore import nanopore_decontamination
from nanodecon.decontaminate_illumina import illumina_decontamination

def decontaminate(arguments):
    set_up_output_and_check_input(arguments)
    if arguments.illumina != None:
        illumina_decontamination(arguments)
    if arguments.nanopore != None:
        nanopore_decontamination(arguments)
    #Plasmid: https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid

    #https://datadryad.org/stash/dataset/doi:10.15146/R33X2J pøasmid d
def set_up_output_and_check_input(arguments):
    if not os.path.exists(arguments.output):
        os.makedirs(arguments.output)

    if arguments.illumina != None:
        if not os.path.isfile(arguments.illumina[0]):
            print ('Input file does not exist')
            sys.exit()

    if arguments.nanopore != None:
        if not os.path.isfile(arguments.nanopore):
            print ('Input file does not exist')
            sys.exit()