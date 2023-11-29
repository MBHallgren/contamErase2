import os
import sys
import logging
import subprocess
import multiprocessing


class KmergenetyperRunner():
    def __init__(self, illumina, nanopore, fasta, arguments):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.check_for_kmergenetyper()
        self.illumina = illumina
        self.nanopore = nanopore
        self.fasta = fasta
        self.arguments = arguments

    def check_for_kmergenetyper(self):
        """Checks if kmergenetyper is installed"""
        #Requires kmergenetyper to be installed in the PATH
        try:
            subprocess.call(["kmergenetyper"], stdout=open(os.devnull, 'wb'))
        except FileNotFoundError:
            self.logger.info("kmergenetyper is not installed correctly directly in the PATH.")
            sys.exit(1)

    def run(self):
        """runs kmergenetyper"""
        cores = multiprocessing.cpu_count()
        kmergenetyper_cmd = 'kmergenetyper -keep'
        if self.illumina != None:
            kmergenetyper_cmd += ' -illumina {}'.format(self.illumina)
        if self.nanopore != None:
            kmergenetyper_cmd += ' -nanopore {}'.format(self.nanopore)
        if self.fasta != None:
            kmergenetyper_cmd += ' -fasta {}'.format(self.fasta)
        if self.arguments != None:
            kmergenetyper_cmd += ' {}'.format(self.arguments)
        print (kmergenetyper_cmd)
        self.logger.info("Running kmergenetyper with the following command: {}".format(kmergenetyper_cmd))
        os.system(kmergenetyper_cmd)