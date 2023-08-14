#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'


import logging

from GPMsDB_tk.calc import CalcHit, CalcScore, CalcScore2, CalcRamdom
from GPMsDB_tk.defaultValues import DefaultValues


class SearchBestHit2(object):
    def __init__(self):
        self.logger = logging.getLogger('GPMsDB_tk')

    def run(self, input_file, reference, peaks, ppm, first, top, score_type,
            adjust, reps_db, all_db, tax, ncbi, com, genes, t_peak, t_use, minimum, no_gen):
        # first search
        self.logger.info('[identify] 1st search.')

        hit, exact = CalcHit(self, peaks, adjust, ppm, 1, reps_db, score_type)

        # second search
        self.logger.info('[identify] 2nd search.')

        hit_all, exact_all = CalcHit(
            self, peaks, adjust, ppm, 1, all_db, score_type)

        # Calculating score
        self.logger.info('[identify] Calculating score.')
        scores = {}
        if no_gen == "limit":
            scores = CalcScore(self, score_type, hit,
                               hit_all, exact, exact_all, genes)
        elif no_gen == "linear":
            scores = CalcScore2(self, score_type, hit,
                                hit_all, exact, exact_all, genes)

        # output
        self.logger.info("Searching done.")

        return scores
