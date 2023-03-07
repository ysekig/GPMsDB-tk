#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2022 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'CC BY-NC-SA 4.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import os
import sys
import logging
import statistics
from GPMsDB_tk.calc import CalcHit, CalcScore, CalcRamdom
from GPMsDB_tk.defaultValues import DefaultValues
from scipy import stats


class SearchBestHit(object):
    def __init__(self):
        self.logger = logging.getLogger('GPMsDB_tk')

    def run(self, input_file, reference, peaks, ppm, first, top, score_type, adjust, reps_db, 
        all_db, tax, ncbi, strain, com, genes, taxonomy, t_peak, t_use, minimum, lineage):
        # first search
        self.logger.info('[identify] 1st search.')

        hit, exact = CalcHit(self, peaks, adjust, ppm, 1, reps_db, score_type)

        ramd = 0
        if lineage == "prokaryotes":
            num1 = DefaultValues.HIT_EXCLUDE_REP_PR
            num2 = DefaultValues.HIT_EXCLUDE_ALL_PR
            num3 = DefaultValues.HIT_EXCLUDE_ALL2_PR
            prob1 =DefaultValues.PROBABIL_RIBOSOMAL_PR
            if reference == "reps":
                prob2 = DefaultValues.PROBABIL_REPS_90_PR
                prob3 = DefaultValues.PROBABIL_REPS_99_PR
            elif reference == "all" or "custom":
                prob2 = DefaultValues.PROBABIL_ALL_90_PR
                prob3 = DefaultValues.PROBABIL_ALL_99_PR
        else:
            num1 = DefaultValues.HIT_EXCLUDE_REP_FG
            num2 = DefaultValues.HIT_EXCLUDE_ALL_FG
            num3 = DefaultValues.HIT_EXCLUDE_ALL2_FG
            prob1 =DefaultValues.PROBABIL_RIBOSOMAL_FG
            if reference == "reps":
                prob2 = DefaultValues.PROBABIL_REPS_90_FG
                prob3 = DefaultValues.PROBABIL_REPS_99_FG
            elif reference == "all" or "custom":
                prob2 = DefaultValues.PROBABIL_ALL_90_FG
                prob3 = DefaultValues.PROBABIL_ALL_99_FG

        result = sorted(hit.items(), key=lambda x: x[1], reverse=True)[
            :int(first)]
        if reference == 'reps':
            result_s = sorted(hit.items(), key=lambda x: x[1], reverse=True)[
                num1:]
            if len(result_s) == 0:
                result_s = sorted(
                    hit.items(), key=lambda x: x[1], reverse=True)
            if len(result_s) <= 100:
                print(
                    'number of candidates from the 1st screening too low: ', input_file)
                ramd = 1
        elif reference == 'all' or 'custom':
            result_s = sorted(hit.items(), key=lambda x: x[1], reverse=True)[
                num2:]
            if len(result_s) == 0:
                result_s = sorted(hit.items(), key=lambda x: x[1], reverse=True)[
                    num3:]
                if len(result_s) <= 100:
                    print(
                        'number of candidates from the 1st screening too low: ', input_file)
                    ramd = 1

        dic = dict(result)

        # second search
        self.logger.info('[identify] 2nd search.')
        all_db_limit = {}
        for h in dic.keys():
            all_db_limit[h] = []
            all_db_limit[h] = all_db[h]

        hit_all, exact_all = CalcHit(
            self, peaks, adjust, ppm, 1, all_db_limit, score_type)

        # Calculating score
        self.logger.info('[identify] Calculating score.')
        scores = {}
        scores = CalcScore(self, score_type, hit, hit_all,
                           exact, exact_all, genes)

        result2 = sorted(scores.items(), key=lambda x: x[1], reverse=True)[
            :int(top)]
        dic2 = dict(result2)

        # random sampling
        if ramd == 0:
            self.logger.info(
                '[identify] Calculating scores from ramdomly selected genomes.')

            ramdom_list = list(result_s)
            hit_all2, exact_all2 = CalcRamdom(
                self, ramdom_list, peaks, adjust, ppm, 1, all_db, score_type)
            scores2 = {}
            scores2 = CalcScore(self, score_type, hit,
                                hit_all2, exact, exact_all2, genes)
            ramdomscore = []
            for cc in scores2.keys():
                ramdomscore.append(scores2[cc])

            mean = statistics.mean(ramdomscore)
            stdev = statistics.stdev(ramdomscore)
        else:
            mean = 0
            stdev = 0

        # output
        self.logger.info("Searching done.")
        dir, filename = os.path.split(input_file)
        self.logger.info('#Search result here: with m/z adjustment of ' +
                         str("{:.1f}".format(adjust)) + " ppm")
        self.logger.info('#Input: ' + str(t_peak) + ' peaks found, ' + str(t_use) +
                         ' peaks used with relative intensity higher than ' + str(round(minimum, 5)))
        self.logger.info('#Tolerance: ' + str(ppm) +
                         ' ppm; Input file: ' + str(filename))
        self.logger.info('#Score type: ' + str(score_type) +
                         '; Reference type: ' + str(reference))
        self.logger.info('#Random sampling score: ' + str(round(mean, 2)
                                                          ) + '; standard dev: ' + str(round(stdev, 2)))
        if not com == '':
            self.logger.info('#' + str(com))
        self.logger.info(
            '#Genome Id\tprotein_hit\tribosomal_hit\tscore\tprobability\tlikelihood(%)\tncbi_name\tncbi_strain\ttaxonomy_' + str(taxonomy))
        a = 0
        for k in dic2:
            if stdev == 0:
                stdev = 0.01
            upper = stats.norm.sf(x=scores[k], loc=(mean), scale=(stdev * 3))

            if hit[k] > prob1:
                likel = "50%"
                if upper < prob2:
                    likel = "90%"
                if upper < prob3:
                    likel = "99%"
            else:
                likel = "<50%"

            try:
                show_tax = str(tax[k].rstrip())
            except KeyError:
                show_tax = "not assigned"
            try:
                ncbi_name = ncbi[k].rstrip()
            except KeyError:
                ncbi_name = "not assigned"
            try:
                strain_name = strain[k].rstrip()
            except KeyError:
                strain_name = "not assigned"

            self.logger.info(str(k) + '\t' + str(hit[k] + hit_all[k]) + '\t' + str(hit[k]) + '\t' + str(round(scores[k], 3)) + '\t' + str(
                "{:.2e}".format(upper)) + '\t' + likel + '\t' + ncbi_name + '\t' + strain_name + '\t' + show_tax)
            if a == 0:
                genome_ref = k
            a += 1

        self.logger.info(
            "Best matched genome is predicted to be: " + genome_ref)

        return genome_ref
