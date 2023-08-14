#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2022 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import sys
import os
import logging
import pickle
import time

from GPMsDB_tk.common import (selectDb, PeakLoader,
                              makeSurePathExists, checkFileExists,
                              checkEmptyDir)
from GPMsDB_tk.defaultValues import DefaultValues
from GPMsDB_tk.adjustmz import AdjustMZ
from GPMsDB_tk.loop import Loop
from GPMsDB_tk.loop_debug import Loop2
from GPMsDB_tk.peakparser import PeakParser
from GPMsDB_tk.plot_peaks import PlotPeaks
from GPMsDB_tk.searchbest import SearchBestHit
from GPMsDB_tk.common import StopWatch, logger_init


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger('GPMsDB_tk')
        self.stopwatch = StopWatch(self.logger)

    def binFiles(self, binFolder, binExtension):
        binFiles = []
        if binFolder is not None:
            all_files = os.listdir(binFolder)
            for f in all_files:
                if f.endswith(binExtension):
                    binFile = os.path.join(binFolder, f)
                    if os.stat(binFile).st_size == 0:
                        self.logger.warning(
                            "Skipping genome %s as it has a size of 0 bytes." % f)
                    else:
                        binFiles.append(binFile)

        if not binFiles:
            self.logger.error(
                "No genomes found. Check the extension (-x) used to identify bins.")
            sys.exit(1)

        return sorted(binFiles)

    def inspect(self, options):
        logger_init(self.logger, None, silent=False)
        self.logger.info(
            '[inspect] Check for peak list file and generate peak plot.')

        checkFileExists(options.input_file)
        makeSurePathExists(options.out_dir)

        p = PlotPeaks()
        outfile = p.run(options.input_file,
                        options.out_dir,
                        options.filetype)

        self.logger.info("File: " + str(outfile) + " saved.")
        self.stopwatch.lap()

    def adjust(self, options):
        logger_init(self.logger, None, silent=options.silent)
        self.logger.info(
            '[adjust] Check for peak list file and generate peak plot.')

        checkFileExists(options.input_file)
        peaks, t_peak, p_use, com = PeakLoader(
            options.input_file, options.minimum)
        if len(list(peaks)) == 0:
            self.logger.error('Peaks not found.')
            sys.exit(1)

        list_peaks = []
        for i in peaks.keys():
            list_peaks.append(i)

        self.logger.info('[adjust] Loading databases.')
        reps_db = DefaultValues.REPS_REPS_DB
        tax_db = DefaultValues.TAX_NCBI
        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax = pickle.load(f)

        p = AdjustMZ()
        p.run(list_peaks,
              options.ppm_range,
              options.number_of_bins,
              reps,
              tax)

        self.stopwatch.lap()

    def identify(self, options):
        logger_init(self.logger, None, silent=options.silent)
        self.logger.info(
            '[identify] Search for the best matched genome based on given protein peak list')

        checkFileExists(options.input_file)
        peaks, t_peak, p_use, com = PeakLoader(
            options.input_file, options.minimum)
        if len(list(peaks)) == 0:
            self.logger.error('Peaks not found.')
            sys.exit(1)

        list_peaks = []
        for i in peaks.keys():
            list_peaks.append(i)

        self.logger.info('[identify] Loading databases.')
        tax_db, reps_db, all_db, strain_list, no_genes = selectDb(
            options.reference, options.taxonomy)

        ncbi_db = DefaultValues.TAX_NCBI
        with open(ncbi_db, 'rb') as f:
            ncbi = pickle.load(f)
        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(all_db, 'rb') as f:
            all = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax = pickle.load(f)
        with open(strain_list, 'rb') as f:
            strain = pickle.load(f)
        with open(no_genes, 'rb') as f:
            genes = pickle.load(f)

        if options.reference == 'custom':
            db_rep_c = DefaultValues.CUSTOM_LIST_R
            db_all_c = DefaultValues.CUSTOM_LIST_O
            no_genes_c = DefaultValues.CUSTOM_LIST_GENES
            strain_list_c = DefaultValues.CUSTOM_LIST_NAME
            tax_db_c = DefaultValues.CUSTOM_LIST_TAX
            with open(strain_list_c, 'rb') as f:
                ncbi_c = pickle.load(f)
            with open(db_rep_c, 'rb') as f:
                reps_c = pickle.load(f)
            with open(db_all_c, 'rb') as f:
                all_c = pickle.load(f)
            with open(tax_db_c, 'rb') as f:
                tax_c = pickle.load(f)
            with open(strain_list_c, 'rb') as f:
                strain_c = pickle.load(f)
            with open(no_genes_c, 'rb') as f:
                genes_c = pickle.load(f)
            ncbi.update(ncbi_c)
            reps.update(reps_c)
            all.update(all_c)
            tax.update(tax_c)
            genes.update(genes_c)

        p = SearchBestHit()
        p.run(options.input_file,
              options.reference,
              list_peaks,
              options.ppm,
              options.first,
              options.top,
              options.score_type,
              options.adjust,
              reps,
              all,
              tax,
              ncbi,
              strain,
              com,
              genes,
              options.taxonomy,
              t_peak,
              p_use,
              options.minimum)

        self.stopwatch.lap()

    def identify_wf(self, options):
        logger_init(self.logger, None, silent=options.silent)
        self.logger.info(
            '[identify_wf] Identify workflow (adjust > identify).')

        checkFileExists(options.input_file)
        peaks, t_peak, p_use, com = PeakLoader(
            options.input_file, options.minimum)
        if len(list(peaks)) == 0:
            self.logger.error('Peaks not found.')
            sys.exit(1)

        list_peaks = []
        for i in peaks.keys():
            list_peaks.append(i)

        self.logger.info('[identify_wf] Loading databases.')
        reps_db = DefaultValues.REPS_REPS_DB
        tax_db = DefaultValues.TAX_NCBI
        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax_adjust = pickle.load(f)

        tax_db, reps_db, all_db, strain_list, no_genes = selectDb(
            options.reference, options.taxonomy)

        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(all_db, 'rb') as f:
            all = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax = pickle.load(f)
        with open(strain_list, 'rb') as f:
            strain = pickle.load(f)
        with open(no_genes, 'rb') as f:
            genes = pickle.load(f)

        if options.reference == 'custom':
            db_rep_c = DefaultValues.CUSTOM_LIST_R
            db_all_c = DefaultValues.CUSTOM_LIST_O
            no_genes_c = DefaultValues.CUSTOM_LIST_GENES
            strain_list_c = DefaultValues.CUSTOM_LIST_NAME
            tax_db_c = DefaultValues.CUSTOM_LIST_TAX
            with open(db_rep_c, 'rb') as f:
                reps_c = pickle.load(f)
            with open(db_all_c, 'rb') as f:
                all_c = pickle.load(f)
            with open(tax_db_c, 'rb') as f:
                tax_c = pickle.load(f)
            with open(strain_list_c, 'rb') as f:
                strain_c = pickle.load(f)
            with open(no_genes_c, 'rb') as f:
                genes_c = pickle.load(f)
            reps.update(reps_c)
            all.update(all_c)
            tax.update(tax_c)
            tax_adjust.update(strain_c)
            genes.update(genes_c)

        p = AdjustMZ()
        if options.auto_adjust == True:
            adjust = p.run(list_peaks,
                           options.ppm_range,
                           options.number_of_bins,
                           reps,
                           tax_adjust)

            self.stopwatch.lap()
        else:
            adjust = options.adjust

        p = SearchBestHit()
        p.run(options.input_file,
              options.reference,
              list_peaks,
              options.ppm,
              options.first,
              options.top,
              options.score_type,
              adjust,
              reps,
              all,
              tax,
              tax_adjust,
              strain,
              com,
              genes,
              options.taxonomy,
              t_peak,
              p_use,
              options.minimum)

        self.stopwatch.lap()

    def identify_bwf(self, options):
        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][identify_bwf] Full identification workflow for a batch of files.", cnvtime))

        checkFileExists(options.input_list)
        makeSurePathExists(options.out_dir)
        checkEmptyDir(options.out_dir)

        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][identify_bwf] Loading databases.", cnvtime))
        reps_db = DefaultValues.REPS_REPS_DB
        tax_db = DefaultValues.TAX_NCBI
        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax_adjust = pickle.load(f)

        tax_db, reps_db, all_db, strain_list, no_genes = selectDb(
            options.reference, options.taxonomy)

        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(all_db, 'rb') as f:
            all = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax = pickle.load(f)
        with open(strain_list, 'rb') as f:
            strain = pickle.load(f)
        with open(no_genes, 'rb') as f:
            genes = pickle.load(f)

        if options.reference == 'custom':
            db_rep_c = DefaultValues.CUSTOM_LIST_R
            db_all_c = DefaultValues.CUSTOM_LIST_O
            no_genes_c = DefaultValues.CUSTOM_LIST_GENES
            strain_list_c = DefaultValues.CUSTOM_LIST_NAME
            tax_db_c = DefaultValues.CUSTOM_LIST_TAX
            with open(db_rep_c, 'rb') as f:
                reps_c = pickle.load(f)
            with open(db_all_c, 'rb') as f:
                all_c = pickle.load(f)
            with open(tax_db_c, 'rb') as f:
                tax_c = pickle.load(f)
            with open(strain_list_c, 'rb') as f:
                strain_c = pickle.load(f)
            with open(no_genes_c, 'rb') as f:
                genes_c = pickle.load(f)
            reps.update(reps_c)
            all.update(all_c)
            tax.update(tax_c)
            tax_adjust.update(strain_c)
            genes.update(genes_c)

        peakdetect = "no"
        filetype = "pdf"
        p = Loop()
        p.run(options.input_list,
              options.out_dir,
              options.auto_adjust,
              options.ppm_range,
              options.number_of_bins,
              options.reference,
              options.ppm,
              options.first,
              options.top,
              options.score_type,
              options.core,
              options.minimum,
              filetype,
              tax_adjust,
              reps,
              all,
              tax,
              strain,
              genes,
              options.taxonomy,
              peakdetect)

        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][identify_bwf] Finished.", cnvtime))

    def peak(self, options):
        logger_init(self.logger, options.out_dir, silent=options.silent)
        self.logger.info(
            '[peak] Annotate protein identity in a peak list based on given genome ID.')

        checkFileExists(options.input_file)
        makeSurePathExists(options.out_dir)

        db = DefaultValues.GENOME_DIR

        p = PeakParser()
        p.run(options.input_file,
              None,
              options.out_dir,
              options.ppm,
              options.genome_ref,
              options.adjust,
              options.filetype,
              db)

        self.stopwatch.lap()

    def peak_wf(self, options):
        self.logger.info(
            '[peak_wf] Full identification/peak annotation workflow (adjust > identify > peak).')
        logger_init(self.logger, None, silent=options.silent)

        checkFileExists(options.input_file)
        makeSurePathExists(options.out_dir)

        peaks, t_peak, p_use, com = PeakLoader(
            options.input_file, options.minimum)
        if len(list(peaks)) == 0:
            self.logger.error('Peaks not found.')
            sys.exit(1)

        list_peaks = []
        for i in peaks.keys():
            list_peaks.append(i)

        self.logger.info('[peak_wf] Loading databases.')
        reps_db = DefaultValues.REPS_REPS_DB
        ncbi = DefaultValues.TAX_NCBI
        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(ncbi, 'rb') as f:
            tax_adjust = pickle.load(f)

        tax_db, reps_db, all_db, strain_list, no_genes = selectDb(
            options.reference, options.taxonomy)

        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(all_db, 'rb') as f:
            all = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax = pickle.load(f)
        with open(strain_list, 'rb') as f:
            strain = pickle.load(f)
        with open(no_genes, 'rb') as f:
            genes = pickle.load(f)

        if options.reference == 'custom':
            db_rep_c = DefaultValues.CUSTOM_LIST_R
            db_all_c = DefaultValues.CUSTOM_LIST_O
            no_genes_c = DefaultValues.CUSTOM_LIST_GENES
            strain_list_c = DefaultValues.CUSTOM_LIST_NAME
            tax_db_c = DefaultValues.CUSTOM_LIST_TAX
            with open(db_rep_c, 'rb') as f:
                reps_c = pickle.load(f)
            with open(db_all_c, 'rb') as f:
                all_c = pickle.load(f)
            with open(tax_db_c, 'rb') as f:
                tax_c = pickle.load(f)
            with open(strain_list_c, 'rb') as f:
                strain_c = pickle.load(f)
            with open(no_genes_c, 'rb') as f:
                genes_c = pickle.load(f)
            reps.update(reps_c)
            all.update(all_c)
            tax.update(tax_c)
            tax_adjust.update(strain_c)
            genes.update(genes_c)

        p = AdjustMZ()
        if options.auto_adjust == True:
            adjust = p.run(list_peaks,
                           options.ppm_range,
                           options.number_of_bins,
                           reps,
                           tax_adjust)

            self.stopwatch.lap()
        else:
            adjust = options.adjust

        p = SearchBestHit()
        best = p.run(options.input_file,
                     options.reference,
                     list_peaks,
                     options.ppm,
                     options.first,
                     options.top,
                     options.score_type,
                     adjust,
                     reps,
                     all,
                     tax,
                     tax_adjust,
                     strain,
                     com,
                     genes,
                     options.taxonomy,
                     t_peak,
                     p_use,
                     options.minimum)

        self.stopwatch.lap()

        db = DefaultValues.GENOME_DIR

        p = PeakParser()
        p.run(options.input_file,
              options.reference,
              options.out_dir,
              options.ppm,
              best,
              adjust,
              options.filetype,
              db)

        self.stopwatch.lap()

    def batch_wf(self, options):
        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][peak_bwf] Full identification/peak annotation workflow for a batch of files.", cnvtime))

        checkFileExists(options.input_list)
        makeSurePathExists(options.out_dir)
        checkEmptyDir(options.out_dir)

        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][peak_bwf] Loading databases.", cnvtime))
        reps_db = DefaultValues.REPS_REPS_DB
        tax_db = DefaultValues.TAX_NCBI
        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax_adjust = pickle.load(f)

        tax_db, reps_db, all_db, strain_list, no_genes = selectDb(
            options.reference, options.taxonomy)

        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(all_db, 'rb') as f:
            all = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax = pickle.load(f)
        with open(strain_list, 'rb') as f:
            strain = pickle.load(f)
        with open(no_genes, 'rb') as f:
            genes = pickle.load(f)

        if options.reference == 'custom':
            db_rep_c = DefaultValues.CUSTOM_LIST_R
            db_all_c = DefaultValues.CUSTOM_LIST_O
            no_genes_c = DefaultValues.CUSTOM_LIST_GENES
            strain_list_c = DefaultValues.CUSTOM_LIST_NAME
            tax_db_c = DefaultValues.CUSTOM_LIST_TAX
            with open(db_rep_c, 'rb') as f:
                reps_c = pickle.load(f)
            with open(db_all_c, 'rb') as f:
                all_c = pickle.load(f)
            with open(tax_db_c, 'rb') as f:
                tax_c = pickle.load(f)
            with open(strain_list_c, 'rb') as f:
                strain_c = pickle.load(f)
            with open(no_genes_c, 'rb') as f:
                genes_c = pickle.load(f)
            reps.update(reps_c)
            all.update(all_c)
            tax.update(tax_c)
            tax_adjust.update(strain_c)
            genes.update(genes_c)

        peakdetect = "yes"

        p = Loop()
        p.run(options.input_list,
              options.out_dir,
              options.auto_adjust,
              options.ppm_range,
              options.number_of_bins,
              options.reference,
              options.ppm,
              options.first,
              options.top,
              options.score_type,
              options.core,
              options.minimum,
              options.filetype,
              tax_adjust,
              reps,
              all,
              tax,
              strain,
              genes,
              options.taxonomy,
              peakdetect)

        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][peak_bwf] Finished.", cnvtime))

    def debug(self, options):
        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][debugging] Full identification workflow for a batch of files.", cnvtime))

        checkFileExists(options.input_list)
        makeSurePathExists(options.out_dir)
        checkEmptyDir(options.out_dir)

        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][debugging] Loading databases.", cnvtime))
        reps_db = DefaultValues.REPS_REPS_DB
        tax_db = DefaultValues.TAX_NCBI
        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax_adjust = pickle.load(f)

        tax_db, reps_db, all_db, strain_list, no_genes = selectDb(
            options.reference, 'gtdb')

        with open(reps_db, 'rb') as f:
            reps = pickle.load(f)
        with open(all_db, 'rb') as f:
            all = pickle.load(f)
        with open(no_genes, 'rb') as f:
            genes = pickle.load(f)
        with open(tax_db, 'rb') as f:
            tax = pickle.load(f)

        if options.reference == 'custom':
            db_rep_c = DefaultValues.CUSTOM_LIST_R
            db_all_c = DefaultValues.CUSTOM_LIST_O
            no_genes_c = DefaultValues.CUSTOM_LIST_GENES
            with open(db_rep_c, 'rb') as f:
                reps_c = pickle.load(f)
            with open(db_all_c, 'rb') as f:
                all_c = pickle.load(f)
            with open(no_genes_c, 'rb') as f:
                genes_c = pickle.load(f)
            reps.update(reps_c)
            all.update(all_c)
            genes.update(genes_c)

        p = Loop2()
        p.run(options.input_list,
              options.out_dir,
              options.auto_adjust,
              options.ppm_range,
              options.number_of_bins,
              options.reference,
              options.ppm,
              options.score_type,
              options.core,
              options.minimum,
              tax_adjust,
              reps,
              all,
              genes,
              options.gene_number,
              tax)

        now = time.ctime()
        cnvtime = time.strptime(now)
        print(time.strftime(
            "[%Y-%m-%d %H:%M:%S][batch_wf] Finished.", cnvtime))

    def parse_options(self, options):
        if options.subparser_name == 'data':
            self.update_DB(options)
        elif options.subparser_name == 'adjust':
            self.adjust(options)
        elif options.subparser_name == 'inspect':
            self.inspect(options)
        elif options.subparser_name == 'identify':
            self.identify(options)
        elif options.subparser_name == 'identify_wf':
            self.identify_wf(options)
        elif options.subparser_name == 'identify_bwf':
            self.identify_bwf(options)
        elif options.subparser_name == 'peak_wf':
            self.peak_wf(options)
        elif options.subparser_name == 'peak_bwf':
            self.batch_wf(options)
        elif options.subparser_name == 'peak':
            self.peak(options)
        elif options.subparser_name == 'debug':
            self.debug(options)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0
