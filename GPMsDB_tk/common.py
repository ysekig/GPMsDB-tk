#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2022 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'CC BY-NC-SA 4.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import os
import errno
import sys
import logging
import time
import ntpath

import GPMsDB_tk
from GPMsDB_tk.defaultValues import DefaultValues


def PeakLoader(inputFile, minimum):
    checkFileExists(inputFile)

    t_peak = 0
    intens = {}
    com = ''
    for line in open(inputFile):
        if "#" in line:
            continue

        lineSplit = line.split()
        if lineSplit[0] == '':
            continue
        elif 'COM=' in lineSplit[0]:
            com = lineSplit[0].strip().replace("COM=", "")
        else:
            try:
                peak = float(lineSplit[0])
            except:
                continue
            try:
                intens[peak] = float(lineSplit[1])
                t_peak += intens[peak]
            except:
                t_peak += 1
                intens[peak] = 1

    p_use = 0
    total = 0
    peaks = {}

    result = sorted(intens.items(), reverse=False)
    dic = dict(result)

    for x in dic:
        total += 1
        if float(intens[x]) / float(t_peak) > minimum:
            peaks[x] = (intens[x]) / float(t_peak)
            p_use += 1

    return peaks, total, p_use, com


def selectDb(ref, tax, lineage):
    if lineage == "prokaryotes":
        if tax == "genome":
            tax_db = DefaultValues.TAX_GT_PR
        elif tax == "silva":
            tax_db = DefaultValues.TAX_SILVA_PR
        elif tax == "ncbi":
            tax_db = DefaultValues.TAX_NCBI_PR

        if ref == 'reps':
            db_rep = DefaultValues.REPS_REPS_DB_PR
            db_all = DefaultValues.REPS_ALL_DB_PR
            no_genes = DefaultValues.REPS_GENE_PR
        if ref == 'all':
            db_rep = DefaultValues.ALL_REPS_DB_PR
            db_all = DefaultValues.ALL_ALL_DB_PR
            no_genes = DefaultValues.GENE_PR
        if ref == 'custom':
            db_rep = DefaultValues.ALL_REPS_DB_PR
            db_all = DefaultValues.ALL_ALL_DB_PR
            no_genes = DefaultValues.GENE_PR
        strain_list = DefaultValues.STRAIN_DB_PR
    else:
        if tax == "genome":
            tax_db = DefaultValues.TAX_GT_FG
        elif tax == "silva":
            tax_db = DefaultValues.TAX_SILVA_FG
        elif tax == "ncbi":
            tax_db = DefaultValues.TAX_NCBI_FG

        if ref == 'reps':
            db_rep = DefaultValues.REPS_REPS_DB_FG
            db_all = DefaultValues.REPS_ALL_DB_FG
            no_genes = DefaultValues.REPS_GENE_FG
        if ref == 'all':
            db_rep = DefaultValues.ALL_REPS_DB_FG
            db_all = DefaultValues.ALL_ALL_DB_FG
            no_genes = DefaultValues.GENE_FG
        if ref == 'custom':
            db_rep = DefaultValues.ALL_REPS_DB_FG
            db_all = DefaultValues.ALL_ALL_DB_FG
            no_genes = DefaultValues.GENE_FG
        strain_list = DefaultValues.STRAIN_DB_FG

    return tax_db, db_rep, db_all, strain_list, no_genes


def checkEmptyDir(inputDir):
    if not os.path.exists(inputDir):
        makeSurePathExists(inputDir)
    else:
        files = os.listdir(inputDir)
        if len(files) != 0:
            logger = logging.getLogger('GPMsDB_tk')
            logger.error('Output directory must be empty: ' + inputDir)
            sys.exit(1)


def checkFileExists(inputFile):
    if not os.path.exists(inputFile):
        logger = logging.getLogger('GPMsDB_tk')
        logger.error('Input file does not exists: ' + inputFile)
        sys.exit(1)


def checkFileExistsNoBreak(inputFile):
    if not os.path.exists(inputFile):
        logger = logging.getLogger('GPMsDB_tk')
        logger.error('Input file does not exists: ' + inputFile)
        return "1"
    return "0"


def checkDirExists(inputDir):
    if not os.path.exists(inputDir):
        logger = logging.getLogger('GPMsDB_tk')
        logger.error('Input directory does not exists: ' + inputDir)
        sys.exit(1)


def makeSurePathExists(path):
    if not path:
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger('GPMsDB_tk')
            logger.error('Specified path does not exist: ' + path)
            sys.exit(1)


def restoreStdOut(outFile, oldStdOut):
    if (outFile != ''):
        try:
            sys.stdout.close()
            sys.stdout = oldStdOut
        except:
            logger = logging.getLogger('GPMsDB_tk')
            logger.error("Error restoring stdout ", outFile)
            sys.exit(1)


def genomeIdFromFilename(filename):
    genId = os.path.basename(filename)
    genId = os.path.splitext(genId)[0]

    return genId


def checkDbDir(inputDir):
    if not os.path.exists(inputDir):
        makeSurePathExists(inputDir)
    else:
        files = os.listdir(inputDir)
        if len(files) != 0:
            logger = logging.getLogger('GPMsDB_tk')
            logger.error(
                'Database (annotation files dir) is empty. Please set GPMsDB files correctly')
            sys.exit(1)


class StopWatch():
    def __init__(self, logger):
        self.time_start = time.time()
        self.time_latest = self.time_start
        self.logger = logger

    def clear(self):
        self.time_start = time.time()
        self.time_latest = self.time_start

    def lap(self):
        now = time.time()
        lap = now - self.time_latest
        total = now - self.time_start
        lap2 = int(lap + 0.5)
        h = lap2 // 3600 
        m = (lap2 - h * 3600) // 60
        s = lap2 - h * 3600 - m * 60
        total2 = int(total + 0.5)
        h2 = total2 // 3600 
        m2 = (total2 - h2 * 3600) // 60
        s2 = total2 - h2 * 3600 - m2 * 60
        self.logger.info(
            f"\n {{lap time: {h:02}:{m:02}:{s:02}, total time: {h2:02}:{m2:02}:{s2:02}}}")
        self.time_latest = now


def version():
    versionFile = open(os.path.join(GPMsDB_tk.__path__[0], 'VERSION'))
    return versionFile.readline().strip()


def logger_init(logger, output_dir=None, filename="GPMsDB-tk.log", silent=False):
    GPMsDB_tk_logger = logger
    GPMsDB_tk_logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")
    stream_logger = logging.StreamHandler(sys.stdout)
    stream_logger.setFormatter(log_format)
    GPMsDB_tk_logger.addHandler(stream_logger)
    if silent:
        GPMsDB_tk_logger.is_silent = True
        stream_logger.setLevel(logging.ERROR)

    if output_dir != None:
        os.makedirs(output_dir, exist_ok=True)
        timestamp_file_logger = logging.FileHandler(
            os.path.join(output_dir, filename), 'a')
        timestamp_file_logger.setFormatter(log_format)
        GPMsDB_tk_logger.addHandler(timestamp_file_logger)

    GPMsDB_tk_logger.info('%s v%s' % ("GPMsDB-tk", version()))
    GPMsDB_tk_logger.info(ntpath.basename(
        sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))
