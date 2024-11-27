#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2022 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import os
import sys

class DefaultValues():
    try:
        GENERIC_PATH = os.environ['GPMsDB_PATH']
    except KeyError:
        print('  ERROR ')
        print("The 'GPMsDB_PATH' environment variable is not defined.")
        print('Please set this variable to your reference data package.' + '\n')
        sys.exit(1)

    NO_THREAD = 4           #number of default threads

    CHECK_RANGE = 1000      #range of torelance (ppm) to check
    NO_BIN = 5              #numbert of bins to be tested for given range of ppm.
    TORELANCE = 200         #torelance (ppm)'

    MIN_PEAK = 0.0002        #minimum peak relative abundance (0.001 as 0.1%%)'

    HIT_RETAIN_FST = 200    #number of hits retained in the 1st screening based on ribosomal proteins
    HIT_SHOW = 20          #number of top hits shown

    HIT_EXCLUDE_REP = 1000    #number of top hits excluded for calculating average scores in representative db
    HIT_EXCLUDE_ALL = 50000    #number of top hits excluded for calculating average scores in all db
    HIT_EXCLUDE_ALL2 = 400    #number of top hits excluded for calculating average scores in all db (2nd criteria)

    GENE_LIMIT = 800    #number of genes under which gene number normalization is canceled
    RIBOSOMAL_WEIGHT = 7      #weighting of ribosomal protein detection compared with other proteins

    GPMsDB_PATH = GENERIC_PATH

    REPS_REPS_DB = os.path.join(GPMsDB_PATH, 'mass', 'ribosomal_reps.db')
    REPS_ALL_DB = os.path.join(GPMsDB_PATH, 'mass', 'all_reps.db')
    ALL_REPS_DB = os.path.join(GPMsDB_PATH, 'mass', 'ribosomal.db')
    ALL_ALL_DB = os.path.join(GPMsDB_PATH, 'mass', 'all.db')
    REPS_GENE = os.path.join(GPMsDB_PATH, 'mass', 'reps_genes.db')
    ALL_GENE = os.path.join(GPMsDB_PATH, 'mass', 'all_genes.db')

    TAX_NCBI = os.path.join(GPMsDB_PATH, 'taxonomy', 'ncbi_name.db')
    STRAIN_DB = os.path.join(GPMsDB_PATH, 'taxonomy', 'ncbi_strain.db')
    TAX_GTDB = os.path.join(GPMsDB_PATH, 'taxonomy', 'gtdb_taxonomy.db')
    TAX_GG = os.path.join(GPMsDB_PATH, 'taxonomy', 'ssu_gg_taxonomy.db')
    TAX_SILVA = os.path.join(GPMsDB_PATH, 'taxonomy', 'ssu_silva_taxonomy.db')
    GENOME_DIR = os.path.join(GPMsDB_PATH, 'genomes')

    CUSTOM_LIST_R = os.path.join(GPMsDB_PATH, 'custom', 'custom_ribosomals.db')
    CUSTOM_LIST_O = os.path.join(GPMsDB_PATH, 'custom', 'custom_others.db')
    CUSTOM_LIST_GENES = os.path.join(GPMsDB_PATH, 'custom', 'custom_genes.db')
    CUSTOM_LIST_NAME = os.path.join(GPMsDB_PATH, 'custom', 'custom_names.db')
    CUSTOM_LIST_TAX = os.path.join(GPMsDB_PATH, 'custom', 'custom_taxonomy.db')
    
    PROBABIL_RIBOSOMAL = 5
    PROBABIL_REPS_90 = 0.00001
    PROBABIL_REPS_99 = 0.000001
    PROBABIL_ALL_90 =  0.000001
    PROBABIL_ALL_99 =  0.0000001
