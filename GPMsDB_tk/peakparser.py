#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import os
import sys
import logging

import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})

from GPMsDB_tk.common import PeakLoader
from GPMsDB_tk.common import makeSurePathExists, checkFileExists


class PeakParser(object):
    def __init__(self):
        self.exr = '_detected.tsv'
        self.figdpi = 200
        self.logger = logging.getLogger('GPMsDB_tk')

    def run(self, input_file, reference, out_dir, ppm, genome_ref, adjust, filetype, db):
        if filetype.lower() == "pdf":
            figdpi = 72
        else:
            figdpi = self.figdpi

        file_name = str(genome_ref + '_annotation.tsv')
        reference_file = os.path.join(db, file_name)
        if not os.path.exists(reference_file):
            self.logger.error('Annotation file not exist for ' + file_name)
            return

        basename = os.path.basename(input_file)
        basename_without_ext = os.path.splitext(
            os.path.basename(input_file))[0]

        outfile = os.path.join(out_dir, basename_without_ext +
                               "_annotationwith_" + str(genome_ref) + "." + filetype)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        intens, t, tp, com = PeakLoader(input_file, 0)
        if len(list(intens)) == 0:
            self.logger.error('Peaks not found.')
            sys.exit(1)

        if com == "":
            com = "not specified"

        data = {}
        for line in open(reference_file):
            lineSplit = line.split('\t')
            if lineSplit[0] == '#' or lineSplit[0] == '':
                continue

            else:
                peak = lineSplit[0].rstrip()
                data[peak] = str(lineSplit[1].rstrip() + "\t" + peak)

        fd = {}
        for i in intens.keys():
            upper_peak = float(i) + (float(i) * ppm / 1000000) + \
                (float(i) * adjust / 1000000)
            lower_peak = float(i) - (float(i) * ppm / 1000000) + \
                (float(i) * adjust / 1000000)
            for k in data.keys():
                if lower_peak < float(k) < upper_peak:
                    try:
                        fd[i].append(data[k])
                    except:
                        fd[i] = []
                        fd[i] = data[k]

        ri = {}
        ri2 = {}
        for i in intens.keys():
            ri[i] = intens[i] * 100
            if i in fd.keys():
                ri2[i] = intens[i] * 100

        self.logger.info("peak m/z\tintensity\tgene annotation\tmw")
        for i in fd:
            self.logger.info(
                str(i) + "\t" + "{:.3f}".format(ri[i]) + "\t" + fd[i])

        # generating graphs
        x = list(ri.keys())
        y = list(ri.values())

        x1 = list(ri2.keys())
        y2 = list(ri2.values())

        title = str("Peak annotation for " + str(basename) + "\n" +
                    "Total number of peaks: " + str(len(intens)) + "\n" +
                    "Reference genome: " + str(genome_ref))

        fig, ax = plt.subplots()
        plt.figure(figsize=(15, 10)).patch.set_facecolor('white')
        #plt.suptitle(title, fontsize = 13, x=fig.subplotpars.left, ha = 'left')
        plt.rcParams["font.size"] = 12
        #plt.rcParams['font.family'] = 'sans-serif'
        #plt.rcParams['font.sans-serif'] = ['Arial']

        ax1 = plt.subplot(2, 1, 1)

        if not len(intens) == 0:
            ax1.stem(x, y, linefmt="k-", basefmt=" ",
                     markerfmt=" ", use_line_collection=True)
        if not len(ri2) == 0:
            ax1.stem(x1, y2, linefmt="C3-", basefmt=" ",
                     markerfmt=" ", use_line_collection=True)

        plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.5)
        ax1.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda x, loc: "{:,}".format(int(x))))
        plt.ylim(0,)
        plt.xlim(0, 15000)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.xlabel("m/z", fontsize=14)
        plt.ylabel("relative intensity (linear, %)", fontsize=14)

        ax2 = plt.subplot(2, 1, 2)

        plt.yscale('log')

        if not len(intens) == 0:
            ax2.stem(x, y, linefmt="k-", basefmt=" ",
                     markerfmt=" ", use_line_collection=True)
        if not len(ri2) == 0:
            ax2.stem(x1, y2, linefmt="C3-", basefmt=" ",
                     markerfmt=" ", use_line_collection=True)

        plt.grid(True, axis="y", color='black', linestyle=':', linewidth=0.5)
        ax2.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda x, loc: "{:,}".format(int(x))))
        plt.xlim(0, 15000)
        plt.xlabel("m/z", fontsize=14)
        plt.ylabel("relative intensity (log, %)", fontsize=14)

        result = sorted(fd.items(), key=lambda x: x[0], reverse=True)
        dic = dict(result)
        y_axis_min, y_axis_max = ax1.get_ylim()

        if len(fd) > 0:
            a = 0
            for k in dic:
                bar = (y_axis_max - ri[k]) / y_axis_max * figdpi * 3
                line = fd[k].split("\t")
                try:
                    if line[0] == "":
                        line[0] = "hypothetical protein"
                except:
                    line[0] = "hypothetical protein"
                gene = str(str(
                    k) + " m/z, " + line[0] + " (" + "{:.3f}".format(float(line[1])) + " theoretical m/z)")
                ax1.annotate(str(gene),
                             xy=(k, ri[k]),
                             xycoords='data', rotation=90,
                             xytext=(890 - a, 610),
                             textcoords='figure points', size=9,
                             arrowprops=dict(arrowstyle="-",
                                             relpos=(0.5, 0),
                                             connectionstyle="arc, angleA=-90, armA=" + "10" +
                                             ",angleB= 90, armB=" + str(bar) +
                                             ",rad=0", linewidth=0.5, color='#808080'),

                             horizontalalignment='center', verticalalignment='bottom')
                a += 10

        plt.savefig(outfile, dpi=figdpi, format=filetype,
                    bbox_inches='tight', pad_inches=0.1)
        plt.close()

        return
