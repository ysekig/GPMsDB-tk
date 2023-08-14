#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'GPL3.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import logging
from GPMsDB_tk.calc import CalcHit


class AdjustMZ():
    def __init__(self):
        self.logger = logging.getLogger('GPMsDB_tk')

    def run(self, peaks, range_ms, num, db, tax):
        self.logger.info(
            "m/z adjustment using the ribosomal protein database.")

        scan = float(range_ms / num)

        c = {}
        d = {}
        q = {}
        r = {}
        result = {}

        for bin in range(0, num + 1):
            c[bin] = {}
            d[bin] = {}
            q[bin] = 0

            c[bin], d[bin] = CalcHit(self, peaks, scan, 200, bin, db, "adjust")

            result_init = sorted(
                c[bin].items(), key=lambda x: x[1], reverse=True)[:int(1)]
            result[bin] = dict(result_init)

            self.logger.info("---itaration " + str(1 + bin) +
                             ": +" + str(bin * scan) + " ppm")
            self.logger.info('\t' + "genome id" + '\t' + "n of hits" +
                             '\t' + "average deviation (ppm)" + "\t" + "organism")
            r[bin] = int(bin * scan)
            for l in result[bin].keys():
                try:
                    tax_line = tax[l]
                except:
                    tax_line = "not defined"
                self.logger.info('\t' + str(l) + '\t' + str(c[bin][l]) + '\t' + '{:.1f}'.format(
                    d[bin][l] / c[bin][l]) + "\t" + str(tax_line))
                q[bin] += c[bin][l]

        for bin2 in range(1, num + 1):
            bin = bin2 + num
            bin3 = -1 * bin2
            r[bin] = bin3 * scan
            c[bin] = {}
            d[bin] = {}
            q[bin] = 0

            c[bin], d[bin] = CalcHit(
                self, peaks, scan, 200, bin3, db, "adjust")

            result_init = sorted(
                c[bin].items(), key=lambda x: x[1], reverse=True)[:int(1)]
            result[bin] = dict(result_init)

            self.logger.info("---itaration " + str(1 + bin) +
                             ": " + str(r[bin]) + " ppm",)
            self.logger.info('\t' + "genome id" + '\t' + "n of hits" +
                             '\t' + "average deviation (ppm)" + "\t" + "organism")
            for l in result[bin]:
                try:
                    tax_line = tax[l]
                except:
                    tax_line = "not defined"

                self.logger.info('\t' + str(l) + '\t' + str(c[bin][l]) + '\t' + '{:.1f}'.format(
                    d[bin][l] / c[bin][l]) + "\t" + str(tax_line))
                q[bin] += c[bin][l]

        final_result = sorted(
            q.items(), key=lambda x: x[1], reverse=True)[:int(1)]
        f_result = dict(final_result)
        self.logger.info("---Final result ")
        for l in f_result:
            a = 1
            for n in result[l]:
                self.logger.info(
                    "\tBest adjustment: " + '{:.1f}'.format(r[l] + (d[l][n] / c[l][n])) + " ppm, " + str(q[l]))
                adjust = r[l] + (d[l][n] / c[l][n])
                if a == 1:
                    break

        return adjust
