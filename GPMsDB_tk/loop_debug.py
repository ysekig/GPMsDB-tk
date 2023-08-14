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
import multiprocessing as mp

from GPMsDB_tk import __version__
from GPMsDB_tk.defaultValues import DefaultValues
from GPMsDB_tk.peakparser import PeakParser
from GPMsDB_tk.searchbest_linear import SearchBestHit2
from GPMsDB_tk.common import PeakLoader
from GPMsDB_tk.adjustmz import AdjustMZ
from GPMsDB_tk.common import makeSurePathExists, checkFileExists, checkFileExistsNoBreak


def version():
    import GPMsDB_tk_ident
    versionFile = open(os.path.join(GPMsDB_tk_ident.__path__[0], 'VERSION'))
    return versionFile.readline().strip()

class Loop2(object):
  def __init__(self):
      self.defaultout = "out"
      self.logger = logging.getLogger('GPMsDB_tk')

  def __workerThread(self, queueIn, queueOut, out_dir, auto_adjust, ppm_range, number_of_bins,
                     reference, ppm, score_type, minimum, tax_adjust, reps, all, genes, no_gen, tax):
      while True:
            in_file = queueIn.get(block=True, timeout=None)
            if in_file == None:
                break

            basename_without_ext = os.path.splitext(os.path.basename(in_file))[0]
            out_file = str(basename_without_ext + ".txt")
            outfile = os.path.join(out_dir, out_file)

            a = checkFileExistsNoBreak(in_file)
            if a == "1":
                queueOut.put(in_file)
                continue

            peaks, t_peak, p_use, com = PeakLoader(in_file, minimum)

            if len(list(peaks)) == 0:
                queueOut.put(in_file)
                continue

            list_peaks = []
            for i in peaks.keys():
                list_peaks.append(i)

            p = AdjustMZ()
            if auto_adjust == True:
                adjust = p.run(list_peaks,
                        ppm_range,
                        number_of_bins,
                        reps,
                        tax_adjust)
            else:
                adjust = 0

            p = SearchBestHit2()
            scores = p.run(in_file,
                     reference,
                     list_peaks,
                     ppm,
                     "",
                     "",
                     score_type,
                     adjust,
                     reps,
                     all,
                     tax,
                     tax_adjust,
                     com,
                     genes,
                     t_peak,
                     p_use,
                     minimum,
                     no_gen)

            result = sorted(scores.items(), key=lambda x:x[1], reverse=True)[:int(10)]
            dic = dict(result)

            with open(outfile, mode='w') as f:
                for n in dic:
                    f.write(n + "\t" + str(scores[n]) + "\t" + tax[n] + "\n")

            queueOut.put(in_file)

  def __writerThread(self, numDataItems, writerQueue):
      processedItems = 0
      while True:
            a = writerQueue.get(block=True, timeout=None)
            if a is None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems, numDataItems, float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
      sys.stdout.flush()
      sys.stdout.write('\n')

  def run(self, input_list, out_dir, auto_adjust, ppm_range, number_of_bins,
          reference, ppm, score_type, core, minimum, tax_adjust, reps, all,
          genes, no_gen, tax):

      peaklist_files = []
      for line in open(input_list):
        if line.startswith("#"):
            continue
        elif line == "":
            continue

        element = line.split("\t")
        peaklist_files.append(element[0].strip())

      print('  Number of unprocessed peak lists: %d' % len(peaklist_files))
      if not os.path.exists(out_dir):
          os.mkdir(out_dir)

      workerQueue = mp.Queue()
      writerQueue = mp.Queue()

      for f in peaklist_files:
        workerQueue.put(f)

      for _ in range(core):
        workerQueue.put(None)

      try:
        workerProc = [mp.Process(target = self.__workerThread, args = (workerQueue, writerQueue,
           out_dir, auto_adjust, ppm_range, number_of_bins, reference, ppm, score_type, minimum,
           tax_adjust, reps, all, genes, no_gen, tax)) for _ in range(core)]
        writeProc = mp.Process(target = self.__writerThread, args = (len(peaklist_files), writerQueue))

        writeProc.start()

        for p in workerProc:
          p.start()

        for p in workerProc:
          p.join()

        writerQueue.put(None)
        writeProc.join()
      except:
        for p in workerProc:
          p.terminate()

        writeProc.terminate()
