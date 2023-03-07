#!/usr/bin/env python

__author__ = 'Yuji Sekiguchi'
__copyright__ = 'Copyright (c) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)'
__credits__ = ['Yuji Sekiguchi']
__license__ = 'CC BY-NC-SA 4.0'
__maintainer__ = 'Yuji Sekiguchi'
__email__ = 'y.sekiguchi@aist.go.jp'
__status__ = 'Development'

import random
from GPMsDB_tk.defaultValues import DefaultValues

cpdef tuple CalcHit(self, list peaks, double scan, int ppm, int bin, dict db, str s_type):
    cdef:
      int b
      str i
      int x
      int n
      double r
      double rexact
      double upper_peak
      double lower_peak
      list result_init
      dict result
      dict c = {}
      dict d = {}
      dict e = {}
      str double_count = "No" #"Yes"
      float peak

    for genome_id, genome_peaks in db.items():
        b = 0
        rexact = 0
        rexact0 = 0
        last_index = 0
        o = len(genome_peaks)
        for index, peak in enumerate(peaks):
              upper_peak = float(peak) + (float(peak) * ppm /1000000) + (float(peak) * scan * bin /1000000)
              lower_peak = float(peak) - (float(peak) * ppm /1000000) + (float(peak) * scan * bin /1000000)
              for n in range(last_index, len(genome_peaks)):
                  genome_peak = genome_peaks[n]
                  if float(lower_peak) < float(genome_peak) < float(upper_peak):
                      b += 1
                      last_index = n + 1
                      if s_type == "ms":
                          r = (abs(float(genome_peak) - float(peak) + (float(peak) * scan * bin /1000000)) / float(peak) * 1000000)
                          rexact += 1 - (0.5 * r / ppm)
                      else:
                          pass
                      if double_count == "Yes":
                          pass
                      elif double_count == "No":
                          break
                  elif float(genome_peak) <= float(lower_peak):
                      last_index = n
                      continue
                  elif float(genome_peak) >= float(upper_peak):
                      break

        c[genome_id] = b
        d[genome_id] = rexact

    return c,d


cpdef tuple CalcRamdom(self, list ramdom_list, list peaks, double scan, int ppm, int bin, dict db, str s_type):
    cdef:
      int b
      str i
      int x
      int n
      int value
      double r
      double rexact
      double upper_peak
      double lower_peak
      list result_init
      dict result
      dict c = {}
      dict d = {}
      dict e = {}
      str gen
      str double_count = "No" #"Yes"
      float peak

    for ind,(genome_id, _) in enumerate(ramdom_list):
        if ind > 199:
          break
        genome_peaks = db[genome_id]
        b = 0
        rexact = 0
        last_index = 0
        o = len(genome_peaks)
        for index, peak in enumerate(peaks):
          upper_peak = float(peak) + (float(peak) * ppm /1000000) + (float(peak) * scan * bin /1000000)
          lower_peak = float(peak) - (float(peak) * ppm /1000000) + (float(peak) * scan * bin /1000000)
          for n in range(last_index, len(genome_peaks)):
              genome_peak = genome_peaks[n]
              if float(lower_peak) < float(genome_peak) < float(upper_peak):
                  b += 1
                  last_index = n + 1
                  if s_type == "ms":
                      r = (abs(float(genome_peak) - float(peak) + (float(peak) * scan * bin /1000000)) / float(genome_peak) * 1000000)
                      rexact += 1 - (0.5 * r / ppm)
                  else:
                      pass
                  if double_count == "Yes":
                      pass
                  elif double_count == "No":
                      break
              elif float(genome_peak) <= float(lower_peak):
                  last_index = n
                  continue
              elif float(genome_peak) >= float(upper_peak):
                  break
          

        c[genome_id] = b
        d[genome_id] = rexact
        e[genome_id] = o

    return c, d


cpdef dict CalcScore(self, str score_type, dict reps, dict all, dict ms_reps, dict ms_all, dict genes):
        cdef:
          str k
          str i
          double j
          dict h = {}
          int total
          int gl = DefaultValues.GENE_LIMIT
          int rw = DefaultValues.RIBOSOMAL_WEIGHT

        if score_type == 'weighted':
          for k in all.keys():
            try:
              if genes[k] < gl:
                  j = gl
              else:
                  j = genes[k]
              h[k] = (rw * float(reps[k]) + float(all[k])) / j
            except KeyError:
              continue

        elif score_type == 'ms':
          for k in all.keys():
            try:
              if genes[k] < gl:
                  j = gl
              else:
                  j = genes[k]
              h[k] = (rw * float(ms_reps[k]) + float(ms_all[k])) / j
            except KeyError:
              continue

        elif score_type == 'unweighted':
          for k in all.keys():
            try:
              if genes[k] < gl:
                  j = gl
              else:
                  j = genes[k]
              h[k] = (float(reps[k]) + float(all[k])) / j
            except KeyError:
              continue

        return h

cpdef dict CalcScore2(self, str score_type, dict reps, dict all, dict ms_reps, dict ms_all, dict genes):
        cdef:
          str k
          str i
          dict h = {}
          int total
          int rw = DefaultValues.RIBOSOMAL_WEIGHT

        if score_type == 'weighted':
          for k in all.keys():
            try:
              h[k] = (rw * float(reps[k]) + float(all[k])) / genes[k]
            except KeyError:
              continue

        elif score_type == 'ms':
          for k in all.keys():
            try:
              h[k] = (rw * float(ms_reps[k]) + float(ms_all[k])) / genes[k]
            except KeyError:
              continue

        elif score_type == 'unweighted':
          for k in all.keys():
            try:
              h[k] = ((float(reps[k]) + float(all[k])) / genes[k])
            except KeyError:
              continue

        return h
