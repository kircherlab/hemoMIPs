#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *20.01.2012

"""

import sys,os
import string
import pysam
from collections import defaultdict

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] BAMfile")
parser.add_option("-1","--verbose1",dest="verbose1",help="Print quality filtered sequences",default=False,action='store_true')
parser.add_option("-2","--verbose2",dest="verbose2",help="Print all completely unknown pairs",default=False,action='store_true')
parser.add_option("--no_A_extend",dest="noA",help="Do not extend short indexes by A",default=False,action='store_true')
parser.add_option("--desperate",dest="desperate",help="Check also for complemented and reverse complemented sequences",default=False,action="store_true")
parser.add_option("--max_lines",dest="maxlines",help="Maximum number of lines in output (default 10)",default=10,type="int")
parser.add_option("--454",dest="is454",help="Check for converted 454 indexes",default=False,action="store_true")
parser.add_option("--microRNA",dest="isMicroRNA",help="Check for microRNA indexes",default=False,action="store_true")
parser.add_option("--noQS",dest="isQS",help="Do not print number of quality filtered over total number of unknown",default=True,action="store_false")
(options, args) = parser.parse_args()

table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv') # COMPLEMENT DNA

def revcompl(seq):
  global table
  return seq.translate(table)[::-1]

def compl(seq):
  return seq.translate(table)[::-1]

if options.is454:
  P7 = ['/mnt/solexa/bin/IndexSeqs/35_nea_ordered_7nt_index.txt']
  P5 = []
elif options.isMicroRNA:
  P7 = ['/mnt/solexa/bin/IndexSeqs/wolfi_smallRNA_5pIndex.txt']
  P5 = []
else:
  P7 = ['/mnt/solexa/bin/IndexSeqs/all_P7_indexes.txt']
  P5 = ['/mnt/solexa/bin/IndexSeqs/96_5p_indexes.txt']

nucleotides = "ACGT"
P5_dict = {}
P7_dict = {}

def is_dna(seq):
  valid = True
  for elem in seq:
    if elem not in nucleotides:
      valid=False
      break
  return valid

for filename in P5:
  if not os.path.exists(filename): continue
  infile = open(filename)
  for line in infile:
    fields = line.split()
    if len(fields) >= 2 and is_dna(fields[0]):
      seq = fields[0].upper()
      name = fields[1]
      if seq in P5_dict:
        if name not in P5_dict[seq]:
          P5_dict[seq]+=":"+name
      else:
        P5_dict[seq]=name
  infile.close()
  if (not options.is454) and ("AGATCTC" not in P5_dict):
    P5_dict["AGATCTC"] = "IS4"

for filename in P7:
  if not os.path.exists(filename): continue
  infile = open(filename)
  for line in infile:
    fields = line.split()
    if len(fields) >= 2 and is_dna(fields[0]):
      seq = fields[0].upper()
      if len(seq) == 6 and not options.noA: seq+="A"
      name = fields[1]
      if seq in P7_dict:
        if name not in P7_dict[seq]:
          P7_dict[seq]+=":"+name
      else:
        P7_dict[seq]=name
  infile.close()
  if options.is454 and ("CAAGGCA" not in P7_dict):
    P7_dict["CAAGGCA"] = "454-B-adapter"

#print P7_dict
#print P5_dict

for filename in args:
  if os.path.exists(filename):
    print "Reading %s..."%filename
    cbamfile = pysam.Samfile(filename, "rb" )

    counter_unknown = 0
    counter_quality = 0
    counter_total = 0
    counter_wrong = 0

    counts = defaultdict(int)
    wrong = defaultdict(int)

    for read in cbamfile:
      counter_total += 1
      if read.is_qcfail: 
        counter_quality += 1
        if not (options.verbose1 or options.verbose2): continue

      is_unknown = False
      is_wrong = False
      index1,index2 = None,None
      for (key,value) in read.tags:
        if key == "RG" and value == "unknown":
          is_unknown = True
          counter_unknown += 1
        elif key == "RG" and value == "wrong":
          is_wrong = True
          counter_wrong += 1
        elif key == "XI": index1 = value
        elif key == "XJ": index2 = value

      if not is_unknown and not is_wrong: continue

      if index1 != None and index2 != None:
        ind1, ind2 = None, None
        if index1 in P7_dict: ind1 = P7_dict[index1]
        elif options.desperate and revcompl(index1) in P7_dict: ind1 = "RevCompl:"+P7_dict[revcompl(index1)]
        elif options.desperate and compl(index1) in P7_dict: ind1 = "Compl:"+P7_dict[compl(index1)]
        elif options.desperate and index1[::-1] in P7_dict: ind1 = "Reverse:"+P7_dict[index1[::-1]]
        if index2 in P5_dict: ind2 = P5_dict[index2]
        elif options.desperate and revcompl(index2) in P5_dict: ind2 = "RevCompl:"+P5_dict[revcompl(index2)]
        elif options.desperate and compl(index2) in P5_dict: ind2 = "Compl:"+P5_dict[compl(index2)]
        elif options.desperate and index2[::-1] in P5_dict: ind2 = "Reverse:"+P5_dict[index2[::-1]]

        if (ind1 == None or ind2 == None) and not options.verbose2: continue
        if is_unknown: counts[(index1,index2,ind1,ind2)] += 1
        if is_wrong: wrong[(index1,index2,ind1,ind2)] += 1
      elif index1 != None:
        ind1 = None
        if index1 in P7_dict: ind1 = P7_dict[index1]
        elif options.desperate and revcompl(index1) in P7_dict: ind1 = "RevCompl:"+P7_dict[revcompl(index1)]
        elif options.desperate and compl(index1) in P7_dict: ind1 = "Compl:"+P7_dict[compl(index1)]
        elif options.desperate and index1[::-1] in P7_dict: ind1 = "Reverse:"+P7_dict[index1[::-1]]

        if ind1 == None and not options.verbose2: continue
        if is_unknown: counts[(index1,ind1)] += 1
        if is_wrong: wrong[(index1,ind1)] += 1

    allcounts = map(lambda (x,y):(y,x),counts.iteritems())
    allcounts.sort()

    count = 0
    print "Perfect unknown indexes:"
    for obs,names in allcounts[::-1]:
      if (options.maxlines > 0) and (count > options.maxlines): break
      count += 1
      if len(names) == 4:
        sys.stdout.write("%7d\t%20s\t%s,%s\n"%(obs,names[0]+","+names[1],names[2],names[3]))
      else:
        sys.stdout.write("%7d\t%10s\t%s\n"%(obs,names[0],names[1]))
    if counter_wrong > 0:
      print "Perfect wrong index pairs:"
      allcounts = map(lambda (x,y):(y,x),wrong.iteritems())
      allcounts.sort()
      count = 0
      for obs,names in allcounts[::-1]:
        if (options.maxlines > 0) and (count > options.maxlines): break
        count += 1
        if len(names) == 4:
          sys.stdout.write("%7d\t%20s\t%s,%s\n"%(obs,names[0]+","+names[1],names[2],names[3]))
        else:
          sys.stdout.write("%7d\t%10s\t%s\n"%(obs,names[0],names[1]))

    if options.isQS: sys.stdout.write("%d unknown; %d wrong; %d out of %d quality filtered (%.2f%%)\n"%(counter_unknown,counter_wrong, counter_quality,counter_total,counter_quality/float(counter_total)*100))
