#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *09.01.2015
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict 
import pysam

parser = OptionParser("%prog [options]")
parser.add_option("-p","--prefix", dest="prefix", help="Prefix for output filenames (default Length)",default="Length")
parser.add_option("-l","--library", dest="library", help="Use library name from RG read header rather than the sample ID",default=False,action="store_true")
parser.add_option("--max_length", dest="max_length", help="Maximum length considered for the output (default 300)",default=300,type="int")
parser.add_option("--binsize", dest="binsize", help="Bin size for length binning (default 5)",default=5,type="int")
(options, args) = parser.parse_args()

if options.library: options.all=True

rgroups = {}
binmiddle = options.binsize/2.0

rarray = [0,10,20,30]
for filename in args:
  if os.path.exists(filename):
    print "Reading %s..."%filename
    rgroups = {}
    cbamfile = pysam.Samfile(filename, "rb" )
    id2lib = {}
    if options.library and 'RG' in cbamfile.header:
      for rgroup in cbamfile.header['RG']:
        if 'LB' in rgroup and 'ID' in rgroup:
          id2lib[rgroup['ID']] = rgroup['LB']
    for read in cbamfile:
        library,count = '',1
        for (key,value) in read.tags:
          if key == "RG":
            if value in id2lib: library = id2lib[value]
            else: library = value
          elif key == "XP":
            count = value
        if library not in rgroups: 
          rgroups[library] = [defaultdict(int),defaultdict(float),defaultdict(int),defaultdict(int),defaultdict(int),defaultdict(int)] #  [Count,Ave,0,10,20,30]

        lseq = None
        if not read.is_paired:
          lseq = len(read.seq)
          lseq = (lseq//options.binsize)*options.binsize+binmiddle
        elif read.is_read1:
          lseq = min(options.max_length,abs(read.isize))
          if lseq == 0: lseq = options.max_length
          lseq = (lseq//options.binsize)*options.binsize+binmiddle

        if lseq != None:
          rgroups[library][0][lseq]+=count
          rgroups[library][1][lseq]+=read.mapq*count
          for ind,val in enumerate(rarray):
            if read.mapq <= val:
              rgroups[library][ind+2][lseq]+=count

for library in rgroups:
  outfile = open("%s_%s.tsv"%(options.prefix.rstrip("_"),library),'w')
  outfile.write('Length\tCount\tMapQ0\tMapQ10\tMapQ20\tMapQ30\tAveMapQ\n')
  for length in sorted(rgroups[library][0].keys()):
    outfile.write("%.1f\t%d\t%d\t%d\t%d\t%d\t%.2f\n"%(length,rgroups[library][0][length],rgroups[library][2][length],rgroups[library][3][length],rgroups[library][4][length],rgroups[library][5][length],rgroups[library][1][length]/float(rgroups[library][0][length])))

