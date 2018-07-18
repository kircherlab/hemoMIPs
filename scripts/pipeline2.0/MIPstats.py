#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *08.04.2014
:Type: tool
:Input: BAM
:Output: TSV

"""

import sys,os
import pysam
from collections import defaultdict
from optparse import OptionParser,OptionGroup

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outfile", dest="outfile", help="Name of output file (def MIPstats.tsv)",default="MIPstats.tsv")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
(options, args) = parser.parse_args()

if options.verbose:
  sys.stderr.write("Opening input and output files.\n")

fileflags = 'wb'

if options.pipe: 
  infile = pysam.Samfile( "-", 'rb' )
  outfile = sys.stdout
else:
  if len(args) > 0 and os.path.exists(args[0]):
    infile = pysam.Samfile( args[0], 'rb' )
    outfile = open(options.outfile,'w')
  elif len(args) > 0:
    sys.stderr.write("Input file (%s) does not exist.\n"%(args[0]))
    sys.exit()
  else:
    sys.stderr.write("No input file specified.\n")
    sys.exit()

allsamples = set()
countDict = {}

for read in infile:
  if (read.is_paired and read.is_read1) or (not read.is_paired): 
    sample,mip = None,None
    for key,value in read.tags:
      if key == "RG": sample = value
      elif key == "ZM": mip = value
    if sample != None and mip != None:
      if mip not in countDict: 
        countDict[mip]=defaultdict(int)
      countDict[mip][sample]+=1
      allsamples.add(sample)

allsamples = list(allsamples)
allsamples.sort()

#sampleTotals = defaultdict(int)
#for mip,counts in countDict.iteritems():
  #for sample in allsamples:
    #sampleTotals[sample]+=counts[sample]

outfile.write("#MIP\t%s\n"%("\t".join(allsamples)))
for mip,counts in sorted(countDict.iteritems()):
  outfile.write("%s\t%s\n"%(mip,"\t".join(map(lambda sample: "%d"%(counts[sample]),allsamples)))) #"%.6f"%(counts[sample]/float(sampleTotals[sample])
