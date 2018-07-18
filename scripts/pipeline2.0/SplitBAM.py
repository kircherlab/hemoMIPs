#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *15.01.2012
"""

import sys, os
from optparse import OptionParser
import pysam

parser = OptionParser("%prog [options] input.bam")
parser.add_option("-p","--PIPE", dest="pipe", help="Read BAM from stdin instead of file",default=False,action="store_true")
parser.add_option("--library", dest="library", help="Use library name from RG read header rather than the sample ID",default=False,action="store_true")
parser.add_option("--switchMainJIndex", dest="switch", help="Switch main read and second index sequence on the fly",default=False,action="store_true")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'Split').",default="Split")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.",default = None)
(options, args) = parser.parse_args()

if options.outdir != None and os.path.isdir(options.outdir):
  options.outdir = options.outdir.rstrip('/')+"/"
elif options.outdir == None:
  options.outdir = ""
else:
  sys.stderr.write('Error: No valid output directory defined\n')
  sys.exit()

if len(options.outprefix) > 0:
  options.outprefix = options.outprefix.rstrip('_')+'_'

filename = None
if options.pipe:
  sys.stderr.write("Reading from STDIN...\n")
  filename = "-"
elif len(args) > 0 and os.path.exists(args[0]):
  filename = args[0]
  sys.stderr.write("Reading %s...\n"%filename)
else:
  sys.stderr.write('Error: No valid input file defined\n')
  sys.exit()

rgroups = {}
id2lib = {}
outfiles = {}
cbamfile = pysam.Samfile(filename, "rb" )
if options.library and 'RG' in cbamfile.header:
  for rgroup in cbamfile.header['RG']:
    if 'LB' in rgroup and 'ID' in rgroup:
      id2lib[rgroup['ID']] = rgroup['LB']

for read in cbamfile:
  library,count = '',1
  if options.switch:
    hseq = read.seq
    hqual = read.qual
    iseq = "*"
    iqual = "*"
    newTags = []
    for (key,value) in read.tags:
      if key == "XJ": 
        iseq = value
        newTags.append(("XJ",hseq))
      elif key == "YJ": 
        newTags.append(("YJ",hqual))
        iqual = value
      else:
        newTags.append((key,value))
    read.seq = iseq
    read.qual = iqual
    read.tags = newTags
  for (key,value) in read.tags:
    if key == "RG":
      if value in id2lib: library = id2lib[value]
      else: library = value
    elif key == "XP":
      count = value
  if library in rgroups:
    if not read.is_paired or read.is_read1:
      rgroups[library][0]+=1
      rgroups[library][1]+=count
    rgroups[library][2]+=len(read.seq)
  else:
    if not read.is_paired or read.is_read1:
      rgroups[library]=[1,count,len(read.seq)]
    else:
      rgroups[library]=[1,count,0]
  if library in outfiles: outfiles[library].write(read)
  else:
    outfiles[library] = pysam.Samfile(options.outdir+options.outprefix+"%s.bam"%(library), "wb", template=cbamfile )
    outfiles[library].write(read)
cbamfile.close()

sys.stderr.write('Library\tReads\tReadsXP\tBases\n')
for library,counts in sorted(rgroups.iteritems()):
  sys.stderr.write("%s\t%d\t%d\t%d\n"%(library,counts[0],counts[1],counts[2]))
