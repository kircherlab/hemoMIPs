#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *10.07.2016
"""

import sys, os
from optparse import OptionParser
import pysam
import string

table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv') # COMPLEMENT DNA

parser = OptionParser("%prog [options] input.bam")
parser.add_option("-p","--PIPE", dest="pipe", help="Read BAM from stdin instead of file",default=False,action="store_true")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'Unaligned').",default="Unaligned")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.",default = None)
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
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
  if options.verbose:
    sys.stderr.write("Reading from STDIN...\n")
  filename = "-"
elif len(args) > 0 and os.path.exists(args[0]):
  filename = args[0]
  if options.verbose:
    sys.stderr.write("Reading %s...\n"%filename)
else:
  sys.stderr.write('Error: No valid input file defined\n')
  sys.exit()

cbamfile = pysam.Samfile(filename, "rb")
outfile = None
if options.pipe:
  outfile = pysam.Samfile( "-", "wb", template = cbamfile)
  if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
else:
  outfilename = options.outdir+options.outprefix+".bam"
  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
  outfile = pysam.Samfile(outfilename, "wb", template=cbamfile )

goodFields = set(['RG','XI','XJ','YI','YJ'])
for read in cbamfile:
  new_tags = []
  for (key,value) in read.tags:
    if key in goodFields: new_tags.append((key,value))

  newRead = pysam.AlignedRead()
  newRead.qname = read.qname
  newRead.is_unmapped = True
  newRead.pos = -1
  newRead.mpos = -1
  newRead.tags = new_tags
  newRead.tid = -1
  newRead.rnext = -1
  
  if read.is_reverse:
    newRead.seq = (read.seq[::-1]).translate(table)
    newRead.qual = read.qual[::-1]
  else:
    newRead.seq = read.seq
    newRead.qual = read.qual

  if read.is_paired: newRead.is_paired = True
  if read.is_read1: newRead.is_read1 = True
  if read.is_read2: newRead.is_read2 = True

  outfile.write(newRead)
cbamfile.close()
