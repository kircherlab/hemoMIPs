#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Report QC-pass reads to BAM

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *26.11.2011
:Type: tool
:Input: BAM
:Output: BAM

"""

import sys,os
import pysam
import gzip

from optparse import OptionParser,OptionGroup
parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from PIPE and write BAM to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'Filter').",default="Filter")
parser.add_option("-f","--filter",dest="filterlst",help="Filter list of read ids",default=None)
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on (default on)",default=True,action="store_false")
(options, args) = parser.parse_args()

if options.outdir != None and not os.path.isdir(options.outdir):
  sys.stderr.write("Output folder does not exist!\n")
  sys.exit()
elif options.outdir == None:
  options.outdir = ""
else:
  options.outdir = options.outdir.rstrip('/')+'/'

filterset = None
if options.filterlst != None and os.path.exists(options.filterlst):
  filterset = set()
  if options.filterlst.endswith(".gz") or options.filterlst.endswith(".gzip"):
    infile = gzip.open(options.filterlst)
  else:
    infile = open(options.filterlst)
  for line in infile:
    filterset.add(line.strip())
  infile.close()
  if options.verbose:
    sys.stderr.write('Loaded %d read IDs for filtering...\n'%(len(filterset)))

## CREATE OUTPUT FILE(S)/STREAM
fileflags = 'wb'

files = args
if options.pipe: files = [None]
outfile = None

for filename in files:
  if filename == None:
    infile = pysam.Samfile( "-", 'rb' )
  else:
    infile = pysam.Samfile( filename, 'rb' )

  if outfile == None:
    if options.verbose: sys.stderr.write("Creating output files/streams...\n")

    if options.pipe:
      outfile = pysam.Samfile( "-", fileflags, template = infile)
      if options.verbose: sys.stderr.write("BAM output on stdout...\n")
    else:
      outfilename = options.outdir+options.outprefix+".bam"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfile = pysam.Samfile( outfilename , fileflags, template = infile)

    for read in infile:
      if (read.qname in filterset):
        outfile.write(read)
        