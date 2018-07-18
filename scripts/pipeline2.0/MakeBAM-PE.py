#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Convert SR reads in BAM file to PE (e.g. for using the data with Picard EstimateLibraryComplexity.jar

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *11.01.2016
:Input: BAM
:Output: BAM

"""

import sys,os
import math
import pysam
import string
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv') # COMPLEMENT DNA

#sys.path.append("/mnt/solexa/bin/")
#from library import is_complex_comp
#from library import is_complex_entropy

from optparse import OptionParser

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'Filter').",default="Filter")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
parser.add_option("-m", "--max_length", dest="max_length_cutoff", help="Maximum length of reads (default 100)",type='int',default=100)
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
(options, args) = parser.parse_args()

if options.outprefix == "":
  sys.stderr.write("Outprefix can not be empty!\n")
  sys.exit()

if options.outdir != None and not os.path.isdir(options.outdir):
  sys.stderr.write("Output folder does not exist!\n")
  sys.exit()
elif options.outdir == None:
  options.outdir = ""
else:
  options.outdir = options.outdir.rstrip('/')+'/'

## CREATE OUTPUT FILE(S)/STREAM
fileflags = 'wb'
if options.SAM: fileflags = 'w'

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

    helper = infile.header
    helper['SQ'] = [{'LN': 0, 'SN': '*'}]

    if options.pipe:
      outfile = pysam.Samfile( "-", fileflags, header=helper)
      if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
    else:
      outfilename = options.outdir+options.outprefix+".bam"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfile = pysam.Samfile( outfilename , fileflags, header = helper)
    
    for read in infile:
      if read.is_paired:
        nread = pysam.AlignedRead()
        nread.qname=read.qname
        seq = read.seq
        qual = read.qual
        nread.seq=seq[:options.max_length_cutoff]
        nread.qual=qual[:options.max_length_cutoff]
        nread.is_reverse = read.is_reverse
        nread.is_read1 = read.is_read1
        nread.is_read2 = read.is_read2
        nread.is_unmapped = True
        nread.mate_is_unmapped = True
        nread.pos = -1
        nread.mpos = -1
        nread.is_paired = True
        tags=[]
        for (key,value) in read.tags:
          if key in ["RG","XI","XJ","YI","YJ"]: tags.append((key,value))
        nread.tags=tags
        outfile.write(nread)
      else:
        seq = read.seq
        qual = read.qual
        seqid = read.qname
        if seqid.startswith("M_"): seqid = seqid[2:]
        tags=[]
        for (key,value) in read.tags:
          if key in ["RG","XI","XJ","YI","YJ"]: tags.append((key,value))

        forward = pysam.AlignedRead()
        forward.qname = seqid
        forward.seq = seq[:options.max_length_cutoff]
        forward.qual = qual[:options.max_length_cutoff]
        forward.is_reverse = read.is_reverse
        forward.is_read1 = True
        forward.is_paired = True
        forward.is_unmapped = True
        forward.mate_is_unmapped = True
        forward.pos = -1
        forward.mpos = -1
        forward.tags=tags
        reverse = pysam.AlignedRead()
        reverse.qname = seqid
        reverse.seq = seq[-options.max_length_cutoff:]
        reverse.qual = qual[-options.max_length_cutoff:]
        reverse.is_reverse = not read.is_reverse
        reverse.is_read2 = True
        reverse.is_paired = True
        reverse.is_unmapped = True
        reverse.mate_is_unmapped = True
        reverse.pos = -1
        reverse.mpos = -1
        reverse.tags=tags
        outfile.write(forward)
        outfile.write(reverse)
