#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *14.05.2014
:Type: tool
:Input: BAM
:Output: BAM

"""

import sys, os
import math
import random

import pysam
from optparse import OptionParser,OptionGroup
import string

from MergeTrimReads import set_adapter_sequences, set_options, set_keys, process_SR, overlap_reads, consensus_reads
table = string.maketrans('TGCA','ACGT') # COMPLEMENT DNA

maxadapter_comp = 30
min_length = 5
notified = set()

def revcompl(seq):
  global table
  seq=seq.translate(table)
  return seq[::-1]

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write it to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'MergeTrim').",default="MergeTrim")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
parser.add_option("--insertOne",dest="insertOne",help="For reads failing regular consensus calling, check whether they can be rescued by inserting 1 base at the beginning of either index read (default off)",default=False,action="store_true")
parser.add_option("--reverse", dest="reverse", help="Second read out is reverse complement (def Off)",default=False,action="store_true")
parser.add_option("--remove",dest="remove",help="Remove reads failing consensus calling (default off)",default=False,action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
(options, args) = parser.parse_args()

adapter_F="ATCTCGTATGCCGTCTTCTGCTTG"
set_adapter_sequences(adapter_F, "", adapter_F, maxadapter_comp)
set_options(1, False, False, False)
set_keys("")

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

################################
## PROCESS FILES
################################

## CREATE OUTPUT FILE(S)/STREAM
fileflags = 'wb'
if options.SAM: fileflags = 'w'

files = args
if options.pipe: files = [None]

outfile = None
count_all = 0
count_nothing = 0
count_failed = 0
count_consensus = 0

for filename in files:
  if filename == None:
    infile = pysam.Samfile( "-", 'rb' )
  else:
    infile = pysam.Samfile( filename, 'rb' )

  cheader = infile.header

  if outfile == None:
    if options.verbose: sys.stderr.write("Creating output files/streams...\n")

    if options.pipe:
      if len(cheader) == 0:
        cheader['HD'] = {'VN': '1.0'}
      outfile = pysam.Samfile( "-", fileflags, header = cheader)
      if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
    else:
      outfilename = options.outdir+options.outprefix+".bam"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfile = pysam.Samfile( outfilename, fileflags, header = cheader)

  for read in infile:
    count_all += 1

    tagSeq1 = ""
    tagQual1 = "*"
    tagSeq2 = ""
    tagQual2 = "*"

    tmp_readTags = []
    for key,value in read.tags:
      if key == "XI":
        tagSeq1 = value
      elif key == "XJ":
        tagSeq2 = value
      elif key == "YI":
        tagQual1 = value
      elif key == "YJ":
        tagQual2 = value
      else:
        tmp_readTags.append((key,value))
        
    if (len(tagSeq1) == 0) or (len(tagSeq2) == 0): 
      count_nothing += 1
      outfile.write(read)
      continue

    if tagQual1 == "*": tagQual1 = len(tagSeq1)*"0"
    if tagQual2 == "*": tagQual2 = len(tagSeq2)*"0"
    
    consseq,consqual = "","*"
    if len(tagSeq1) == len(tagSeq2):
      if options.reverse:
        flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,revcompl(tagSeq2),tagQual2[::-1])
      else:
        flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,tagSeq2,tagQual2)
      if len(consseq) == 0 and options.insertOne: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads("N"+tagSeq1[:-1],"1"+tagQual1[:-1],revcompl(tagSeq2),tagQual2[::-1])
        else:
          flag,consseq,consqual = consensus_reads("N"+tagSeq1[:-1],"1"+tagQual1[:-1],tagSeq2,tagQual2)
      if len(consseq) == 0 and options.insertOne: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,revcompl("N"+tagSeq2[:-1]),("1"+tagQual2[:-1])[::-1])
        else:
          flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,"N"+tagSeq2[:-1],"1"+tagQual2[:-1])
        
        
    elif len(tagSeq1) > len(tagSeq2) :
      # Try merging with different length
      if options.reverse:
        flag,consseq,consqual = overlap_reads(tagSeq1,tagQual1,tagSeq2,tagQual2)
      else:
        flag,consseq,consqual = overlap_reads(tagSeq1,tagQual1,revcompl(tagSeq2),tagQual2[::-1])
      if len(consseq) == 0: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads(tagSeq1[:len(tagSeq2)],tagQual1[:len(tagSeq2)],revcompl(tagSeq2),tagQual2[::-1])
        else:
          flag,consseq,consqual = consensus_reads(tagSeq1[:len(tagSeq2)],tagQual1[:len(tagSeq2)],tagSeq2,tagQual2)
      if len(consseq) == 0 and options.insertOne: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads("N"+tagSeq1[:(len(tagSeq2)-1)],"1"+tagQual1[:(len(tagSeq2)-1)],revcompl(tagSeq2),tagQual2[::-1])
        else:
          flag,consseq,consqual = consensus_reads("N"+tagSeq1[:(len(tagSeq2)-1)],"1"+tagQual1[:(len(tagSeq2)-1)],tagSeq2,tagQual2)
      if len(consseq) == 0 and options.insertOne: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads(tagSeq1[:len(tagSeq2)],tagQual1[:len(tagSeq2)],revcompl("N"+tagSeq2[:-1]),("1"+tagQual2[:-1])[::-1])
        else:
          flag,consseq,consqual = consensus_reads(tagSeq1[:len(tagSeq2)],tagQual1[:len(tagSeq2)],"N"+tagSeq2[:-1],"1"+tagQual2[:-1])

    else: 
      # Try merging with different length
      if options.reverse:
        flag,consseq,consqual = overlap_reads(tagSeq1,tagQual1,tagSeq2,tagQual2)
      else:
        flag,consseq,consqual = overlap_reads(tagSeq1,tagQual1,revcompl(tagSeq2),tagQual2[::-1])
      if len(consseq) == 0: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,revcompl(tagSeq2[:len(tagSeq1)]),tagQual2[:len(tagSeq1)][::-1])
        else:
          flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,tagSeq2[:len(tagSeq1)],tagQual2[:len(tagSeq1)])
      if len(consseq) == 0 and options.insertOne: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,revcompl("N"+tagSeq2[:(len(tagSeq1)-1)]),("0"+tagQual2[:(len(tagSeq1)-1)])[::-1])
        else:
          flag,consseq,consqual = consensus_reads(tagSeq1,tagQual1,"N"+tagSeq2[:(len(tagSeq1)-1)],"0"+tagQual2[:(len(tagSeq1)-1)])
      if len(consseq) == 0 and options.insertOne: 
        if options.reverse:
          flag,consseq,consqual = consensus_reads("N"+tagSeq1[:-1],"1"+tagQual1[:-1],revcompl(tagSeq2[:len(tagSeq1)]),tagQual2[:len(tagSeq1)][::-1])
        else:
          flag,consseq,consqual = consensus_reads("N"+tagSeq1[:-1],"1"+tagQual1[:-1],tagSeq2[:len(tagSeq1)],tagQual2[:len(tagSeq1)])

    if consseq != '': 
      count_consensus += 1
      tmp_readTags.append(("XI",consseq))
      tmp_readTags.append(("YI",consqual))
      read.tags = tmp_readTags
      outfile.write(read)
    else:
      if flag == '': count_failed += 1
      if not options.remove:
        outfile.write(read)

    if options.verbose and (count_all % 10000 == 0) and (count_all > 0):
        sys.stderr.write("<==> Total %d; Consensus %d; Kept %d; Adapter dimers/failed %d\n"%(count_all,count_consensus,count_nothing,count_failed))

  if not options.pipe:
    infile.close()

sys.stderr.write("Total %d; Consensus %d; Kept %d; Adapter dimers/failed %d\n"%(count_all,count_consensus,count_nothing,count_failed))
if outfile != None:
  outfile.close()
