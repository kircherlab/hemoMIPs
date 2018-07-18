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

# ADD FLAG TO QUALITY FILTER REAG TAG
def add_quality_flag(tags,qtype):
  qc_tag_helper = []
  foundzq = False
  if tags != None:
    for tag,value in tags:
      if tag == "ZQ": 
        foundzq = True
        qc_tag = list(set(list(value+qtype)))
        qc_tag.sort()
        qc_tag_helper.append((tag,"".join(qc_tag)))
      else: qc_tag_helper.append((tag,value))
  if not foundzq: qc_tag_helper.append(("ZQ",qtype))
  return qc_tag_helper

# CREATE NEW BAM READ FROM ORIGINAL READS AND MERGED SEQUENCE/QUALITY
def create_new_read(new_seq,new_qual,read1,read2):
  global notified
  new = pysam.AlignedRead()
  if not read1.qname.startswith("M_"): new.qname = "C_"+read1.qname
  elif not read1.qname.startswith("M_"): new.qname = "CM_"+read1.qname

  else: new.qname = read1.qname
  new.seq = new_seq
  new.qual = new_qual
  new.is_unmapped = True
  new.pos = -1
  new.mpos = -1

  new.is_qcfail = read1.is_qcfail and read2.is_qcfail
  if read1.tags != None: htags = dict(read1.tags)
  else: htags = {}
  if (len(new_seq) < min_length):
    new.is_qcfail = True
    if "ZQ" in htags: htags["ZQ"]+="L"
    else: htags["ZQ"]="L"
  stags = set()
  new_tags = []
  if read2.tags != None:
    for tag,value in read2.tags:
      stags.add(tag)
      if tag == "NM" or tag == "MD": continue
      elif tag in htags and value != htags[tag]: # NEW TAG DIFF VALUE
        if tag == "ZQ":
          qc_tag = list(set(list(value+htags[tag])))
          qc_tag.sort()
          new_tags.append((tag,"".join(qc_tag)))
        else:
          if tag not in notified:
            sys.stderr.write("Do not know how to combine %s BAM tags. Information of one of the reads will get lost during merging.\n"%tag)
            notified.add(tag)
      elif tag in htags and value == htags[tag]: # SAME TAG AND VALUE
        new_tags.append((tag,value))
      else: # NEW TAG
        new_tags.append((tag,value))
  for tag,value in htags.iteritems():
    if tag not in stags: new_tags.append((tag,value))
  new.tags = new_tags
  return new


options_adapter_F="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG"

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write it to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'MergeTrim').",default="MergeTrim")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
parser.add_option("-f", "--adapter", dest="adapter_F", help="Adapter that is observed for the additional read (def. Multiplex: %s)"%options_adapter_F[:maxadapter_comp],default=options_adapter_F)
parser.add_option("-t","--trimCutoff",dest="trimCutoff",help="Lowest number of adapter bases to be observed for Single Read trimming (default 1)",default=1,type="int")
parser.add_option("--readInTag", dest="readInTag", help="Generate consensus using the read in the following BAM tag (def XJ)",default="XJ")
parser.add_option("--qualInTag", dest="qualInTag", help="Generate consensus using the qualities in the following BAM tag (def YJ)",default="YJ")
parser.add_option("--reverse", dest="reverse", help="Additional read is reverse complement (def Off)",default=False,action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
(options, args) = parser.parse_args()

set_adapter_sequences(options.adapter_F, "", options.adapter_F, maxadapter_comp)
set_options(options.trimCutoff, False, False, False)
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
count_consensus_fix = 0

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
    tagSeq = ""
    tagQual = "*"
    
    for key,value in read.tags:
      if key == options.readInTag:
        tagSeq = value
      elif key == options.qualInTag:
        tagQual = value[::-1]

    if len(tagSeq) == 0: 
      count_nothing += 1
      outfile.write(read)
      continue
    if tagQual == "*": tagQual = len(tagSeq)*"0"
    
    read1 = read.seq
    qual1 = read.qual
    if qual1 == "*": qual1 = len(read1)*"0"
    
    consseq,consqual = "","*"
    flag,newseq,newqual = process_SR(tagSeq,tagQual)
    #sys.stderr.write("%s %s %s | %s %s\n"%(flag,newseq,newqual,tagSeq,tagQual))
    if flag == 'D': 
      count_failed += 1
      read.tags = add_quality_flag(read.tags,flag)
    else:
      if newseq != '':
        if len(newseq) == len(read1):
          if options.reverse:
            flag,consseq,consqual = consensus_reads(read1,qual1,revcompl(newseq),newqual[::-1])
          else:
            flag,consseq,consqual = consensus_reads(read1,qual1,newseq,newqual)
          if len(consseq) != 0: count_consensus += 1
        else: 
          # Try fixed value triming and consensus first
          if options.reverse:
            flag,consseq,consqual = consensus_reads(read1,qual1,revcompl(tagSeq[:len(read1)]),tagQual[:len(read1)][::-1])
          else:
            flag,consseq,consqual = consensus_reads(read1,qual1,tagSeq[:len(read1)],tagQual[:len(read1)])
          # Try merging with different length
          if len(consseq) == 0: 
            if options.reverse:
              flag,consseq,consqual = overlap_reads(read1,qual1,newseq,newqual)
            else:
              flag,consseq,consqual = overlap_reads(read1,qual1,revcompl(newseq),newqual[::-1])
            if len(consseq) != 0: count_consensus += 1
          else:
            if len(consseq) != 0: count_consensus_fix += 1
      elif len(tagSeq) >= len(read1) :
        if options.reverse:
          flag,consseq,consqual = consensus_reads(read1,qual1,revcompl(tagSeq[:len(read1)]),tagQual[:len(read1)][::-1])
        else:
          flag,consseq,consqual = consensus_reads(read1,qual1,tagSeq[:len(read1)],tagQual[:len(read1)])
        if len(consseq) != 0: count_consensus_fix += 1
      
    if consseq != '': 
      read.seq = consseq
      read.qual = consqual
      read.is_unmapped = True
      read.pos = -1
      read.mpos = -1
      read.cigar = None
      read.mapq = 0
      read.tid = 0
      read.mrnm = 0
      
      new_tags = []
      for tag,value in read.tags:
        if tag != options.readInTag and tag != options.qualInTag and tag != "NM" and tag != "MD": new_tags.append((tag,value))
      read.tags = new_tags
      outfile.write(read)
    else:
      if flag == '': count_failed += 1
      outfile.write(read)

    if options.verbose and (count_all % 10000 == 0) and (count_all > 0):
        sys.stderr.write("<==> Total %d; Consensus %d; ConsensusFix %d; Kept %d; Adapter dimers/failed %d\n"%(count_all,count_consensus,count_consensus_fix,count_nothing,count_failed))

  if not options.pipe:
    infile.close()

sys.stderr.write("Total %d; Consensus %d; ConsensusFix %d; Kept %d; Adapter dimers/failed %d\n"%(count_all,count_consensus,count_consensus_fix,count_nothing,count_failed))
if outfile != None:
  outfile.close()
