#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *03.05.2013
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

from MergeTrimReads import overlap_reads

min_length = 5

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
  new = pysam.AlignedRead()
  if not read1.qname.startswith("M_"): new.qname = "M_"+read1.qname
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
          sys.stderr.write("Do not know how to combine %s BAM tags. Information of one of the reads will get lost during merging.\n"%tag)
      elif tag in htags and value == htags[tag]: # SAME TAG AND VALUE
        new_tags.append((tag,value))
      else: # NEW TAG
        new_tags.append((tag,value))
  for tag,value in htags.iteritems():
    if tag not in stags: new_tags.append((tag,value))
  new.tags = new_tags
  return new

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write it to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'MergeTrim').",default="MergeTrim")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
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
count_overlap = 0
count_nothing = 0

for filename in files:
  if filename == None:
    infile = pysam.Samfile( "-", 'rb' )
  else:
    infile = pysam.Samfile( filename, 'rb' )

  cheader = infile.header
  if ('HD' not in cheader) or ('SO' not in cheader['HD']) or (cheader['HD']['SO'] != 'queryname'):
    sys.stderr.write("Input BAM has wrong sorting order. Skipping input...\n")
    continue

  if outfile == None:
    if options.verbose: sys.stderr.write("Creating output files/streams...\n")

    if options.pipe:
      outfile = pysam.Samfile( "-", fileflags, header = cheader)
      if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
    else:
      outfilename = options.outdir+options.outprefix+".bam"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfile = pysam.Samfile( outfilename, fileflags, header = cheader)

  bread = None
  for read in infile:
    if read.is_paired and bread == None:
      bread = read
      continue
    elif read.is_paired and bread != None:
      if read.qname == bread.qname:
        count_all += 1
        # SAVE SEQUENCES IN MERGING VARIABLES
        if bread.is_read1:
          read1 = bread.seq
          qual1 = bread.qual
          read2 = read.seq
          qual2 = read.qual
        else:
          read1 = read.seq
          qual1 = read.qual
          read2 = bread.seq
          qual2 = bread.qual
        # IF WE HAVE NO QUALITY SCORES, ASSUME EQUAL QUALITY SCORES OF 15
        if qual1 == "*": qual1 = len(read1)*"0"
        if qual2 == "*": qual2 = len(read2)*"0"

        flag,newseq,newqual = overlap_reads(read1,qual1,read2,qual2)

        if newseq != '':
          count_overlap += 1
          new = create_new_read(newseq,newqual,bread,read)
          outfile.write(new)
        else:
          count_nothing += 1
          outfile.write(bread)
          outfile.write(read)
        bread = None
    else: # Single End read -- do nothing
      count_all += 1
      count_nothing += 1
      outfile.write(read)

    if options.verbose and (count_all % 10000 == 0) and (count_all > 0):
        sys.stderr.write("<==> Total %d; Overlapped %d; Kept PE/SR %d\n"%(count_all,count_overlap,count_nothing))

  if not options.pipe:
    infile.close()

sys.stderr.write("Total %d; Overlapped %d; Kept PE/SR %d\n"%(count_all,count_overlap,count_nothing))
if outfile != None:
  outfile.close()
