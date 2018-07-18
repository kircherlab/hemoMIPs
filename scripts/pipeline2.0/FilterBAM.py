#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *26.09.2011
:Type: tool
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

def is_complex_comp(seq,cutoff):
  counts = []
  for x in 'ACGT':
    counts.append(seq.count(x))
  total = float(sum(counts))
  return (max(counts) <= cutoff*total)

def is_complex_entropy(seq,cutoff):
  counts = []
  for x in 'ACGT':
    counts.append(seq.count(x))
  total = float(sum(counts))
  entropy = 0
  for elem in counts:
    if (total > 0) and (elem/total <> 0):
      entropy -= elem/total*math.log(elem/total,2)
  return (entropy(seq) >= cutoff)

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

def remove_quality_flag(tags,qtype):
  qc_tag_helper = []
  foundzq = False
  if tags != None:
    for tag,value in tags:
      if tag == "ZQ": 
        qc_tag = list(set(list(value))-set([qtype]))
        qc_tag.sort()
        qc_tag = "".join(qc_tag)
        if len(qc_tag) > 0:
          qc_tag_helper.append((tag,qc_tag))
          foundzq = True
      else: qc_tag_helper.append((tag,value))
  return foundzq,qc_tag_helper


from optparse import OptionParser,OptionGroup

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'Filter').",default="Filter")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")

group = OptionGroup(parser, "Complexity filter options")
group.add_option("-e","--entropy",dest="entropy",help="Apply sequence entropy filter",default=False,action="store_true")
group.add_option("-f","--frequency",dest="frequency",help="Apply base frequency filter",default=False,action="store_true")
group.add_option("", "--comp_cutoff", dest="comp_cutoff", help="Entropy value [0.0-2.0] or fraction [0.0-1.0] of most frequent base accepted (default 0.85)",default=0.85,type="float")
parser.add_option_group(group)

group = OptionGroup(parser, "Quality filter options")
group.add_option("-q", "--quality", dest="qual_fixed", help="Apply quality cutoff as fixed value",default=False,action="store_true")
group.add_option("-t", "--trim", dest="qual_trim", help="Trim back sequences instead of filtering",default=False,action="store_true")
group.add_option("-a", "--average", dest="qual_average", help="Apply quality cutoff as average value (ignores qual_number)",default=False,action="store_true")
group.add_option("", "--qual_cutoff", dest="qual_cutoff", help="Quality score cutoff (default 15)",default=15,type="int")
group.add_option("", "--qual_number", dest="qual_number", help="Maximum number of bases below cutoff (default = 0)",type='int',default=0)
parser.add_option_group(group)

group = OptionGroup(parser, "Length filter options")
group.add_option("-l", "--min_length", dest="min_length_cutoff", help="Minimum length cutoff (default 0)",type='int',default=0)
group.add_option("-m", "--max_length", dest="max_length_cutoff", help="Maximum length cutoff (default -1=OFF)",type='int',default=-1)
parser.add_option_group(group)

group = OptionGroup(parser, "Clip bases of reads (applied prior to other filters, considers reverse flag, discards available alignments)")
group.add_option("", "--clip", dest="clip", help="Clip number of bases (default 0, negative values = number of bases to trim from read end).",type='int',default=0)
group.add_option("", "--position", dest="position", help="First base position to remove (1-based)",type='int',default=0)
parser.add_option_group(group)

(options, args) = parser.parse_args()



if options.entropy and options.frequency:
  sys.stderr.write("Only one type of complexity filter allowed!\n")
  sys.exit()
elif options.frequency and (options.comp_cutoff < 0 or options.comp_cutoff > 1):
  sys.stderr.write("Cutoff value out of range [0.0-1.0] for base frequency filter: %.2f!\n"%options.comp_cutoff)
  sys.exit()
elif options.entropy and (options.comp_cutoff < 0 or options.comp_cutoff > 2):
  sys.stderr.write("Cutoff value out of range [0.0-2.0] for entropy filter: %.2f!\n"%options.comp_cutoff)
  sys.exit()

if (options.qual_fixed and options.qual_trim) or (options.qual_fixed and options.qual_average) or (options.qual_trim and options.qual_average):
  sys.stderr.write("Only one type of quality filter allowed!\n")
  sys.exit()

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

    if options.pipe:
      outfile = pysam.Samfile( "-", fileflags, template = infile)
      if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
    else:
      outfilename = options.outdir+options.outprefix+".bam"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfile = pysam.Samfile( outfilename , fileflags, template = infile)

    for read in infile:
      ### READ CLIPPING
      if options.clip > 0 and options.position > 0:
        #sys.stderr.write("CL1...")
        seq = read.seq
        qual = read.qual
        if read.is_reverse:
          read.seq = seq[::-1].translate(table)
          read.qual = qual[::-1]
          read.is_reverse = False

        read.seq = seq[:options.position-1]+seq[options.position-1+options.clip:]
        read.qual = qual[:options.position-1]+qual[options.position-1+options.clip:]
        # REMOVE MAPPING RESULTS -- UGLY!
        read.cigar = None
        read.mapq = 0
        read.tid = 0
        read.mrnm = 0
        read.is_unmapped = True
        if read.is_paired:
          read.mate_is_unmapped = True
          read.is_proper_pair = False
        read.pos = -1
        read.mpos = -1
        new_tags = []
        for tag,value in read.tags:
          if tag != "NM" and tag != "MD": new_tags.append((tag,value))
        read.tags = new_tags
      elif options.clip < 0 and options.position == 0:
        #sys.stderr.write("CL2...")
        qual = read.qual
        seq = read.seq
        if read.is_reverse:
          read.seq = seq[::-1].translate(table)
          read.qual = qual[::-1]
          read.is_reverse = False

        if len(seq) > abs(options.clip):
          read.seq = seq[:options.clip]
          read.qual = qual[:options.clip]
        # REMOVE MAPPING RESULTS -- UGLY!
        read.cigar = None
        read.mapq = 0
        read.tid = 0
        read.mrnm = 0
        read.is_unmapped = True
        if read.is_paired:
          read.mate_is_unmapped = True
          read.is_proper_pair = False
        read.pos = -1
        read.mpos = -1
        new_tags = []
        for tag,value in read.tags:
          if tag != "NM" and tag != "MD": new_tags.append((tag,value))
        read.tags = new_tags

      ### QUALITY FILTER CLIPPING
      if options.qual_trim:
        #sys.stderr.write("QT...")
        seq = read.seq
        qual = read.qual
        if read.is_reverse:
          qual = qual[::-1]
          seq = seq[::-1]
        tol = options.qual_number
        for pos,qs in enumerate(map(lambda x:ord(x)-33,qual)):
          if qs < options.qual_cutoff: tol -= 1
          if tol < 0:
            seq = seq[:pos]
            qual = qual[:pos]
            break
        if read.is_reverse:
          read.seq = seq[::-1]
          read.qual = qual[::-1]
        else:
          read.seq = seq
          read.qual = qual
        # REMOVE MAPPING RESULTS -- UGLY!
        read.cigar = None
        read.mapq = 0
        read.tid = 0
        read.mrnm = 0
        read.is_unmapped = True
        if read.is_paired:
          read.mate_is_unmapped = True
          read.is_proper_pair = False
        read.pos = -1
        read.mpos = -1
        new_tags = []
        for tag,value in read.tags:
          if tag != "NM" and tag != "MD": new_tags.append((tag,value))
        read.tags = new_tags
        read.is_qcfail,read.tags = remove_quality_flag(read.tags,'Q')
      ### ACTUAL QUALITY FILTERS
      elif options.qual_average or options.qual_fixed:
        #sys.stderr.write("QA...")
        if (options.qual_average and (sum(map(lambda x:ord(x)-33,read.qual))/float(len(read.qual)) < options.qual_cutoff)) or (options.qual_fixed and (len(filter(lambda x:ord(x)-33 < options.qual_cutoff,read.qual)) > options.qual_number)):
          read.is_qcfail = True
          read.tags = add_quality_flag(read.tags,'Q')
        else:
          read.is_qcfail,read.tags = remove_quality_flag(read.tags,'Q')

      ### COMPLEXITY FILTERS
      if options.entropy or options.frequency:
        #sys.stderr.write("EF...")
        if (options.entropy and not is_complex_comp(read.seq,options.comp_cutoff)) or (options.frequency and not is_complex_entropy(read.seq,options.comp_cutoff)):
          read.is_qcfail = True
          read.tags = add_quality_flag(read.tags,'C')
        else:
          read.is_qcfail,read.tags = remove_quality_flag(read.tags,'C')

      ### READ LENGTHS FILTERS
      if (options.max_length_cutoff > -1) or (options.min_length_cutoff > 0):
        #sys.stderr.write("LF...")
        if (len(read.seq) < options.min_length_cutoff) or ((options.max_length_cutoff > -1) and (len(read.seq) > options.max_length_cutoff)):
          read.is_qcfail = True
          read.tags = add_quality_flag(read.tags,'L')
        else:
          read.is_qcfail,read.tags = remove_quality_flag(read.tags,'L')

      #sys.stderr.write("Writing.\n")
      outfile.write(read)
