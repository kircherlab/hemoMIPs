#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
Redefine BAM read group based on index sequences and quality scores available
from the BAM file by comparing them against a list of known indexes. If
requested seperate in multiple BAM files.

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *04.01.2012
:Type: tool
:Input: BAM
:Output: BAM

"""

import sys
import os
import string
table = string.maketrans('.','N')

import pysam
from collections import defaultdict
from library import read_fastq
from optparse import OptionParser,OptionGroup

index_length1 = None
index_length2 = None
ireadlength1 = None
ireadlength2 = None

next_bases1 = 'ATCTCGTATGCCGTCTTCTGCTTG'
next_bases2 = 'GTGTAGATCTCGGTGGTCGCCGTATCATT'
nucleotides = ["A","C","G","T","N"]

INF_QUALITY = 200

def get_min_qual(qualstring,qualN):
  cmin = INF_QUALITY
  allowZero = 1
  for pos,i in enumerate(map(lambda x:ord(x)-33,qualstring)):
    if (len(qualN) > pos) and qualN[pos] == "N": continue
    if allowZero > 0 and i == 0: allowZero -= 1
    elif i < cmin: cmin = i
  return cmin

def setN(cindex,setNstring):
  if (cindex == None) or (setNstring == None) or len(setNstring) != len(cindex):
    return cindex
  else:
    cur = list(cindex)
    for cpos,base in enumerate(setNstring.upper()):
      if base == "N": cur[cpos]="N"
    return "".join(cur)

def extreme_read_index_file(filename,flength1=None,flength2=None):
  global next_bases1
  global next_bases2
  global nucleotides
  mlindex1 = None
  mlindex2 = None
  count1 = 0
  count2 = 0
  isDouble = False

  # FIRST SCAN OF INDEX FILE
  infile = open(filename)
  for line in infile:
    if len(line.strip()) > 0: # CHECK FOR INDEX IN FIRST COLUMN
      if line.startswith('#'): continue
      fields = line.split()
      if len(fields) >= 2:
        cindex = fields[0].upper()
        # SKIP NON VALID LINES
        valid = True
        for elem in cindex:
          if elem not in nucleotides:
            valid=False
            break
        if not valid: continue
        if (mlindex1 == None) or (len(cindex) > mlindex1): mlindex1 = len(cindex)
        count1 += 1
      if len(fields) >= 3: # CHECK FOR DOUBLE INDEX FILE
        cindex = fields[1].upper()
        # SKIP NON VALID LINES
        valid = True
        for elem in cindex:
          if elem not in nucleotides:
            valid=False
            break
        if (not valid) and (mlindex2 == None): continue
        elif (not valid):
          sys.stderr.write("Unexpected line in index file, assuming non-double index: %s\n"%line.strip()[:40])
          sys.exit()
        if (mlindex2 == None) or (len(cindex) > mlindex2): mlindex2 = len(cindex)
        count2 += 1
        isDouble = True
  infile.close()

  if (count1 != count2) and (mlindex2 != 0):
    sys.stderr.write("Inconsistent double index input file, assuming non-double index.\n")
    mlindex2 = 0
    isDouble = False

  if (flength1 > 0) and (flength1 > mlindex1): mlindex1 = flength1
  if (flength2 > 0) and (flength2 > mlindex2): mlindex2 = flength2

  if not isDouble:
    sys.stderr.write("Extreme identification requires double index data!")
    sys.exit()

  pairs = {}
  index_ones = defaultdict(int)
  infile = open(filename) # READ IN INDEX SEQUENCES
  for line in infile:
    if line.startswith('#'): continue
    variant_helper =[]
    if len(line.strip()) > 0:
      fields = line.split()
      if len(fields) > 2:
        cindex1 = fields[0].upper()
        cindex1 = cindex1+next_bases1[:(mlindex1-len(cindex1))]
        cindex2 = fields[1].upper()
        cindex2 = cindex2+next_bases2[:(mlindex2-len(cindex2))]
        samplename = fields[2]
        pairs[samplename] = (cindex1,cindex2)
        index_ones[cindex1] += 1
  infile.close()

  names = {}
  to_fix = {}
  for samplename,(cindex1,cindex2) in pairs.iteritems():
    if index_ones[cindex1] == 1: names[cindex1] = (samplename,None)
    else:
      if cindex1 in to_fix: to_fix[cindex1].append([samplename,cindex2])
      else: to_fix[cindex1]= [[samplename,cindex2]]

  for cindex1,sample_indexes in to_fix.iteritems():
    masked_indexes = list(sample_indexes)
    for ind in range(len(sample_indexes[0][1])):
      dbases = set()
      for name,cindex in sample_indexes: dbases.add(cindex[ind])
      if len(dbases) != len(sample_indexes):
        for sInd in range(len(sample_indexes)): masked_indexes[sInd][1]=masked_indexes[sInd][1][:ind]+"N"+masked_indexes[sInd][1][ind+1:]
    names[cindex1] = (None,masked_indexes)

  return mlindex1,mlindex2,names,None,None

def read_index_file(filename,skipFirst=True,oneDist=True,nbases=True,flength1=None,flength2=None):
  global next_bases1
  global next_bases2
  global nucleotides
  global options
  mlindex1 = None
  mlindex2 = None
  count1 = 0
  count2 = 0
  isDouble = False

  # FIRST SCAN OF INDEX FILE
  infile = open(filename)
  for line in infile:
    if len(line.strip()) > 0: # CHECK FOR INDEX IN FIRST COLUMN
      if line.startswith('#'): continue
      fields = line.split()
      if len(fields) >= 2:
        cindex = fields[0].upper()
        # SKIP NON VALID LINES
        valid = True
        for elem in cindex:
          if elem not in nucleotides:
            valid=False
            break
        if not valid: continue
        if (mlindex1 == None) or (len(cindex) > mlindex1): mlindex1 = len(cindex)
        count1 += 1
      if len(fields) >= 3: # CHECK FOR DOUBLE INDEX FILE
        cindex = fields[1].upper()
        # SKIP NON VALID LINES
        valid = True
        for elem in cindex:
          if elem not in nucleotides:
            valid=False
            break
        if (not valid) and (mlindex2 == None): continue
        elif (not valid):
          sys.stderr.write("Unexpected line in index file, assuming non-double index: %s\n"%line.strip()[:40])
          sys.exit()
        if (mlindex2 == None) or (len(cindex) > mlindex2): mlindex2 = len(cindex)
        count2 += 1
        isDouble = True
  infile.close()
  if (count1 != count2) and (mlindex2 != 0):
    sys.stderr.write("Inconsistent double index input file, assuming non-double index.\n")
    mlindex2 = 0
    isDouble = False

  if (flength1 > 0) and (flength1 > mlindex1): mlindex1 = flength1
  if (flength2 > 0) and (flength2 > mlindex2): mlindex2 = flength2

  res1 = {}
  res2 = {}
  names = {}
  names[None]="conflict"
  infile = open(filename) # READ IN INDEX SEQUENCES
  for line in infile:
    if line.startswith('#'): continue
    variant_helper =[]
    if len(line.strip()) > 0:
      fields = line.split()
      if len(fields) >= 2:
        cindex1 = fields[0].upper()
        # SKIP NON VALID LINES
        valid = True
        for elem in cindex1:
          if elem not in nucleotides:
            valid=False
            break
        if not valid: continue

        # REMOVE SOME BASES OF FIRST INDEX AND REPLACE BY N
        if (options.setN != None) and (len(options.setN) == mlindex1):
          cur = list(cindex1+next_bases1[:(mlindex1-len(cindex1))])
          nvariants = []
          for cpos,base in enumerate(options.setN.upper()):
            if base == "N":
              cur[cpos]="N"
          cindex1 = "".join(cur)

        if not options.skip_error1: variant_helper.append((cindex1,res1,mlindex1,next_bases1))
        hcindex = cindex1+next_bases1[:(mlindex1-len(cindex1))]
        if (hcindex in res1) and (res1[hcindex]!=cindex1):
          sys.stderr.write("Error: Index (1) causes a conflict %s with variant of %s\n"%(hcindex,res1[hcindex]))
          res1[hcindex]=None
        else: res1[hcindex]=cindex1

        cindex2 = None
        if isDouble:
          cindex2 = fields[1].upper()

          # REMOVE SOME BASES OF SECOND INDEX AND REPLACE BY N
          if (options.setN_2nd != None) and (len(options.setN_2nd) == mlindex2):
            cur = list(cindex2+next_bases2[:(mlindex2-len(cindex2))])
            for cpos,base in enumerate(options.setN_2nd.upper()):
              if base == "N":
                cur[cpos]="N"
            cindex2 = "".join(cur)

          if not options.skip_error2: variant_helper.append((cindex2,res2,mlindex2,next_bases2))
          hcindex = cindex2+next_bases2[:(mlindex2-len(cindex2))]
          if (hcindex in res2) and (res2[hcindex]!=cindex2):
            sys.stderr.write("Error: Index (2) causes a conflict %s with variant of %s\n"%(hcindex,res2[hcindex]))
            res2[hcindex]=None
          else: res2[hcindex]=cindex2

          names[(cindex1,cindex2)]=fields[2] #.split("-")[-1]
        else:
          names[(cindex1,cindex2)]=fields[1] #.split("-")[-1]

        for cindex,res,mlindex,next_bases in variant_helper: # ERROR CORRECTION
          # SKIP FIRST BASE
          if skipFirst:
            cur = cindex[1:]+next_bases[:(mlindex+1-len(cindex))]
            if (cur in res) and (res[cur]!=cindex) and (res[cur]!=None):
              sys.stderr.write("Error: Skipping first base in index causes a conflict %s->%s with variant of %s\n"%(cindex,cur,res[cur]))
              res[cur]=None
            elif (cur not in res): res[cur]=cindex
          # CREATE 1nt MUTANTS
          if oneDist:
            for cpos in range(mlindex):
              for base in nucleotides:
                cur = list(cindex+next_bases[:(mlindex-len(cindex))])
                cur[cpos]=base
                cur="".join(cur)
                if (cur in res) and (res[cur]!=cindex) and (res[cur]!=None):
                  sys.stderr.write("Error: Creating mutants of indexes causes a conflict %s->%s with variant of %s\n"%(cindex,cur,res[cur]))
                  res[cur]=None
                  #sys.exit()
                elif (cur not in res): res[cur]=cindex
          # CREATE Ns
          if nbases:
            for cpos in range(mlindex):
              cur = list(cindex+next_bases[:(mlindex-len(cindex))])
              cur[cpos]="N"
              cur="".join(cur)
              if (cur in res) and (res[cur]!=cindex) and (res[cur]!=None):
                sys.stderr.write("Error: creating mutants of indexes causes a conflict %s->%s with variant of %s\n"%(cindex,cur,res[cur]))
                res[cur]=None
              elif (cur not in res): res[cur]=cindex
      else:
        sys.stderr.write("Error: Line is not valid: %s\n"%line.strip())
        sys.exit()

  return mlindex1,mlindex2,names,res1,res2


parser = OptionParser("%prog [options] seq_files")
parser.add_option("-p","--PIPE",dest="pipe",help="Read from and write to PIPE",default=False,action="store_true")
parser.add_option("", "--SAM", dest="SAM", help="Input/output is SAM not BAM.",default=False,action="store_true")
parser.add_option("-i","--index",dest="index",help="File describing index sequences used")
parser.add_option("--bases_after_index",dest="nextbases1",help="Nucleotides in adapter downstream of index (default %s)"%next_bases1,default=next_bases1)
parser.add_option("--bases_after_2ndindex",dest="nextbases2",help="Nucleotides in adapter downstream of second index (default %s)"%next_bases2,default=next_bases2)
parser.add_option("--nextera",dest="nextera",help="Order of reads like Nextera",default=False, action="store_true")

group = OptionGroup(parser, "Output options")
group.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
group.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'remultiplex').",default="remultiplex")
group.add_option("--separate_files",dest="separate_files",help="Create output file for each index",default=False,action="store_true")
group.add_option("--remove",dest="remove",help="Remove reads not matching any know index",default=False,action="store_true")
group.add_option("--summary",dest="summary",help="Show summary for reads writen to separate files",default=False,action="store_true")
group.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
parser.add_option_group(group)

group = OptionGroup(parser, "Quality filtering of index read")
group.add_option("-q","--quality",dest="quality",help="Minimum quality score of bases in both indexes (default 0:off)",type='int',default=0)
group.add_option("--qualityN",dest="qualityN",help="Consider subset of bases for quality score filter (e.g. NyyyyyN removes first and last base)",default="")
group.add_option("--2nd_qualityN",dest="qualityN_2nd",help="Consider subset of bases for quality score filter (e.g. NyyyyyN removes first and last base)",default="")
group.add_option("--setN",dest="setN",help="Set bases from index read to N (e.g. NyyyyyN removes first and last base)")
group.add_option("--2nd_setN",dest="setN_2nd",help="Set bases from second index read to N (e.g. NyyyyyN removes first and last base)")
group.add_option("--extreme",dest="extreme",help="Consider second index mostly failed. Only check for index informative positions",default=False,action="store_true")
parser.add_option_group(group)

group = OptionGroup(parser, "Sequencing error correction")
group.add_option("--no_error_1st",dest="skip_error1",help="Deactivate error correction for index 1",default=False,action="store_true")
group.add_option("--no_error_2nd",dest="skip_error2",help="Deactivate error correction for index 2",default=False,action="store_true")
group.add_option("--no_skip_first_base",dest="skip_first_base",help="Do not consider skipping first base of index",default=True,action="store_false")
group.add_option("--no_mutants",dest="mutants",help="Do not consider index mutants within one base distance",default=True,action="store_false")
group.add_option("--no_Ns",dest="Nbases",help="Do not allow one N in index sequences",default=True,action="store_false")
parser.add_option_group(group)
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

fileflags = 'b'
if options.SAM: fileflags = ''

next_bases1=options.nextbases1.upper()
next_bases2=options.nextbases2.upper()

options.qualityN = options.qualityN.upper()
options.qualityN_2nd = options.qualityN_2nd.upper()

files = args
if options.pipe: files = [None]

first = True
SAMheader = { 'HD': {'VN': '1.4','SO':'queryname'}, 'SQ': [{'LN': 0, 'SN': '*'}] }
names,indexes1,indexes2={},{},{}
isdoubleIndex = False

outfiles = {}
for filename in files:
  if filename == None:
    infile = pysam.Samfile( "-", 'r'+fileflags )
    header= infile.header
    if "SQ" in header:
      SAMheader['SQ'] = header['SQ']
  else:
    infile = pysam.Samfile( filename, 'r'+fileflags )
    header= infile.header
    if "SQ" in header:
      SAMheader['SQ'] = header['SQ']

  for read in infile:

    oindseq1,oindqual1,oindseq2,oindqual2 = None,None,None,None
    for key,value in read.tags:
      if key == "XI": oindseq1 = value
      elif key == "YI": oindqual1 = value
      elif key == "XJ": oindseq2 = value
      elif key == "YJ": oindqual2 = value

    if first:
      first = False
      ## CREATE INDEX VARIANTS
      if (options.index != None) and os.path.exists(options.index):
        if options.verbose: sys.stderr.write("Calculating index variants...\n")
        if not options.extreme:
          index_length1,index_length2,names,indexes1,indexes2 = read_index_file(options.index,skipFirst=options.skip_first_base,oneDist=options.mutants,nbases=options.Nbases,flength1=None if oindseq1 == None else len(oindseq1),flength2=None if oindseq2 == None else len(oindseq2))
          if options.verbose: sys.stderr.write("Created %d (%d for second index read) putative index variants for %d index sequences.\n"%(len(indexes1),len(indexes2),len(names)-1))
          isdoubleIndex = (len(indexes2) > 0)
        else:
          index_length1,index_length2,names,indexes1,indexes2 = extreme_read_index_file(options.index,flength1=None if oindseq1 == None else len(oindseq1),flength2=None if oindseq2 == None else len(oindseq2))
          isdoubleIndex = True
        if not options.extreme:
          names["unknown"]="unknown"
          if isdoubleIndex: names["wrong"]="wrong"
        else:
          names["unknown"]=("unknown",None)
      else:
        sys.stderr.write("Need valid index file. No file defined or file does not exist.\n")
        sys.exit()

      if ireadlength2 == 0: ireadlength2 = None
      
      rgs=[]
      for v in names.values():
          rgs.append({'ID': v, 'PL':'Illumina', 'LB': v, 'SM': v})
      SAMheader['RG'] = rgs    

      if len(outfiles) == 0:
        if options.verbose: sys.stderr.write("Creating output files/streams...\n")
        if options.pipe:
          if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
          if len(names) > 0:
            ostream = pysam.Samfile( "-", 'w'+fileflags, header = SAMheader)
            for tag,name in names.iteritems():
              outfiles[tag] = [ostream,0,None]
          else:
            outfiles[None] = [pysam.Samfile( "-", 'w'+fileflags, header = SAMheader),0,None]
        else:
          if options.separate_files and len(names) > 0:
            for tag,name in names.iteritems():
              outfilename = options.outdir+options.outprefix+"_%s"%name+".bam"
              if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
              outfiles[tag] = [pysam.Samfile( outfilename , 'w'+fileflags, header = SAMheader),0,outfilename]
          elif len(names) > 0:
            outfilename = options.outdir+options.outprefix+".bam"
            if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
            ostream = pysam.Samfile( outfilename , 'w'+fileflags, header = SAMheader)
            for tag,name in names.iteritems():
              outfiles[tag] = [ostream,0,None]
          else:
            outfilename = options.outdir+options.outprefix+".bam"
            if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
            outfiles[None] = [pysam.Samfile( outfilename , 'w'+fileflags, header = SAMheader),0,None]
   
    if oindseq1 == None:
      if not options.remove:
        outfile[None].write(read)
    else:
      is_qcfail = True

      cind = None
      tags = []
      if not options.extreme:
        indseq1 = setN(oindseq1,options.setN)
        indseq2 = setN(oindseq2,options.setN_2nd)

        minqual1 = INF_QUALITY
        if options.quality > 0: minqual1 = get_min_qual(oindqual1,options.qualityN)
        minqual2 = INF_QUALITY
        if isdoubleIndex and options.quality > 0: minqual2 = get_min_qual(oindqual2,options.qualityN_2nd)

        if ireadlength2 == None: indseq2 = None

        if (len(indexes1) == 0): 
          cind = None
          is_qcfail = False
        elif (not isdoubleIndex): # Single Index
          if (indseq1 in indexes1): cind = (indexes1[indseq1],None)
          else: cind = "unknown"
          if cind == (None,None): cind = None
          if (minqual1 >= options.quality): is_qcfail = False
        else: #Double index
          if (minqual1 >= options.quality) and (minqual2 >= options.quality): is_qcfail = False
          if (indseq1 in indexes1) and (indseq2 in indexes2):
            assign1 = indexes1[indseq1]
            assign2 = indexes2[indseq2]
            if (assign1 == None) or (assign2 == None): cind = None
            elif (assign1,assign2) not in names: cind = "wrong"
            else: cind = (assign1,assign2)
          else: cind = "unknown"
        if cind in names: tags.append(("RG",names[cind]))
      else:
        cind = "unknown"
        minqual1 = INF_QUALITY
        if options.quality > 0: minqual1 = get_min_qual(oindqual1,options.qualityN)
        if (minqual1 >= options.quality): is_qcfail = False

      if cind in names: new_rg = names[cind]
      elif cind != None: new_rg = cind
      else: new_rg = "unknown"

      new_tags = []
      found_rg,found_zq = False,False
      for key,value in read.tags:
        if key == 'RG':
          new_tags.append(('RG',new_rg))
          found_rg = True
        elif key == "ZQ":
          if "I" in value and is_qcfail:
            new_tags.append((key,value))
            found_zq = True
          elif "I" not in value and is_qcfail:
            value = list(value)
            value.append("I")
            value.sort()
            value = "".join(value)
            new_tags.append((key,value))
            read.is_qcfail = True
            found_zq = True
          elif "I" in value and not is_qcfail:
            value = value.replace("I","")
            if len(value) == 0:
              read.is_qcfail = False
          else:
            new_tags.append((key,value))
            found_zq = True
        else:
          new_tags.append((key,value))
      if not found_rg: new_tags.append(('RG',new_rg))
      if not found_zq and is_qcfail: new_tags.append(('ZQ','I'))
      read.tags = new_tags

      if len(outfiles) > 1:
        outfiles[cind][1]+=1
        if (new_rg != "unknown") or (not options.remove): outfiles[cind][0].write(read)
      else:
        outfiles[None][1]+=1
        if (new_rg != "unknown") or (not options.remove): outfiles[None][0].write(read)

  # CREATE SUMMARY AND CLEAN UP NOT USED FILES
  summary = []
  closed = False
  if options.verbose and len(outfiles) > 1: sys.stderr.write("Cleaning up not needed output files...")
  for tag,value in outfiles.iteritems():
    if value[2] != None and not closed: 
      value[0].close()
      close = True
    if value[1] == 0 and value[2] != None:
      try:
        if options.verbose: sys.stderr.write("Removing: %s\n"%(value[2]))
        os.remove(value[2])
      except: sys.stderr.write("Error: Removing output files %s %s\n"(value[2]))
      if options.summary and tag in names: summary.append("  Observed %10d reads for %s (%s)\n"%(value[1],str(names[tag]),str(tag)))
    elif options.summary:
      if (len(outfiles) > 1) and tag in names: summary.append("  Observed %10d reads for %s (%s)\n"%(value[1],str(names[tag]),str(tag)))
      else: summary.append("  Created %10d sequences\n"%value[1])
  summary.sort()
  summary.reverse()
  if options.summary: sys.stderr.write("Summary:\n")
  for line in summary:
    sys.stderr.write(line)

  if filename != None:
    infile.close()
