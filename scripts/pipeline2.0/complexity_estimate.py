#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
Estimates the complexity of a sequencing library. Supports various input formats and
additional filters. Creates an report pdf

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *14.04.2009
:Type: tool
:Input: eland, bowtie, sam
:Output: pdf

options:
        --type			    Input file type (eland,bowtie,sam)
    -e, --entropy                   Entropy cutoff to be used for sequence filtering (default = 0.8)
    -l, --no_length                 Do not consider the length of the alignment as additional information
    -f, --min_length                Do not consider reads shorter than N nt
    -t, --no_mt                     Do not consider hits to MT
    -m, --merged_only               Consider only merged reads from sam/bowtie files
    -c, --coverage_filter           Do not consider very frequent coordinates in estimate (in percentiles: 0-100,default 100=OFF)
        --outfile_coverage_filter   Print complexity filtered sequences to file (default OFF)
    -r, --repeat_mask               Do not consider reads of which at least FRACTION overlaps a repeat masked region (does not work in combination with no_length!)
        --twoBitToFa                Path to twoBitToFa binary (default /mnt/solexa/bin/blat64bit/twoBitToFa)
        --twoBitReference           2bit file of soft masked reference genome (default /mnt/expressions/martin/sequence_db/genomes/hg18.2bit)
        --tmp                       Path to tmp directory (default /tmp/)
    -o, --outfile                   Name of the resulting PDF file (default LibCompEst.pdf)
        --keep                      Keep R input file and log

"""

import sys
import os
import math
import random
import time
from library import entropy
from library import read_fastq
from subprocess import Popen, PIPE
from optparse import OptionParser

entropy_cutoff = 0.8
masked_cutoff = 0.5
merged_only = False

def frac_masked(seq):
  countLower = 0
  for x in seq:
    if x.islower(): countLower+=1
  if len(seq) > 0:
    return countLower/float(len(seq))
  else:
    return None

def get_masking_from_genome(idlist,genome2bit,binary,tmpfolder):
  global masked_cutoff
  tmpfile = tmpfolder+"/"+str(time.time()).replace(".","_")
  outfile = open(tmpfile+".ids",'w')
  for elem in idlist:
    outfile.write(elem.rstrip('+-')+"\n")
  outfile.close()
  #print "%s -seqList=%s.ids %s %s.fa"%(binary,tmpfile,genome2bit,tmpfile)
  print "Extracting sequences from the reference genome..."
  proc = Popen("%s -seqList=%s.ids %s %s.fa"%(binary,tmpfile,genome2bit,tmpfile),shell=True)
  proc.wait()
  os.remove(tmpfile+".ids")
  removelist = set()
  print "Evaluating repeat masking status..."
  for cid,cseq,cqual in read_fastq(tmpfile+".fa"):
    if frac_masked(cseq) >= masked_cutoff:
      #print frac_masked(cseq),cseq
      if cid+"+" in idlist: removelist.add(cid+"+")
      elif cid+"-" in idlist: removelist.add(cid+"-")
      else: print cid+" not found"
  os.remove(tmpfile+".fa")
  return removelist


def read_bowtie(filename):
  global merged_only
  infile = open(filename)
  idx = 0
  for line in infile:
    idx += 1
    fields = line.split()
    try:
      if (len(fields) > 6) and (fields[6]=='0'):
        if (not merged_only) or (merged_only and fields[0].startswith("M_")):
          yield fields[2],fields[3],fields[1],fields[4],line
        else:
          yield idx,None,None,None,None
      else:
        yield idx,None,None,None,None
    except (IOError,ValueError,IndexError):
      print "Error parsing BOWTIE output file in line:",line.strip()
  infile.close()

def read_eland(filename):
  infile = open(filename)
  idx = 0
  for line in infile:
    idx += 1
    try:
      fields = line.split()
      if len(fields) > 3:
        if fields[2] in ["U0","U1","U2"]: #OLD ELAND FILE UNIQUELY PLACED
          if fields[8] == "F": fields[8] = "+"
          else: fields[8] = "-"
          yield fields[6].split(".")[0],fields[7],fields[8],fields[1],line
        elif fields[2] in ["QC","NM"]:
          pass
        elif fields[2] in ["R0","R1","R2"]:
          yield idx,None,None,None,None
        else: #NEW ELAND FILE
          counts = fields[2].split(":")
          if (len(counts) == 3) and ( (counts[0] == "1") or
                                     ((counts[0] == "0") and (counts[1] == "1")) or 
                                     ((counts[0] == "0") and (counts[1] == "0") and (counts[2] == "1"))): ## UNIQUELY PLACED
            help_clean = (fields[3].split(",")[0]).strip("ACGT0123456789").split(":")
            strand = "-"
            if help_clean[1].endswith('F'): strand = "+"
            yield help_clean[0].split(".")[0],help_clean[1][:-1],strand,fields[1],line
          else:
            yield idx,None,None,None,None
    except (IOError,ValueError,IndexError):
      print "Error parsing ELAND output file in line:",line.strip()
  raise StopIteration
  infile.close()

def read_sam(filename):
  global merged_only
  strand_flag = 16
  first_read_flag = 64
  fid = 0
  fflag = 1
  fchr = 2
  fstart = 3
  fmqual = 4
  fcigar = 5
  finsert_size = 8
  fseq = 9
  fqual = 10
  infile = open(filename)
  idx = 0
  #count_lines = 0
  for line in infile:
    #count_lines+=1
    fields = line.split()
    try:
      insert_size = int(fields[finsert_size])
      if insert_size < 0: continue
      if (fields[fchr] == "*" or fields[fcigar] == "*"): continue
      if (len(fields) > fqual) and (fields[fmqual]!="0"):
        idx += 1
        if int(fields[fmqual]) >= 30:
          if (not merged_only) or (merged_only and fields[fid].startswith("M_")):
            strand = "+"
            if int(fields[fflag]) & strand_flag == strand_flag: strand = "-"
            if insert_size > 0:
              yield fields[fchr],fields[fstart],strand,fields[fseq]+"N"*max(0,insert_size-len(fields[fseq])),line
            else:
              yield fields[fchr],fields[fstart],strand,fields[fseq],line
          else:
            yield idx,None,None,None,None
        else:
          yield idx,None,None,None,None
    except (IOError,ValueError,IndexError):
      print "Error parsing SAM output file in line:",line.strip()
  #print idx,count_lines
  infile.close()


parser = OptionParser()
parser.add_option("--type", dest="filetype", help="Input file type (eland,bowtie,sam)",choices=["eland","bowtie","sam"],default="bowtie")
parser.add_option("-e", "--entropy", dest="entropy", help="Entropy cutoff to be used for sequence filtering (default = %.2f)"%entropy_cutoff,default=entropy_cutoff,type="float")
parser.add_option("-l", "--no_length", dest="nolength", help="Do not consider the length of the alignment as additional information",default=False,action="store_true")
parser.add_option("-f", "--min_length", dest="min_length", help="Do not consider reads shorter than N nt",type="int")
parser.add_option("-t", "--no_mt", dest="no_mthits", help="Do not consider hits to MT",default=False,action="store_true")
parser.add_option("-m", "--merged_only", dest="merged_only", help="Consider only merged reads from sam/bowtie files",default=False,action="store_true")
parser.add_option("-c", "--coverage_filter", dest="coverage", help="Do not consider very frequent coordinates in estimate (in percentiles: 0-100,default 100=OFF)",default=100,type="int")
parser.add_option("--outfile_coverage_filter",dest="cov_outfile", help="Print complexity filtered sequences to file (default OFF)")
parser.add_option("-r", "--repeat_mask", dest="repeat_mask", help="Do not consider reads of which at least FRACTION overlaps a repeat masked region (does not work in combination with no_length!)",type="float")
parser.add_option("--twoBitToFa", dest="twoBitToFa", help="Path to twoBitToFa binary (default /mnt/solexa/bin/blat64bit/twoBitToFa)",default="/mnt/solexa/bin/blat64bit/twoBitToFa")
parser.add_option("--twoBitReference", dest="twoBitReference", help="2bit file of soft masked reference genome (default /mnt/expressions/martin/sequence_db/genomes/hg18.2bit)",default="/mnt/expressions/martin/sequence_db/genomes/hg18.2bit")
parser.add_option("--tmp", dest="tmp", help="Path to tmp directory (default /tmp/)",default="/tmp/")
parser.add_option("--keep", dest="keep", help="Keep R input file and log",default=False,action="store_true")
parser.add_option("-o", "--outfile", dest="outfile", help="Name of the resulting PDF file (default LibCompEst.pdf)",default="LibCompEst.pdf")
(options, args) = parser.parse_args()

file_reader = None
if options.filetype == "bowtie": file_reader = read_bowtie
elif options.filetype == "sam": file_reader = read_sam
elif options.filetype == "eland": file_reader = read_eland
else:
  print "Unimplemented type for filter file."
  sys.exit()

entropy_cutoff = options.entropy
if (options.coverage < 0) or (options.coverage > 99): options.coverage=None
merged_only = options.merged_only

if (options.repeat_mask != None) and (options.repeat_mask <= 1.0) and (options.repeat_mask >= 0.0):
  masked_cutoff = options.repeat_mask
if (options.repeat_mask != None):
  if not os.path.exists(options.twoBitReference):
    print "Need valid reference genome file"
    sys.exit()
  if not os.path.exists(options.twoBitToFa):
    print "Need valid path to twoBitToFa binary"
    sys.exit()
  if not os.path.isdir(options.tmp):
    print "Need valid path to tmp directory"
    sys.exit()

count = 0
count_mt = 0
count_length = 0
count_reads = 0
count_complex = 0
all_coords = []
filenames = []
runs = []
coverage = {}
for elem in args:
  if os.path.exists(elem) or elem == "/dev/stdin":
      filenames.append(elem.split("/")[-1])
      for prun in elem.split("/"):
        if (len(prun) > 6) and (prun[6] == "_") and (("SOLEXA" in prun) or ("HWI" in prun) or ("BIOLAB" in prun) or ("SL-" in prun)):
          if prun not in runs: runs.append(prun)
      print "Reading",elem
      for chrom,start,strand,seq,line in file_reader(elem):
        count_reads += 1
        if seq != None:
          count += 1
          if entropy(seq) > entropy_cutoff:
            count_complex += 1
            cid = chrom+":"+start+strand
            if not options.nolength:
              cid = chrom+":"+start+"-"+str(int(start)+len(seq))+strand
            if "chrM" in cid: 
              count_mt+=1
              continue
            if ((options.min_length != None) and (len(seq) < options.min_length)):
              count_length+=1
              continue
            all_coords.append(cid)
            if cid in coverage:
              coverage[cid]+=1
            else:
              coverage[cid]=1

if len(all_coords) == 0:
  print "No input data. Stopping..."
  sys.exit()

count_repeat = 0
if (options.repeat_mask != None):
  set_coords = set(all_coords)
  print "Checking repeat masking status for %d coordinates (may take some time)..."%len(set_coords)
  removelist = get_masking_from_genome(set_coords,options.twoBitReference,options.twoBitToFa,options.tmp)
  del set_coords
  print "Removing %d coordinates based on repeat masking status of the genome"%len(removelist)
  nall_coords = []
  ncoverage = {}
  for elem in all_coords:
    if elem not in removelist:
      nall_coords.append(elem)
      if elem in ncoverage:
        ncoverage[elem]+=1
      else:
        ncoverage[elem]=1
    else:
      count_repeat+=1
  all_coords = nall_coords
  coverage = ncoverage

print "Shuffling data..."
random.shuffle(all_coords)
stops=[]
for elem in [1,2,3,4,5,6,7,8,9]:
  stops.append(elem*len(all_coords)//10)
stops.append(None)

sampled = []
unique = []
outtable = "\n\nSampled     Unique\n"
for stop in stops:
  sampled.append(len(all_coords[:stop]))
  unique.append(len(set(all_coords[:stop])))
  outtable+="%10d %10d\n"%(sampled[-1],unique[-1])

R_string = "library(gplots)\npdf(\"%s\",width=8,height=12)\n"%options.outfile
if len(runs) > 0:
  text = "Summary report for library complexity estimate\n==============================================\n\nRun:\n  "+"\n  ".join(runs)+"\n\nInput files:\n  "+"\n  ".join(filenames)+"\n\n"
else:
  text = "Summary report for library complexity estimate\n==============================================\n\nInput files:\n  "+"\n  ".join(filenames)+"\n\n"
text += "Number of mapped sequences: %d\n"%(count_reads)
text += "Number of uniquely placed sequences: %d (%0.2f%%)\n"%(count,count/float(count_reads)*100)
text += "Number passed sequence complexity filter: %d (%0.2f%%)\n"%(count_complex,count_complex/float(count_reads)*100)
if count_mt > 0: text += "Number of filtered MT sequences: %d (%0.2f%%)\n"%(count_mt,(count_complex-count_mt)/float(count_reads)*100)
if count_length > 0: text += "Number of length filtered sequences (< %d nt): %d (%0.2f%%)\n"%(options.min_length,count_length,(count_complex-count_mt-count_length)/float(count_reads)*100)
if count_repeat > 0: text += "Number of repeat filtered sequences (%.2f%%): %d (%0.2f%%)\n"%(masked_cutoff*100.0,count_repeat,(count_complex-count_mt-count_length-count_repeat)/float(count_reads)*100)
if sampled[-1] != unique[-1]:
  if sampled[-1] < unique[-1]*1.2:
    R_string+= "textplot(%r,halign='left',valign='top',fixed.width=T,tab.width=8,mar=c(2,2,2,2))\n"%(text+outtable+"\n\nWARNING: Library probably not sequenced deep enough!\n")
  else:
    R_string+= "textplot(%r,halign='left',valign='top',fixed.width=T,tab.width=8,mar=c(2,2,2,2))\n"%(text+outtable)
    R_string+= "M <- NA\n"
  R_string+= "sampled = c("+", ".join(map(str,sampled))+")\n"
  R_string+= "sunique = c("+", ".join(map(str,unique))+")\n"
  R_string+= "tsampled = %d\n"%sampled[-1]
  R_string+= "tsunique = %d\n"%unique[-1]
  R_string+= "cfit<-nls(sunique ~ 0+ M*(1-exp(-sampled/M)),start = list(M=%d),control=list(warnOnly=T))\n"%unique[-1]
  R_string+= "summary(cfit)\n"
  R_string+= "print(paste('Sampled',as.character(sunique[length(sunique)]),sep=' '))\n"
  R_string+= "tdata <- c(unlist(sampled),%d*2,%d*3,%d*4,%d*5,%d*6,%d*7,%d*8,%d*9,%d*10)\n"%tuple([sampled[-1]]*9)
  R_string+= "par(mfrow=c(2,1))\n"
  R_string+= "textplot(paste('Output of fitting:\\n==================\\n\\n',gsub('\\U2019','',gsub('\\U2018','',paste(capture.output(summary(cfit)),sep='',collapse='\\n'))),sep=''),halign='left',valign='center',fixed.width=T,tab.width=8,mar=c(0,0,0,0))\n"
  R_string+= "plot(tdata,coef(cfit)*(1-exp(-tdata/coef(cfit))),type=\"l\",col=2,xlab=\"Sampled\",ylab=\"Uniquely sampled\",main=paste(\"Estimated number of unique molecules in library M = \",as.character(round(coef(cfit))),sep=\"\"),mar=c(2,2,2,2))\n"
  R_string+= "points(sampled,sunique)\n"
  R_string+= "par(mfrow=c(1,1))\n"
else:
  R_string+= "textplot(%r,halign='left',valign='top',fixed.width=T,tab.width=8,mar=c(2,2,2,2))\n"%(text+outtable+"\n\nERROR: Library not sequenced deep enough for estimate!")

if options.coverage != None:
  old_max_unique = unique[-1]
  print "Applying coverage filter to data..."
  cov_cutoff = None
  values = map(lambda (x,y):(y,x),coverage.iteritems())
  values.sort()
  pmin = values[0][0]
  p5 = values[int(len(values)*(0.05))][0]
  p25 = values[int(len(values)*(0.25))][0]
  p50 = values[int(len(values)*(0.50))][0]
  p75 = values[int(len(values)*(0.75))][0]
  p95 = values[int(len(values)*(0.95))][0]
  pmax = values[-1][0]
  mean = sum(map(lambda (x,y):x,values))/float(len(values))
  cov_hist = []
  count_hist = []
  for ccov,ccoord in values:
    if len(cov_hist) == 0:
      cov_hist.append(ccov)
      count_hist.append(1)
    else:
      if ccov == cov_hist[-1]:
        count_hist[-1]+=1
      else:
        cov_hist.append(ccov)
        count_hist.append(1)
  cov_stat = "\n\nCoverage:\n- Min, Mean, Max:\n   %6d %6.2f %6d\n- 5th, 25th, 50th, 75th, 95th percentile:\n   %6d %6d %6d %6d %6d"%(pmin,mean,pmax,p5,p25,p50,p75,p95)
  R_string+= "cov_hist = c("+", ".join(map(str,cov_hist))+")\n"
  R_string+= "count_hist = c("+", ".join(map(str,count_hist))+")\n"
  cov_cutoff = values[int(len(values)*(options.coverage/100.0))][0]
  rem_coords = set()
  for ccov,ccoord in values[int(len(values)*(options.coverage/100.0)):]:
    if ccov > cov_cutoff:
      rem_coords.add(ccoord)
  nall_coords=[]
  if options.cov_outfile != None:
    print "Extracting coverage filtered sequences..."
    outfile = open(options.cov_outfile,'w')
    for elem in args:
      if os.path.exists(elem):
        for chrom,start,strand,seq,line in file_reader(elem):
          if seq != None:
            cid = chrom+":"+start+strand
            if not options.nolength:
              cid = chrom+":"+start+"-"+str(int(start)+len(seq))+strand
            if cid in rem_coords:
              outfile.write(line)
    outfile.close()
  for elem in all_coords:
    if elem not in rem_coords:
      nall_coords.append(elem)
  all_coords=nall_coords
  text += "Number passed coverage filter [coverage <= %d]: %d (%0.2f%%)"%(cov_cutoff,len(all_coords),len(all_coords)/float(count_reads)*100)

  stops=[]
  for elem in [1,2,3,4,5,6,7,8,9]:
    stops.append(elem*len(all_coords)//10)
  stops.append(None)

  sampled = []
  unique = []
  outtable = "\n\nSampled     Unique\n"
  for stop in stops:
    sampled.append(len(all_coords[:stop]))
    unique.append(len(set(all_coords[:stop])))
    outtable+="%10d %10d\n"%(sampled[-1],unique[-1])

  R_string+= "sampled = c("+", ".join(map(str,sampled))+")\n"
  R_string+= "sunique = c("+", ".join(map(str,unique))+")\n"
  R_string+= "tsampled = %d\n"%sampled[-1]
  R_string+= "tsunique = %d\n"%unique[-1]

  R_string+= "layout(matrix(1:3,3,1),heights=c(30,30,30))\n"
  R_string+= "textplot(%r,halign='left',valign='top',fixed.width=T,tab.width=8,mar=c(0,2,2,2))\n"%(text+cov_stat)
  R_string+= "plot(count_hist~cov_hist,xlab=\"Coverage\",ylab=\"Number of molecules\",type='h',lwd=3,mar=c(0,10,0,10))\n"
  if sampled[-1] == unique[-1]:
    R_string+= "textplot(%r,halign='left',valign='top',fixed.width=T,tab.width=8,mar=c(2,2,0,2))\n"%(outtable+"\n\nERROR: Library not sequenced deep enough for estimate!")
    R_string+= "M <- NA\n"
  elif sampled[-1] >= unique[-1]*1.2:
      R_string+= "textplot(%r,halign='left',valign='top',fixed.width=T,tab.width=8,mar=c(2,2,0,2))\n"%(outtable)
  else:
      R_string+= "textplot(%r,halign='left',valign='top',fixed.width=T,tab.width=8,mar=c(2,2,0,2))\n"%(outtable+"\n\nWARNING: Library probably not sequenced deep enough!\n")
  if sampled[-1] != unique[-1]:
    R_string+= "cfit<-nls(sunique ~ 0+ M*(1-exp(-sampled/M)),start = list(M=%d),control=list(warnOnly=T))\n"%unique[-1]
    R_string+= "summary(cfit)\n"
    R_string+= "print(paste('Sampled',as.character(sunique[length(sunique)]),sep=' '))\n"
    R_string+= "tdata <- c(unlist(sampled),%d*2,%d*3,%d*4,%d*5,%d*6,%d*7,%d*8,%d*9,%d*10)\n"%tuple([sampled[-1]]*9)
    R_string+= "par(mfrow=c(2,1))\n"
    R_string+= "textplot(paste('Output of fitting:\\n==================\\n\\n',gsub('\\U2019','',gsub('\\U2018','',paste(capture.output(summary(cfit)),sep='',collapse='\\n'))),'\\nCorrected estimate for M (considering coverage filtering):\\n',as.character(round(coef(cfit)*%f)),sep=''),halign='left',valign='center',fixed.width=T,tab.width=8,mar=c(0,2,2,2))\n"%(old_max_unique/float(unique[-1]))
    R_string+= "plot(tdata,coef(cfit)*(1-exp(-tdata/coef(cfit))),type=\"l\",col=2,xlab=\"Sampled\",ylab=\"Uniquely sampled\",main=paste(\"Estimated number of unique molecules in library M = \",as.character(round(coef(cfit))),sep=\"\"),mar=c(2,2,0,2))\n"
    R_string+= "points(sampled,sunique)\n"

R_string+= "dev.off()\n"

print "Handing over to R..."
proc=Popen("/usr/bin/env R --vanilla --quiet",stdin=PIPE,stderr=PIPE,stdout=PIPE,shell=True)
Rlog = proc.communicate(input=R_string)[0]
proc.wait()
if options.keep:
  outfile =  open(".".join(options.outfile.split(".")[:-1])+".R",'w')
  outfile.write(R_string)
  outfile.close()

  outfile =  open(".".join(options.outfile.split(".")[:-1])+".log",'w')
  outfile.write(Rlog)
  outfile.close()

#print R_string
