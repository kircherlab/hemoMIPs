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
import string
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv') # COMPLEMENT DNA
import gzip

def write_fastq(out,ostream,read):
  if (out < 3):
    if out == 0: ostream[-2] += 1
    else: ostream[-1] += 1
    if read.is_reverse:
      ostream[out].write("@%s\n%s\n+\n%s\n"%(read.qname,(read.seq[::-1]).translate(table),read.qual[::-1]))
    else:
      ostream[out].write("@%s\n%s\n+\n%s\n"%(read.qname,read.seq,read.qual))
  elif out == 3:
    seq = ""
    qual = ""
    for key,value in read.tags:
      if key == "XI": seq = value
      elif key == "YI": qual = value
    ostream[out].write("@%s\n%s\n+\n%s\n"%(read.qname,seq,qual))
  elif out == 4:
    seq = ""
    qual = ""
    for key,value in read.tags:
      if key == "XJ": seq = value
      elif key == "YJ": qual = value
    ostream[out].write("@%s\n%s\n+\n%s\n"%(read.qname,seq,qual))
    
    
from optparse import OptionParser,OptionGroup
parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from PIPE",default=False,action="store_true")
parser.add_option("--PIPE2",dest="pipe2",help="Read BAM from PIPE and write FastQ to PIPE",default=False,action="store_true")
parser.add_option("--RG",dest="RG",help="Generate one output file per read group",default=False,action="store_true")
parser.add_option("--RG2",dest="RG2",help="Generate one output file per read group and exclude unknown, wrong, conflict",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'FastQ').",default="FastQ")
parser.add_option("", "--SR", dest="SingleReads", help="Do not output PE data.",default=False,action="store_true")
parser.add_option("-a","--all",dest="all",help="Report all reads (i.e. not only QC-Passed reads)",default=False,action="store_true")
parser.add_option("-e","--include_incomplete_pairs",dest="incl_inc_pairs",help="Include incomplete pairs",default=False,action="store_true")
parser.add_option("", "--no_clean", dest="no_clean", help="Do not remove empty output files",default=False,action="store_true")
parser.add_option("", "--one_length", dest="one_length", help="Do only keep reads of the same length as the first entry",default=False,action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
parser.add_option("-f","--filter",dest="filterlst",help="Filter list of read ids",default=None)
parser.add_option("-I", dest="keepI", help="Write contents of BAM field XI/YI to separate I1 file",default=False,action="store_true")
parser.add_option("-J", dest="keepJ",help="Write contents of BAM field XI/YI to separate I2 file",default=False,action="store_true")
parser.add_option("--gzip", dest="gzip",help="Use GZIP compression for output files",default=False,action="store_true")
(options, args) = parser.parse_args()

filterset = None
if options.filterlst != None and os.path.exists(options.filterlst):
  filterset = set()
  infile = open(options.filterlst)
  for line in infile:
    filterset.add(line.strip())
  infile.close()
  sys.stderr.write('Loaded %d read IDs for filtering...\n'%(len(filterset)))

if options.RG2: options.RG = True

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

if options.pipe2:
  if options.incl_inc_pairs:
    sys.stderr.write("Filtering of incomplete pairs not available with --PIPE2.\n")
  if options.RG:
    sys.stderr.write("Output files per RG not available with --PIPE2.\n")
  options.pipe = True
  options.incl_inc_pairs = True
  options.RG = False

files = args
if options.pipe: files = [None]
outfiles = {}

readl1,readl2 = None,None

for filename in files:
  if filename == None:
    if options.verbose: sys.stderr.write("Reading binary BAM from STDIN...\n")
    infile = pysam.Samfile( "-", 'rb' )
  elif os.path.exists(filename):
    infile = pysam.Samfile( filename, 'rb' )
  else: continue

  if options.pipe2:
    if options.verbose: sys.stderr.write("Writing to STDOUT...\n")
    outfiles[None] = [sys.stdout,sys.stdout,sys.stdout,None,None,None,None,None,None,None,0,0]
  elif options.RG:
    if options.verbose: sys.stderr.write("Creating output files per RG...\n")
  else:
    outfiles[None] = [None,None,None,None,None,None,None,None,None,None,0,0]
    if options.verbose: sys.stderr.write("Creating output files...\n")
    if options.gzip: 
      outfilename = options.outdir+options.outprefix+"_SR.fq.gz"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfiles[None][0] = gzip.open(outfilename,'w')
      outfiles[None][5] = outfilename
    else: 
      outfilename = options.outdir+options.outprefix+"_SR.fq"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfiles[None][0] = open(outfilename,'w')
      outfiles[None][5] = outfilename
    
    if not options.SingleReads:
      if options.gzip:
        outfilename = options.outdir+options.outprefix+"_r1.fq.gz"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][1] = gzip.open(outfilename,'w')
        outfiles[None][6] = outfilename
        outfilename = options.outdir+options.outprefix+"_r2.fq.gz"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][2] = gzip.open(outfilename,'w')
        outfiles[None][7] = outfilename
      else:
        outfilename = options.outdir+options.outprefix+"_r1.fq"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][1] = open(outfilename,'w')
        outfiles[None][6] = outfilename
        outfilename = options.outdir+options.outprefix+"_r2.fq"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][2] = open(outfilename,'w')
        outfiles[None][7] = outfilename
    if options.keepI:
      if options.gzip:
        outfilename = options.outdir+options.outprefix+"_I1.fq.gz"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][3] = gzip.open(outfilename,'w')
        outfiles[None][8] = outfilename
      else:
        outfilename = options.outdir+options.outprefix+"_I1.fq"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][3] = open(outfilename,'w')
        outfiles[None][8] = outfilename
    if options.keepJ:
      if options.gzip:
        outfilename = options.outdir+options.outprefix+"_I2.fq.gz"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][4] = gzip.open(outfilename,'w')
        outfiles[None][9] = outfilename        
      else:
        outfilename = options.outdir+options.outprefix+"_I2.fq"
        if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
        outfiles[None][4] = open(outfilename,'w')
        outfiles[None][9] = outfilename

  if not options.incl_inc_pairs and not options.SingleReads and not options.all: # SEPERATE EXCL INCL PAIRS CASE FOR SPEED REASONS
    incomplete_pairs = {}
    for read in infile:
      RG = None
      if options.RG:
        RG = "unknown"
        for (key,value) in read.tags:
          if key == "RG":
            RG = value
            if RG not in outfiles:
              outfiles[RG] = [None,None,None,None,None,None,None,None,None,None,0,0]
              if options.RG2 and RG in ["unknown","wrong","conflict"]: break
              if options.gzip:
                outfilename = options.outdir+options.outprefix+"_%s_SR.fq.gz"%(RG)
                outfiles[RG][0] = gzip.open(outfilename,'w')
                outfiles[RG][5] = outfilename
              else:
                outfilename = options.outdir+options.outprefix+"_%s_SR.fq"%(RG)
                outfiles[RG][0] = open(outfilename,'w')
                outfiles[RG][5] = outfilename
              if not options.SingleReads:
                if options.gzip:
                  outfilename = options.outdir+options.outprefix+"_%s_r1.fq.gz"%(RG)
                  outfiles[RG][1] = gzip.open(outfilename,'w')
                  outfiles[RG][6] = outfilename
                  outfilename = options.outdir+options.outprefix+"_%s_r2.fq.gz"%(RG)
                  outfiles[RG][2] = gzip.open(outfilename,'w')
                  outfiles[RG][7] = outfilename
                else:
                  outfilename = options.outdir+options.outprefix+"_%s_r1.fq"%(RG)
                  outfiles[RG][1] = open(outfilename,'w')
                  outfiles[RG][6] = outfilename
                  outfilename = options.outdir+options.outprefix+"_%s_r2.fq"%(RG)
                  outfiles[RG][2] = open(outfilename,'w')
                  outfiles[RG][7] = outfilename
              if options.keepI:
                if options.gzip:
                  outfilename = options.outdir+options.outprefix+"_%s_I1.fq.gz"%(RG)
                  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
                  outfiles[RG][3] = gzip.open(outfilename,'w')
                  outfiles[RG][8] = outfilename                  
                else:
                  outfilename = options.outdir+options.outprefix+"_%s_I1.fq"%(RG)
                  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
                  outfiles[RG][3] = open(outfilename,'w')
                  outfiles[RG][8] = outfilename
              if options.keepJ:
                if options.gzip:
                  outfilename = options.outdir+options.outprefix+"_%s_I2.fq.gz"%(RG)
                  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
                  outfiles[RG][4] = gzip.open(outfilename,'w')
                  outfiles[RG][9] = outfilename                  
                else:
                  outfilename = options.outdir+options.outprefix+"_%s_I2.fq"%(RG)
                  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
                  outfiles[RG][4] = open(outfilename,'w')
                  outfiles[RG][9] = outfilename
            break
      if options.RG2 and RG in ["unknown","wrong","conflict"]: continue
      if read.is_paired:
        if read.qname in incomplete_pairs:
          oread = incomplete_pairs[read.qname]
          if not oread.is_qcfail and not read.is_qcfail:
            if filterset == None or (read.qname in filterset):
              if oread.is_read1:
                if readl1 == None and options.one_length:
                  readl1 = len(oread.seq)
                  readl2 = len(read.seq)
                if (readl1 == None) or ((readl1 == len(oread.seq)) and (readl2 == len(read.seq))):
                  write_fastq(1,outfiles[RG],oread)
                  write_fastq(2,outfiles[RG],read)
                  if options.keepI:
                    write_fastq(3,outfiles[RG],oread)
                  if options.keepJ:
                    write_fastq(4,outfiles[RG],oread)
              else:
                if readl1 == None and options.one_length:
                  readl1 = len(read.seq)
                  readl2 = len(oread.seq)
                if (readl1 == None) or ((readl1 == len(read.seq)) and (readl2 == len(oread.seq))):
                  write_fastq(1,outfiles[RG],read)
                  write_fastq(2,outfiles[RG],oread)
                  if options.keepI:
                    write_fastq(3,outfiles[RG],read)
                  if options.keepJ:
                    write_fastq(4,outfiles[RG],read)
          del incomplete_pairs[read.qname]
        else:
          incomplete_pairs[read.qname] = read
      elif not read.is_qcfail:
        if filterset == None or (read.qname in filterset):
          if readl1 == None and options.one_length:
            readl1 = len(read.seq)
          elif (readl1 == None) or (readl1 == len(read.seq)):
            write_fastq(0,outfiles[RG],read)
            if options.keepI:
              write_fastq(3,outfiles[RG],read)
            if options.keepJ:
              write_fastq(4,outfiles[RG],read)
  else:
    for read in infile:
      if not read.is_qcfail or options.all:
        RG = None
        if options.RG:
          for (key,value) in read.tags:
            if key == "RG":
              RG = value
              if RG not in outfiles:
                outfiles[RG] = [None,None,None,None,None,None,None,None,None,None,0,0]
                if options.RG2 and RG in ["unknown","wrong","conflict"]: break
                if options.gzip:
                  outfilename = options.outdir+options.outprefix+"_%s_SR.fq.gz"%(RG)
                  outfiles[RG][0] = gzip.open(outfilename,'w')
                  outfiles[RG][5] = outfilename
                else:
                  outfilename = options.outdir+options.outprefix+"_%s_SR.fq"%(RG)
                  outfiles[RG][0] = open(outfilename,'w')
                  outfiles[RG][5] = outfilename
                if not options.SingleReads:
                  if options.gzip:
                    outfilename = options.outdir+options.outprefix+"_%s_r1.fq.gz"%(RG)
                    outfiles[RG][1] = gzip.open(outfilename,'w')
                    outfiles[RG][6] = outfilename
                    outfilename = options.outdir+options.outprefix+"_%s_r2.fq.gz"%(RG)
                    outfiles[RG][2] = gzip.open(outfilename,'w')
                    outfiles[RG][7] = outfilename
                  else:                     
                    outfilename = options.outdir+options.outprefix+"_%s_r1.fq"%(RG)
                    outfiles[RG][1] = open(outfilename,'w')
                    outfiles[RG][6] = outfilename
                    outfilename = options.outdir+options.outprefix+"_%s_r2.fq"%(RG)
                    outfiles[RG][2] = open(outfilename,'w')
                    outfiles[RG][7] = outfilename
                if options.keepI:
                  if options.gzip:
                    outfilename = options.outdir+options.outprefix+"_%s_I1.fq.gz"%(RG)
                    outfiles[RG][3] = gzip.open(outfilename,'w')
                    outfiles[RG][8] = outfilename                  
                  else:
                    outfilename = options.outdir+options.outprefix+"_%s_I1.fq"%(RG)
                    outfiles[RG][3] = open(outfilename,'w')
                    outfiles[RG][8] = outfilename
                if options.keepJ:
                  if options.gzip:
                    outfilename = options.outdir+options.outprefix+"_%s_I2.fq.gz"%(RG)
                    outfiles[RG][4] = gzip.open(outfilename,'w')
                    outfiles[RG][9] = outfilename                                      
                  else:
                    outfilename = options.outdir+options.outprefix+"_%s_I2.fq"%(RG)
                    outfiles[RG][4] = open(outfilename,'w')
                    outfiles[RG][9] = outfilename                  
              break
        if options.RG2 and RG in ["unknown","wrong","conflict"]: continue
        if options.pipe2: 
          if filterset == None or (read.qname in filterset):
            if readl1 == None and options.one_length:
              readl1 = len(read.seq)
            elif (readl1 == None) or (readl1 == len(read.seq)):
              write_fastq(0,outfiles[RG],read)
              if options.keepI:
                write_fastq(3,outfiles[RG],read)
              if options.keepJ:
                write_fastq(4,outfiles[RG],read)
        elif read.is_paired:
          if options.SingleReads: continue
          elif read.is_read1: 
            if filterset == None or (read.qname in filterset):
              if readl1 == None and options.one_length:
                readl1 = len(read.seq)
              if (readl1 == None) or (readl1 == len(read.seq)):
                write_fastq(1,outfiles[RG],read)
                if options.keepI:
                  write_fastq(3,outfiles[RG],read)
                if options.keepJ:
                  write_fastq(4,outfiles[RG],read)
          else: 
            if filterset == None or (read.qname in filterset):
              if readl2 == None and options.one_length:
                readl2 = len(read.seq)
              if (readl2 == None) or (readl2 == len(read.seq)):
                write_fastq(2,outfiles[RG],read)
        else:
          if filterset == None or (read.qname in filterset):
            if readl1 == None and options.one_length:
              readl1 = len(read.seq)
            if (readl1 == None) or (readl1 == len(read.seq)):
              write_fastq(0,outfiles[RG],read)
              if options.keepI:
                write_fastq(3,outfiles[RG],read)
              if options.keepJ:
                write_fastq(4,outfiles[RG],read)

if not options.pipe2:
  for cid,(ostreamSR,ostreamR1,ostreamR2,ostreamI1,ostreamI2,filenameSR,filenameR1,filenameR2,filenameI2,filenameI2,cSR,cPE) in outfiles.iteritems():
    if ostreamSR != None:
      ostreamSR.close()
      if cSR == 0 and not options.no_clean: os.remove(filenameSR)
    if ostreamR1 != None and ostreamR2 != None:
      ostreamR1.close()
      ostreamR2.close()
      if cPE == 0 and not options.no_clean:
        os.remove(filenameR1)
        os.remove(filenameR2)
    if ostreamI1 != None and ostreamI2 != None:
      ostreamI1.close()
      ostreamI2.close()
