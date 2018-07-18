#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Align reads stored as BAM using BWA

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *10.02.2012
:Type: tool
:Input: BAM
:Output: BAM

"""

import sys,os
import pysam
import subprocess
import random
import traceback

import string
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv') # COMPLEMENT DNA

def join_bam_entry(entry,additional):
  if additional.is_duplicate: entry.is_duplicate = True
  if additional.is_qcfail: entry.is_qcfail = True
  tags = entry.tags
  if tags == None: tags = []
  keys = set(map(lambda (x,y): x,tags))
  for key,value in additional.tags:
    if (key not in keys) and (key in ['ZQ','RG','YI','YJ','XI','XJ']): tags.append((key,value))
  if len(tags) > 0:
    entry.tags = tags
  else:
    entry.tags = None
  return entry

def write_fastq(out,ostream,read):
  global table
  if out == 0: ostream[-2] += 1
  else: ostream[-1] += 1
  if read.is_reverse:
    ostream[out].write("@%s\n%s\n+\n%s\n"%(read.qname,(read.seq[::-1]).translate(table),read.qual[::-1]))
  else:
    ostream[out].write("@%s\n%s\n+\n%s\n"%(read.qname,read.seq,read.qual))

defaultIndexName = "bwa-0.4.9"
from optparse import OptionParser,OptionGroup
parser = OptionParser("%prog [options] BAMfile")
parser.add_option("--mock", dest="mock", help="Don't do anything, do a mock run",action="store_true",default=False)
parser.add_option("-p","--PIPE",dest="pipe",help="Write BAM to PIPE",default=False,action="store_true")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
parser.add_option("-c", "--cores", dest="cores", help="Number of cores to be used in the alignment step (default 1)",default=1,type="int")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output file in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output file (default 'Align').",default="Align")
parser.add_option("--bwa", dest="bwa_bin", help="Path to bwa binary (default /net/shendure/vol1/home/mkircher/bin/bwa)",default="/net/shendure/vol1/home/mkircher/bin/bwa")
parser.add_option("--bwa_genome", dest="bwa_genome", help="Path to genome (folder containing %s* and whole_genome.fa/.fa.fai) (default '')"%(defaultIndexName),default=defaultIndexName)
parser.add_option("--bwa_params", dest="bwa_params", help="Parameters passed to bwa (default '')",default="")
parser.add_option("--indexName", dest="indexName", help="Prefix of genome index (default '%s')"%(defaultIndexName),default=defaultIndexName)
parser.add_option("--only_aligned", dest="only_aligned", help="Skip not aligned reads in BAM conversion",default=False,action="store_true")
parser.add_option("-v","--verbose",dest="verbose",help="Print progress messages",default=False,action='store_true')
(options, args) = parser.parse_args()

if options.cores < 1: options.cores = 1

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

options.bwa_params = options.bwa_params.replace('_',' ')

if not os.path.isfile(options.bwa_genome+"/%s.bwt"%options.indexName) or not os.path.isfile(options.bwa_genome+"/whole_genome.fa.fai"):
  sys.stderr.write("Check genome index and fasta index.\n")
  sys.exit()

filename = args[0]
if len(args) == 0 or not os.path.exists(args[0]):
  sys.stderr.write("Input file is not available: %s\n"%(filename))
  sys.exit()

infile = pysam.Samfile( filename, 'rb' )
cheader = infile.header
outfilename = "-"
fileflags = 'wb'
if options.SAM: fileflags = 'w'

try:
  SRproc_,SRproc,PEproc,forward,reverse = None,None,None,None,None

  #FORWARD READ ALIGNMENT
  cjob = options.bwa_bin+" aln -b1 %s -t %d %s/%s %s"%(options.bwa_params,options.cores,options.bwa_genome,options.indexName,filename)
  if options.mock:
    sys.stderr.write(cjob+"\n")
  else:
    forward = subprocess.Popen(cjob.split(),bufsize=0,stdout=subprocess.PIPE)

  #REVERSE READ ALIGNMENT
  cjob = options.bwa_bin+" aln -b2 %s -t %d %s/%s %s"%(options.bwa_params,options.cores,options.bwa_genome,options.indexName,filename)
  if options.mock:
    sys.stderr.write(cjob+"\n")
  else:
    reverse = subprocess.Popen(cjob.split(),bufsize=0,stdout=subprocess.PIPE)

  #RUN PE ALIGNMENT STEP
  cjob = options.bwa_bin+" sampe -a 800 %s/%s /dev/fd/%d /dev/fd/%d %s %s"%(options.bwa_genome,options.indexName,forward.stdout.fileno(),reverse.stdout.fileno(),filename,filename)
  if options.mock:
    sys.stderr.write(cjob+"\n")
  else:
    PEproc = subprocess.Popen(cjob.split(),bufsize=0,stdout=subprocess.PIPE)
    forward.stdout.close()
    reverse.stdout.close()


  #RUN SR ALIGNMENT STEP
  cjob1 = options.bwa_bin+" aln -b0 %s -t %d %s/%s %s"%(options.bwa_params,options.cores,options.bwa_genome,options.indexName,filename)
  cjob2 = options.bwa_bin+" samse %s/%s - %s"%(options.bwa_genome,options.indexName,filename)
  if options.mock:
    sys.stderr.write(cjob1+" | "+cjob2+"\n")
  else:
    SRproc_ = subprocess.Popen(cjob1.split(),bufsize=0,stdout=subprocess.PIPE)
    SRproc = subprocess.Popen(cjob2.split(),bufsize=0,stdin=SRproc_.stdout,stdout=subprocess.PIPE)
    SRproc_.stdout.close()

  if PEproc != None: PEalign = pysam.Samfile( "/dev/fd/%d"%(PEproc.stdout.fileno()), 'r')
  else: PEalign = None
  if SRproc != None: SRalign = pysam.Samfile( "/dev/fd/%d"%(SRproc.stdout.fileno()) , 'r')
  else: SRalign = None
  SQs = []
  SRread,PEread = None,None
  if PEalign != None:
    cheader['SQ'] = PEalign.header['SQ']
    try: PEread = PEalign.next()
    except: PEread = None
  elif SRalign != None:
    cheader['SQ'] = SRalign.header['SQ']

  if SRalign != None:
    try: SRread = SRalign.next()
    except: SRread = None

  # BAM OUTPUT STREAM
  if options.pipe:
    bam_outfile = pysam.Samfile( outfilename, fileflags, header = cheader)
    if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
  else:
    outfilename = options.outdir+options.outprefix+("_"+options.bwa_params.replace(" ","_")+"_"+options.bwa_genome.strip('/').split('/')[-1]).replace("__","_")+".bam"
    if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
    bam_outfile = pysam.Samfile( outfilename, fileflags, header = cheader)

  readbuffer = {}
  for read in infile:
    if read.is_paired and (PEread != None) and (PEread.qname == read.qname) and (PEread.is_read1 == read.is_read1):
      if (not options.only_aligned) or (not PEread.is_unmapped):
        bam_outfile.write(join_bam_entry(PEread,read))
      PEread = None
    elif not read.is_paired and (SRread != None) and (SRread.qname == read.qname):
      if (not options.only_aligned) or (not SRread.is_unmapped):
        bam_outfile.write(join_bam_entry(SRread,read))
      SRread = None
    elif read.is_paired and (PEread != None):
      readbuffer[read.qname,read.is_read1] = read
      if (PEread.qname,PEread.is_read1) in readbuffer:
        read = readbuffer.pop((PEread.qname,PEread.is_read1))
        if (not options.only_aligned) or (not PEread.is_unmapped):
          bam_outfile.write(join_bam_entry(PEread,read))
        PEread = None
    else:
      sys.stderr.write("Error retrieving original BAM entry! Should not happen.\n")
      raise RuntimeError
    if SRalign != None and SRread == None:
      try: SRread = SRalign.next()
      except: SRread = None
    if PEalign != None and PEread == None:
      try: PEread = PEalign.next()
      except: PEread = None
  infile.close()

  # CLEAN UP PIPES AND FILE STREAMS
  bam_outfile.close()

  if SRproc_ != None: SRproc_.wait()
  if SRproc != None: SRproc.wait()
  if forward != None: forward.wait()
  if reverse != None: reverse.wait()
  if PEproc != None: PEproc.wait()

  if PEalign != None: PEalign.close()
  if SRalign != None: SRalign.close()
except:
  sys.stderr.write("Unexpected termination of BWA2BAM. Cleaning up files...\n")
  exc_type, exc_value, exc_traceback = sys.exc_info()
  sys.stderr.write("%s\n"%str(exc_value))
  traceback.print_tb(exc_traceback)
  if outfilename != "-":
    try:
      os.remove(outfilename)
    except:
      sys.stderr.write("Error removing output file%s\n"%(outfilename))
