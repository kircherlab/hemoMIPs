#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Extract reads to FastQ for BWA MEM alignment

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *14.05.2014
:Type: tool
:Input: BAM
:Output: FastQ

"""

import sys,os
import pysam
import string
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv') # COMPLEMENT DNA

#import re
#A = re.compile("^[!-~]$")                                                       # """Printable character"""
#i = re.compile("^[-+]?[0-9]+$")                                                 # """Singed 32-bit integer"""
#f = re.compile("^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$")                      # """Single-precision foating number"""
#H = re.compile("^[0-9A-F]+$")                                                   # """Byte array in the Hex format"""
#B = re.compile("^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$")         # """Integer or numeric array"""
#Z = re.compile("^[ !-~]+$")                                                     # """Printable string, including space"""
def format_tags(read):
  fields = []
  #global A,i,f,H,B,Z
  for (key,value) in read.tags:
    if type(value) == type(1):
      fields.append("%s:i:%d"%(key,value))
    elif type(value) == type(1.1):
      fields.append("%s:f:%.6f"%(key,value))
    elif type(value) == type("") and len(value) == 1:
      fields.append("%s:A:%s"%(key,value))
    else: # type(value) == type(""):
      fields.append("%s:Z:%s"%(key,value))
  return "\t".join(fields)

def write_fastq(read):
  if read.is_reverse:
    sys.stdout.write("@%s %s\n%s\n+\n%s\n"%(read.qname,format_tags(read),(read.seq[::-1]).translate(table),read.qual[::-1]))
  else:
    sys.stdout.write("@%s %s\n%s\n+\n%s\n"%(read.qname,format_tags(read),read.seq,read.qual))

from optparse import OptionParser,OptionGroup
parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from PIPE",default=False,action="store_true")
parser.add_option("", "--SR", dest="SingleReads", help="Only output Single Reads",default=False,action="store_true")
parser.add_option("", "--R1", dest="ForwardRead", help="Only output Forward Reads",default=False,action="store_true")
parser.add_option("", "--R2", dest="ReverseRead", help="Only output Reverse Reads",default=False,action="store_true")
parser.add_option("-a","--all",dest="all",help="Report all reads (i.e. not only QC-Passed reads, default On)",default=True,action="store_false")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
(options, args) = parser.parse_args()

files = args
if options.pipe: files = [None]
outfiles = {}

if options.SingleReads and options.ForwardRead:
  options.ForwardRead = False
  sys.stderr.write("Only reporting Single Reads.\n")
if (options.SingleReads and options.ReverseRead) or (options.ReverseRead and options.ForwardRead):
  options.ReverseRead = False
  sys.stderr.write("Only reporting Forward Reads.\n")
if not options.SingleReads and not options.ForwardRead and not options.ReverseRead:
  options.SingleReads = True
  sys.stderr.write("Only reporting Single Reads.\n")

for filename in files:
  if filename == None:
    if options.verbose: sys.stderr.write("Reading binary BAM from STDIN...\n")
    infile = pysam.Samfile( "-", 'rb' )
  elif os.path.exists(filename):
    infile = pysam.Samfile( filename, 'rb' )
  else: continue

  incomplete_pairs = {}
  for read in infile:
    if read.is_paired and (options.ForwardRead or options.ReverseRead):
      if read.qname in incomplete_pairs:
        oread = incomplete_pairs[read.qname]
        if options.all or (not oread.is_qcfail and not read.is_qcfail):
          if (oread.is_read1 and options.ForwardRead) or (not oread.is_read1 and options.ReverseRead):
            write_fastq(oread)
          elif (read.is_read1 and options.ForwardRead) or (not read.is_read1 and options.ReverseRead):
            write_fastq(read)
        del incomplete_pairs[read.qname]
      else:
        incomplete_pairs[read.qname] = read
    elif (not read.is_paired) and options.SingleReads and (options.all or (not read.is_qcfail)):
       write_fastq(read)
