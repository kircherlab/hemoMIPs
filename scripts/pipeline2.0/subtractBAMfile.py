#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *23.01.2016
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict 

import pysam
import gzip
from bx.intervals.intersection import Intersecter, Interval
import random

def lengthClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  sum = 0
  for (op,count) in cigar:
    if op in [4,5,6]: sum += count
  return sum

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def updateBuffer(chrom,start,end):
  global filter_file
  global bufferPos
  global options
  newBuffer = set()
  maxStart = start
  if options.veryVerbose: sys.stderr.write("LenBuffer-0: %d (%s:%d-%d)\n"%(len(bufferPos),chrom,start,end))
  # Check what we currently got buffered
  for (rstart,rend,valClipped) in bufferPos:
    if rend >= start:
      newBuffer.add((rstart,rend,valClipped))
    if rstart > maxStart: maxStart = rstart
  #if options.debug: del bufferPos
  bufferPos = newBuffer
  # Read new reads
  if options.veryVerbose: sys.stderr.write("LenBuffer-1: %d\n"%(len(bufferPos)))
  for read in filter_file.fetch(chrom,maxStart-1,end):
    if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
    valClipped = lengthClipped(read.cigar)
    
    if read.is_paired:
      if read.mate_is_unmapped: continue
      if read.rnext != read.tid: continue
      if read.isize == 0: continue
      rstart = min(read.pos,read.pnext)+1 # 1-based
      rend = rstart+abs(read.isize)-1 # end included
    else:
      rstart = read.pos+1 # 1-based
      rend = rstart+aln_length(read.cigar)-1 # end included
    
    bufferPos.add((rstart,rend,valClipped))
    if options.wiggle != 0:
      for wiggle in range(1,options.wiggle+1):
        bufferPos.add((rstart-wiggle,rend,valClipped))
        bufferPos.add((rstart+wiggle,rend,valClipped))
        bufferPos.add((rstart,rend-wiggle,valClipped))
        bufferPos.add((rstart,rend+wiggle,valClipped))
        bufferPos.add((rstart-wiggle,rend-wiggle,valClipped))
        bufferPos.add((rstart+wiggle,rend+wiggle,valClipped))
        bufferPos.add((rstart+wiggle,rend-wiggle,valClipped))
        bufferPos.add((rstart-wiggle,rend+wiggle,valClipped))
      
  if options.veryVerbose: sys.stderr.write("LenBuffer-2: %d\n"%(len(bufferPos)))

parser = OptionParser("%prog [options] inputFile subtractFile")
parser.add_option("-b","--buffer", dest="buffer", help="Buffer region +/- Xbp (default 5000)",default=5000,type="int")
parser.add_option("-u","--update", dest="update", help="Proportion of buffer passed before updating (default 0.9)",default=0.9,type="float")
parser.add_option("-w","--wiggle", dest="wiggle", help="Number of bases wiggle on the substracted reads (default 0)",default=0,type="int")
#parser.add_option("-d","--debug", dest="debug", help="Test run (default Off)",default=False,action="store_true")
parser.add_option("-r","--region", dest="region", help="Region to be looked up (def None)",default="")
parser.add_option("-R","--regionList", dest="regionList", help="List of region to be looked up (def None)",default="")
parser.add_option("-o","--outfile", dest="outfile", help="Output file (default output.bam)",default="output.bam")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("-V","--veryVerbose", dest="veryVerbose", help="Turn on very verbose debug output",default=False,action="store_true")
(options, args) = parser.parse_args()

if (options.update > 1) or (options.update < 0):
  sys.stderr.write("Warning: invalid update value, setting to 0.5\n")
  options.update = 0.5

if options.buffer < 1000:
  sys.stderr.write("Warning: invalid buffer value, setting to 1000\n")
  options.update = 1000
  
if len(args) < 2 or not os.path.exists(args[0]) or not os.path.exists(args[1]) or not os.path.exists(args[0]+".bai") or not os.path.exists(args[1]+".bai"):
  sys.stderr.write("Error: invalid input files or index files missing\n")
  sys.exit()

if options.veryVerbose: 
  options.verbose = True

chrom,start,end = None,None,None
outchrom = None
regionList = [(chrom,start,end)]
if options.region.strip() != "":
  try:
    chrom = options.region.split(':')[0]
    start,end = map(int,options.region.split(':')[1].replace(",","").split("-"))
    outchrom = chrom
    outchrom = outchrom.replace("chr","")
    if outchrom.startswith('gi|'): # gi|224589803|ref|NC_000012.11|
      NCfield = outchrom.split("|")[-2]
      if NCfield.startswith("NC"):
        outchrom = "%d"%(int(NCfield.split(".")[0].split("_")[-1]))
    if options.verbose: 
      sys.stderr.write("Coordinates parsed: Chrom %s Start %d End %d\n"%(chrom,start,end))
    regionList = [(chrom,start,end)]
  except:
    sys.stderr.write("Error: Region not defined in a correct format!\n")
    sys.exit()

if os.path.exists(options.regionList):
  if (regionList[0][0] == None): regionList = []
  regionFile = open(options.regionList)
  for line in regionFile:
    fields = line.rstrip().split(":")
    if len(fields) == 2:
      cchrom = fields[0].replace("chr","")
      fields = fields[1].split("-")
      if (len(fields) == 2) and fields[0].isdigit() and fields[1].isdigit():
        regionList.append((cchrom,int(fields[0]),int(fields[1])))
  regionFile.close()

if (chrom != None and start != None) or (options.regionList != "" and os.path.exists(options.regionList)):
  if not (os.path.exists(args[0].replace(".bam",".bai")) or os.path.exists(args[0]+".bai")):
    sys.stderr.write("Error: Coordinate ranges require indexed BAM files\n")
    sys.exit()

if options.verbose:
  if len(regionList) == 1 and regionList[0][0] == None:
    sys.stderr.write("Complete BAM file defined as input.\n"%(len(regionList)))
  else:
    sys.stderr.write("%d regions defined as input.\n"%(len(regionList)))

input_file = pysam.Samfile( args[0], "rb" )
output_file = pysam.Samfile( options.outfile, "wb", template=input_file)
filter_file = pysam.Samfile( args[1], "rb" )
bufferPos = set() 

ltid, lchrom = None,None
lstart, lend = None,None
lmiddle = None

countExcluded = 0
countKept = 0

iFile = input_file

for (chrom,start,end) in regionList:
  if chrom != None and start != None:
    try:
      iFile = input_file.fetch(chrom,start,end)
    except:
      continue
    if options.verbose: sys.stderr.write("Started new region: %s:%d-%d\n"%(chrom,start,end))
    bufferPos = set() 
    ltid, lchrom = None,None
    lstart, lend = None,None
    lmiddle = None
    
  for read in iFile:
    if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
    valClipped = lengthClipped(read.cigar)

    if ltid != read.tid:
      ltid = read.tid
      lchrom = input_file.getrname(ltid)
      lstart = max(1,read.pos-options.buffer//10)
      lend = read.pos+options.buffer
      lmiddle = round((lend-read.pos)*options.update)+read.pos
      #print lchrom,lstart,lend,lmiddle
      bufferPos = set()
      updateBuffer(lchrom,lstart,lend)
    if read.pos > lmiddle:
      lstart = max(1,read.pos-options.buffer//10)
      lend = read.pos+options.buffer
      lmiddle = round((lend-read.pos)*options.update)+read.pos
      updateBuffer(lchrom,lstart,lend)   
    
    if read.is_paired:
      if read.mate_is_unmapped: continue
      if read.rnext != read.tid: continue
      if read.isize == 0: continue
      rstart = min(read.pos,read.pnext)+1 # 1-based
      rend = rstart+abs(read.isize)-1 # end included
    else:
      rstart = read.pos+1 # 1-based
      rend = rstart+aln_length(read.cigar)-1 # end included

    if chrom != None and rstart < start: continue
    if chrom != None and rstart > end: continue

    if (rstart,rend,valClipped) not in bufferPos:
      output_file.write(read)
      countKept += 1
    else:
      countExcluded += 1

input_file.close()
output_file.close()
filter_file.close()

#if options.verbose:
if countExcluded+countKept > 0:
  sys.stderr.write("Finished. Kept: %d / Excluded: %d [%.2f%%]\n"%(countKept,countExcluded,countExcluded/float(countExcluded+countKept)*100))