#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *06.05.2013
:Type: tool
:Input: BAM
:Output: BAM

"""

import sys,os
import math
import pysam
import string
from optparse import OptionParser,OptionGroup

from bx.intervals.intersection import IntervalTree, Interval

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

def identifyTrim(cigarlist,refbases,reverse=False):
  strand = 1
  if reverse: strand = -1
  new_cigar = []
  trim = 0
  keep = False
  for operation,length in cigarlist[::strand]:
    # Count bases that need to be cut off the read
    if not keep and (operation == 0 or operation == 1 or operation == 4 or operation >= 7):
      if operation == 1 or operation == 4:
        trim += length
        length = 0
      else:
        if length > refbases:
          trim += refbases
        else:
          trim += length
    # Count reference bases (to identify trim point)
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: 
      if length > refbases:
        keep = True
        length-=refbases
        refbases = 0
      else:
        refbases-=length
    if keep and length > 0 and operation != 5:
      new_cigar.append((operation,length))
  return trim*strand,new_cigar[::strand]

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def makeNewRead(trim,loff,newCigar,read,trim2=None,MIPid=None):
  a = pysam.AlignedRead()
  a.qname = read.qname

  a.pos = read.pos+loff
  if trim > 0:
    a.seq = read.seq[trim:trim2]
    a.qual = read.qual[trim:trim2]
  else:
    a.seq = read.seq[trim2:trim]
    a.qual = read.qual[trim2:trim]

  a.flag = read.flag
  # Convert PE to SR
  a.is_paired = False
  a.is_read1 = False
  a.is_read2 = False
  a.rname = read.rname
  a.mapq = read.mapq
  a.cigar = newCigar
  #a.mrnm = read.mrnm
  #a.mpos = read.mpos # WILL BE OFF
  #a.isize = read.isize # WILL BE OFF
  
  tags = []
  for key,value in read.tags:
    if key != "NM" and key != "MD" and key != "ZM":
      tags.append((key,value))
    elif key == "ZM" and MIPid == None:
      tags.append((key,value))
  if MIPid != None:
    tags.append(("ZM",MIPid))
  a.tags = tags
  return a


parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outfile", dest="outfile", help="Name of output file (def trimMIParms.bam)",default="trimMIParms.bam")
parser.add_option("-d", "--design", dest="design", help="MIP design file (default /net/shendure/vol1/home/mkircher/hemophilia/hemomips_design_file_updated.txt)",default="/net/shendure/vol1/home/mkircher/hemophilia/hemomips_design_file_updated.txt")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
(options, args) = parser.parse_args()


if options.verbose:
  sys.stderr.write("Evaluating design file.\n")

ligArms = {}
extArms = {}

if os.path.exists(options.design):
  # DESIGN FILE ARM RANGES ARE 1-BASED, CLOSED INTERVALS
  fchrom,ligStart,ligEnd,extStart,extEnd,fstrand,fmipname = 2,7,8,3,4,17,-1
  infile = open(options.design)
  for line in infile:
    fields = line.rstrip().split("\t")
    if len(fields) > 8:
      if line.startswith(">") or line.startswith('#'):
        for ind,elem in enumerate(fields):
          if elem == "chr" or elem == "chrom": fchrom = ind
          elif elem == "ext_probe_start": extStart = ind
          elif elem == "ext_probe_stop": extEnd = ind
          elif elem == "lig_probe_start": ligStart = ind
          elif elem == "lig_probe_stop": ligEnd = ind
          elif elem == "probe_strand": fstrand = ind
          elif elem == "mip_name": fmipname = ind
      else:
        chrom,lstart,lend,estart,eend,strand,mipname = fields[fchrom],int(fields[ligStart])-1,int(fields[ligEnd]),int(fields[extStart])-1,int(fields[extEnd]),fields[fstrand],fields[fmipname]
        if chrom not in ligArms: 
          ligArms[chrom] = IntervalTree()
          extArms[chrom] = IntervalTree()
        if strand == "+":
          ligArms[chrom].add_interval( Interval(lend, lend+1, value = (lend-lstart,mipname) ) )
          extArms[chrom].add_interval( Interval(estart+1, estart+2, value = (eend-estart,mipname) ) )
        else:
          ligArms[chrom].add_interval( Interval(lstart, lstart+1, value = (lend-lstart,mipname) ) )
          extArms[chrom].add_interval( Interval(eend, eend+1, value = (eend-estart,mipname) ) )
  infile.close()

else:
  sys.stderr.write("Design file (%s) not available.\n"%(options.design))
  sys.exit()

if options.verbose:
  sys.stderr.write("Opening input and output files.\n")

fileflags = 'wb'

if options.pipe: 
  infile = pysam.Samfile( "-", 'rb' )
  outfile = pysam.Samfile( "-", fileflags, template = infile)
else:
  if len(args) > 0 and os.path.exists(args[0]):
    infile = pysam.Samfile( args[0], 'rb' )
    outfile = pysam.Samfile( options.outfile , fileflags, template = infile)
  elif len(args) > 0:
    sys.stderr.write("Input file (%s) does not exist.\n"%(args[0]))
    sys.exit()
  else:
    sys.stderr.write("No input file specified.\n")
    sys.exit()
  
references = dict(map(lambda (x,y):(x,y),enumerate(infile.references)))

ltid,lchrom = None,None
count = 0
ltrim,ltrim2,loff,lnewcigar,lcigar,lpos,lstrand = None,None,None,None,None,None,None
ltchecked = None,None

for read in infile:
  cigar = read.cigar
  
  if cigar != None:
    alength = aln_length(cigar)
    start = read.pos
    end = read.pos+alength-1
    
    if read.tid == ltid:
      chrom = lchrom
    else:
      ltid = read.tid
      chrom = infile.getrname(ltid)
      lchrom = chrom

    if chrom in ligArms:
      found = False
      seq = read.seq
      qual = read.qual
      strand = read.is_reverse
      
      tchecked,mstart,mlength,mstart2,mlength2 = None,None,None,None,None
      mname,mname2 = None,None
      
      if read.is_paired:
        if (read.is_read1 and not strand): # EXPECT LIGATION ARM FOR READ START
          tchecked = "LigArmStart"
          for match in ligArms[chrom].find(start-1,start+2):
            mstart,mlength,mname = match.start,match.value[0],match.value[1]
            mlength += (mstart-start)
            found = True
            break
        elif (read.is_read1 and strand): # EXPECT LIGATION ARM FOR READ END
          tchecked = "LigArmEnd"
          for match in ligArms[chrom].find(end-1,end+2):
            mstart,mlength,mname = match.start,match.value[0],match.value[1]
            mlength -= (mstart-end)
            found = True
            break
        elif (read.is_read2 and strand): # EXPECT EXTENSION ARM FOR READ END
          tchecked = "ExtArmEnd"
          for match in extArms[chrom].find(end-1,end+2):
            mstart,mlength,mname = match.start,match.value[0],match.value[1]
            mlength -= (mstart-end)
            found = True
            break
        elif (read.is_read2 and not strand): # EXPECT EXTENSION ARM FOR READ START
          tchecked = "ExtArmStart"
          for match in extArms[chrom].find(start-1,start+2):
            mstart,mlength,mname = match.start,match.value[0],match.value[1]
            mlength += (mstart-start)
            found = True
            break
      else:
        if not strand: # EXPECT LIGATION ARM FOR READ START AND EXTENSION ARM FOR READ END
          tchecked = "ExtArmEnd"
          #sys.stderr.write('%d %d\n'%(start,end))
          for match in extArms[chrom].find(end-1,end+2):
            #sys.stderr.write('Found ExtArmEnd\n')
            mstart,mlength,mname = match.start,match.value[0],match.value[1]
            mlength -= (mstart-end)
            found = True
            break
          if found:
            tchecked += ";LigArmStart"
            found = False
            for match in ligArms[chrom].find(start-1,start+2):
              #sys.stderr.write('Found LigArmStart\n')
              mstart2,mlength2,mname2 = match.start,match.value[0],match.value[1]
              mlength2 += (mstart2-start)
              found = True
              break
          else:
            tchecked = None
            found = False
        else: # EXPECT LIGATION ARM FOR READ END AND EXTENSION ARM FOR READ START
          tchecked = "LigArmEnd"
          #sys.stderr.write('%d %d\n'%(start,end))
          for match in ligArms[chrom].find(end-1,end+2):
            mstart,mlength,mname = match.start,match.value[0],match.value[1]
            mlength -= (mstart-end)
            found = True
            break
          if found:
            tchecked += ";ExtArmStart"
            found = False
            for match in extArms[chrom].find(start-1,start+2):
              mstart2,mlength2,mname2 = match.start,match.value[0],match.value[1]
              mlength2 += (mstart2-start)
              found = True
              break
          else:
            found = False
        
      if found:
        mipID = None
        if mname == mname2 or mname2 == None: mipID = mname
        elif mname == None: mipID = mname2
        else:
          helper = [mname,mname2]
          helper.sort()
          mipID = "-".join(helper)

        if (tchecked == ltchecked) and (lpos == start) and (lchrom == chrom) and (lstrand == strand) and (lcigar == cigar):
          read = makeNewRead(ltrim,loff,lnewcigar,read,ltrim2,MIPid=mipID)
        else:
          if tchecked == "LigArmEnd;ExtArmStart" or "ExtArmEnd;LigArmStart":
            ltrim, lnewcigar = identifyTrim(cigar,mlength,reverse=True)
            ltrim2, lnewcigar = identifyTrim(lnewcigar,mlength2)
            loff = mlength2
          elif tchecked.endswith("End"):
            ltrim, lnewcigar = identifyTrim(cigar,mlength,reverse=True)
            loff = 0
          else:
            ltrim, lnewcigar = identifyTrim(cigar,mlength)
            loff = mlength
          read = makeNewRead(ltrim,loff,lnewcigar,read,ltrim2,MIPid=mipID)
          
          #if (lpos != start) or (lstrand != strand) or (lchrom != chrom):
            #sys.stderr.write("Found (%s): %s %d (%s) to %d %d\n"%(tchecked,chrom,start,"-" if strand else "+",mstart,mlength))

        outfile.write(read)
        count += 1
        #sys.stderr.write("%s %d-%d %d %s %s: %d %s %d | %s %s %s %s\n"%(chrom,start,end,strand,cigar,tchecked,ltrim,ltrim2,loff,mstart,mlength,mstart2,mlength2))
        lpos = start
        lstrand = strand
        lcigar = cigar
        ltchecked = tchecked
        #if count > 500: sys.exit()

      elif options.verbose:
        if (lpos != start) or (lstrand != strand) or (lchrom != chrom):
          sys.stderr.write("Read(s) not matched with design (%s): %s %d (%s)\n"%(tchecked,chrom,start,"-" if strand else "+"))
        lpos = start
        lstrand = strand
        lcigar = cigar
        tchecked = None
