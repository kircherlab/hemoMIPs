#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *17.07.2015

"""

import sys,os
from optparse import OptionParser
import pysam
import gzip
from collections import defaultdict

def sameSign(val1,val2):
  if val1 > 0 and val2 > 0: return True
  elif val1 < 0 and val2 < 0: return True
  else: return False

def sharedPrefix(s1,s2):
  minLength = min(len(s1),len(s2))
  shared = 0
  for ind in range(minLength-1):
    if s1[ind] == s2[ind]: shared+=1
    else: break
  if minLength == 1:
    return max(0,shared-1)
  else:
    return shared

def sharedSuffix(s1,s2):
  minLength = min(len(s1),len(s2))-1
  shared = 0
  for ind in range(minLength*-1,0)[::-1]:
    if s1[ind] == s2[ind]: shared+=1
    else: break
  return shared

parser = OptionParser()
parser.add_option("-b","--BAM",dest="BAM",help="BAM file",default="realign_all_samples.bam")
parser.add_option("-s","--sites",dest="sites",help="VCF file of variants (only InDels are extracted, def realign_all_samples.vcf.gz)",default="realign_all_samples.vcf.gz")
parser.add_option("-o","--outfile",dest="outfile",help="Write output to file (def STDOUT)",default='')
(options, args) = parser.parse_args()

outfile = sys.stdout if options.outfile == "" else open(options.outfile,'w')

if not os.path.exists(options.BAM) or not (os.path.exists(options.BAM+".bai") or os.path.exists(options.BAM.replace(".bam",".bai"))):
  sys.stderr.write("Error: input BAM file not found\n")
  sys.exit()
  
if not os.path.exists(options.sites):
  sys.stderr.write("Error: input VCF file not found\n")
  sys.exit()

variants = []
infile = gzip.open(options.sites)
for line in infile:
  if line.startswith("#"): continue
  else:
    fields = line.rstrip().split("\t")
    if len(fields[3]) != len(fields[4]):
      ref = fields[3]
      for alt in fields[4].upper().split(','):
        if len(ref) == len(alt): continue
        else:
          trimValue = sharedPrefix(ref,alt)
          if trimValue != 0:
            nref = ref[trimValue:]
            nalt = alt[trimValue:]
          else:
            nref,nalt = ref,alt
          trimValue2 = sharedSuffix(nref,nalt)
          if trimValue2 != 0:
            nref = nref[:-trimValue2]
            nalt = nalt[:-trimValue2]
          if len(nalt) == len(nref) and len(ref) == 1: continue
          elif (trimValue == 0):
            variants.append((fields[0],int(fields[1]),nref,nalt))
          else:
            variants.append((fields[0],int(fields[1])+trimValue,nref,nalt))
infile.close()

#print variants
#print variants[2:3]
#sys.exit()
#variants = [("X",154158201,"T","TG")]
#variants = [("X",138645010,"TATATATAATATATATATAAA","T")]
#variants = [("X",154158192,"C","CTTT")]

infile = pysam.Samfile(options.BAM, "rb" )
for chrom,pos,ref,alt in variants:
  event = len(alt)-len(ref)
  res = defaultdict(int)
  start = pos-1
  end = pos+max(len(ref),len(alt))-1
  for pileupcolumn in infile.pileup(chrom,start,end,max_depth=10**9):
    if start-5 <= pileupcolumn.pos <= end+5:
      for pileupread in pileupcolumn.pileups: 
        #print pileupread.indel,event,sameSign(pileupread.indel,event)
        if sameSign(pileupread.indel,event):
          RG = None
          for key,value in pileupread.alignment.tags:
            if key == "RG": RG=value
          res[RG]+=1
  #res = filter(lambda (x,y):y>=3,res.iteritems())
  res = map(lambda (x,y):"%d:%s"%(y,x),res.iteritems())
  outfile.write("%s_%d_%s/%s %s\n"%(chrom,pos,ref,alt," ".join(res)))
infile.close()

if options.outfile != "": outfile.close()