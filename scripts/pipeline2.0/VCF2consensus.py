#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *10.04.2014
"""

import sys, os
import pysam

import random
from optparse import OptionParser
from collections import defaultdict

###################################
# Example uses:
#
# samtools mpileup -F 0.05 -DBgu -f reference.fa bamfile.bam | bcftools view -cg - | \
# ./VCF2consensus.py -n consensus -e 16458 > outfile.fa
#
# samtools mpileup -F 0.05 -DBgu -f reference.fa bamfile.bam | bcftools view -cg - | \
# ./VCF2consensus.py -n consensus -e 16458 -c 3 -s 0.66 > outfile.fa
#
# samtools mpileup -F 0.05 -DBgu -f reference.fa bamfile.bam | bcftools view -cg - | \
# ./VCF2consensus.py -n consensus -e 16458 -c 10 -s 0.9 > outfile.fa
#
###################################

# VCF FIELDS
fchr = 0
fpos = 1
fdbsnp = 2
fref_allele = 3
falt_allele = 4
fgeno_qual = 5
fflag = 6
finfo = 7
fformat = 8
fvalues = 9

parser = OptionParser()
parser.add_option("-s", "--minAlleleSupport", dest="minAlleleSupport", help="Minimum requiered support for an allele to be not printed as N (default 0, e.g. 0.66|0.9)",default=0,type="float")
parser.add_option("-i", "--minAlleleSupportInDel", dest="minAlleleSupportInDel", help="Minimum requiered support for an allele to be not printed as N (default 0, e.g. 0.66|0.9)",default=0,type="float")
parser.add_option("-c", "--minCoverage", dest="minCoverage", help="Minimum requierd coverage before a site is considered (default 0)",default=0,type="int")
parser.add_option("-n", "--name", dest="name", help="Name of consensus sequence for fasta output (default consensus)",default="consensus")
parser.add_option("-e","--endFill", dest="endFill", help="Fill up end of sequence if last contig position is shorter than N bases (default 0)",default=0,type="int")
parser.add_option("--cutoff_QUAL", dest="QUALcutoff", help="VCF QUAL cutoff (default 30)",default=30,type="float")
parser.add_option("--cutoff_coverage", dest="CoverageCutoff", help="Minimum coverage cutoff (default 8)",default=8,type="float")
(options, args) = parser.parse_args()

if options.minAlleleSupport < 0:
  options.minAlleleSupport = 0
elif options.minAlleleSupport > 1:
  options.minAlleleSupport = 1

if options.minAlleleSupportInDel < 0:
  options.minAlleleSupportInDel = 0
elif options.minAlleleSupportInDel > 1:
  options.minAlleleSupportInDel = 1

header = None
skip = 0
consensus = ""
chrom,lpos = None,0
pdepth,depth = None,None
for line in sys.stdin:
  if line.startswith("##"): continue # VCF specification lines
  elif line.startswith("#"): # HEADER
    header = line[1:].rstrip().split("\t")
  else:
    fields = line.rstrip().split("\t")
    if chrom == None: chrom = fields[fchr]
    if chrom != fields[fchr]:
      sys.stderr.write("Error: Can not handle multiple contigs!\n")
      break

    #N gap filling
    pos = int(fields[fpos])
    if lpos+1 < pos:
      pdepth = None
      fillN = pos-(lpos+1)
      if skip > fillN: 
        skip -= fillN
        fillN = 0
      elif fillN > skip: 
        fillN -= skip
        skip = 0
      if fillN > 0:
        consensus += fillN*'N'
        #sys.stderr.write("Info: Filling %d Ns before site %s:%s\n"%(fillN,fields[fchr],fields[fpos]))
    lpos=pos

    ref,alts = fields[fref_allele],fields[falt_allele]
    alleles = [ref]+alts.split(',')
    info = dict(map(lambda x: tuple(x.split('=')) if ('=' in x) else (x,True), fields[finfo].split(';')))
    values = dict(zip(fields[fformat].split(":"),fields[fvalues].split(":")))
    
    pdepth=depth
    depth = 0.0
    if 'DP' in values: # Use filtered read depth from sample
      depth = float(values['DP'])
    elif 'DP' in info: # revert to raw read depth
      depth = float(info['DP'])
      
    if "INDEL" in info and alts == ".": continue
    elif "INDEL" in info:
      if pdepth != None: depth = pdepth
      if depth < options.minCoverage: continue
      
      if "," in alts: 
        if ('GT' in values) and (values['GT'] in ["1/1","0/1","0/0"]): 
          alts = alts.split(",")[0]
          alleles = [ref]+alts.split(',')
        else:
          sys.stderr.write("Error: Can not handle InDel with two alternative calls (Site %s:%s)!"%(fields[fchr],fields[fpos]))
          sys.exit()

      countRef,countAlt=0,0
      if ('DP4' in info):
        counts = map(int,info['DP4'].split(","))
        countRef,countAlt = counts[0]+counts[1], counts[2]+counts[3]
        #sys.stderr.write("InDel (%d): %d\t%d\t%d\t%.4f\t%.4f\t%.4f\n"%(pos,countRef,countAlt,depth,countRef/depth,countAlt/depth,(depth-countRef-countAlt)/depth))
      else:
        if options.minAlleleSupportInDel > 0:
          sys.stderr.write("Warning: Site (%s:%s) does not have a DP4 field. Can not do minimum support filtering!\n"%(fields[fchr],fields[fpos]))
          
      if len(ref) <= len(alts): # INSERTION
        insertedseq = alts.replace(ref,"",1)
        #totalseq = alts[1:]
        isOK = False
        if ('DP4' in info):
          if (countRef < countAlt) and (countAlt/depth >= options.minAlleleSupportInDel): isOK = True
        elif ('GT' in values) and (values['GT'] == "1/1"): isOK = True
        if isOK: 
          #consensus+=totalseq
          #skip += (len(totalseq)-len(insertedseq))
          consensus+=insertedseq
          #sys.stderr.write("INFO: INSERTION ADDED!\n")
        #else:
          #sys.stderr.write("INFO: INSERTION DID NOT PASS!\n")
        continue
      else: # DELETION
        deleted = ref.replace(alts,"",1)
        isOK = False
        if ('DP4' in info):
          if (countRef < countAlt) and (countAlt/depth >= options.minAlleleSupportInDel): isOK = True
        elif ('GT' in values) and (values['GT'] == "1/1"): isOK = True
        if isOK: 
          skip += len(deleted)
          #sys.stderr.write("INFO: DELETION IN CONSENSUS!\n")
        #else:
          #sys.stderr.write("INFO: DELETION DID NOT PASS!\n")
        continue
    else:
      if skip > 0: 
        skip-=1
        continue
      
      if depth < options.minCoverage or depth == 0: 
        consensus+='N'
        continue
      
      countRef,countAlt=0,0
      if ('DP4' in info):
        counts = map(int,info['DP4'].split(","))
        countRef,countAlt = counts[0]+counts[1], counts[2]+counts[3]
        #sys.stderr.write("SNV (%d): %d\t%d\t%d\t%.4f\t%.4f\t%.4f\n"%(pos,countRef,countAlt,depth,countRef/depth,countAlt/depth,(depth-countRef-countAlt)/depth))
      else:
        if options.minAlleleSupport > 0:
          sys.stderr.write("Warning: Site (%s:%s) does not have a DP4 field. Can not do minimum support filtering!\n"%(fields[fchr],fields[fpos]))

      if alts == ".":
        if ('DP4' not in info) or (countRef/depth >= options.minAlleleSupport): consensus += ref
        else: consensus+='N'
      else:
        if "," in alts: 
          if ('GT' in values) and (values['GT'] in ["1/1","0/1","0/0"]): 
            alts = alts.split(",")[0]
            alleles = [ref]+alts.split(',')
          else:
            sys.stderr.write("Warning: Site (%s:%s) with multiple alleles -- allele count of second alternative allele is not available!\n"%(fields[fchr],fields[fpos]))
        if ('GT' in values) and len(set(values['GT'].split("/"))) == 1: # HOMOZYGOTE CALL
            alleleID = int(values['GT'].split("/")[0])
            if alleleID == 0:
              if ('DP4' not in info) or (countRef/depth >= options.minAlleleSupport): consensus += alleles[alleleID]
              else: consensus+='N'
            elif alleleID == 1:
              if ('DP4' not in info) or (countAlt/depth >= options.minAlleleSupport): consensus += alleles[alleleID]
              else: consensus+='N'
            else:
              if ('DP4' not in info) or ((depth-countRef-countAlt)/depth >= options.minAlleleSupport): consensus += alleles[alleleID]
              else: consensus+='N'
        elif 'DP4' in info:
          if countRef > countAlt: 
            if countRef/depth > options.minAlleleSupport: consensus+=ref
            else: consensus+='N'
          elif countRef < countAlt: 
            if countAlt/depth > options.minAlleleSupport: consensus+=alleles[1]
            else: consensus+='N'
          else:
            if options.minAlleleSupport > 0.5: consensus+='N'
            else:
              if ref == "C" and alleles[1] == "T": consensus+="C"
              elif ref == "T" and alleles[1] == "C": consensus+="C"
              elif ref == "G" and alleles[1] == "A": consensus+="G"
              elif ref == "A" and alleles[1] == "G": consensus+="G"
              else: consensus+=random.choice(alleles[:2])
        else:
          sys.stderr.write("Warning: Site (%s:%s) with no GT or DP4 information available. Reporting reference allele.\n"%(fields[fchr],fields[fpos]))
          consensus+=ref

if lpos+1 < options.endFill:
  fillN = options.endFill-(lpos+1)
  if skip > fillN: 
    skip -= fillN
    fillN = 0
  elif fillN > skip: 
    fillN -= skip
    skip = 0
  consensus += fillN*'N'

if len(consensus) > 0:
  print ">%s"%(options.name)
  while len(consensus) > 0:
    print consensus[:60]
    consensus = consensus[60:]
