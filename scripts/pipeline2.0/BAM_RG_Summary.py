#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *21.09.2011
"""

import sys, os
from optparse import OptionParser
import pysam

parser = OptionParser("%prog [options]")
parser.add_option("-a","--all", dest="all", help="Provide individual counts for all read groups",default=False,action="store_true")
parser.add_option("-q","--mapq", dest="mapq", help="Only consider reads with a map quality score above N (default 0)",default=0,type="int")
parser.add_option("-l","--minlength", dest="minlength", help="Only consider reads with a minimum length of N (default 0)",default=0,type="int")
parser.add_option("-c","--contig", dest="contig", help="Give separate counts for this contig (default '')",default="")
parser.add_option("--library", dest="library", help="Use library name from RG read header rather than the sample ID",default=False,action="store_true")
parser.add_option("--XI", dest="XI", help="Use XI field intead of RG",default=False,action="store_true")
parser.add_option("--XJ", dest="XJ", help="Use XJ field intead of RG",default=False,action="store_true")
(options, args) = parser.parse_args()

if options.library: options.all=True

have_XP = False
contig_tid = None
for filename in args:
  if os.path.exists(filename):
    print "Reading %s..."%filename
    rgroups = {}
    cbamfile = pysam.Samfile(filename, "rb" )
    id2lib = {}
    if options.library and 'RG' in cbamfile.header:
      for rgroup in cbamfile.header['RG']:
        if 'LB' in rgroup and 'ID' in rgroup:
          id2lib[rgroup['ID']] = rgroup['LB']
    for name in cbamfile.references:
      if name == options.contig: contig_tid = cbamfile.gettid(name)
    for read in cbamfile:
      if read.mapq >= options.mapq and len(read.seq) >= options.minlength:
        library,count = '',1
        for (key,value) in read.tags:
          if (not options.XI and not options.XJ and key == "RG") or (options.XJ and key == "XJ") or (options.XI and key == "XI"):
            if value in id2lib: library = id2lib[value]
            else: library = value
          elif key == "XP":
            have_XP = True
            count = value
        if library in rgroups:
          if not read.is_paired or read.is_read1:
            rgroups[library][0]+=1
            rgroups[library][1]+=count
          rgroups[library][2]+=len(read.seq)
          if read.tid == contig_tid: rgroups[library][3]+=1
        else:
          if not read.is_paired or read.is_read1:
            rgroups[library]=[1,count,len(read.seq),1 if read.tid == contig_tid else 0]
          else:
            rgroups[library]=[1,count,0,1 if read.tid == contig_tid else 0]
    cbamfile.close()
    tcount_reads,tcount_seq,tcounts_inc_dups,tcount_contig = 0,0,0,0
    header = False
    to_sort = rgroups.keys()
    to_sort.sort()
    for library in to_sort:
      counts = rgroups[library]
      if options.all:
        if have_XP and contig_tid != None: 
          if not header:
            header=True
            print "Library\tIndependentReads\tReads\tBases\tReadsOn:%s"%options.contig
          print "%s\t%d\t%d\t%d\t%d"%(library,counts[0],counts[1],counts[2],counts[3])
        elif have_XP and contig_tid == None: 
          if not header:
            header=True
            print "Library\tIndependentReads\tReads\tBases"
          print "%s\t%d\t%d\t%d"%(library,counts[0],counts[1],counts[2])
        elif not have_XP and contig_tid != None: 
          if not header:
            header=True
            print "Library\tReads\tBases\tReadsOn:%s"%options.contig
          print "%s\t%d\t%d\t%d"%(library,counts[0],counts[2],counts[3])
        elif not have_XP and contig_tid == None: 
          if not header:
            header=True
            print "Library\tReads\tBases"
          print "%s\t%d\t%d"%(library,counts[0],counts[2])
      tcount_reads += counts[0]
      tcounts_inc_dups += counts[1]
      tcount_seq += counts[2]
      tcount_contig += counts[3]
    if have_XP and contig_tid != None: print "Total Reads:\t%d\tReadsInclDups:\t%d\tBases:\t%d\tReads on %s:\t%d"%(tcount_reads,tcounts_inc_dups,tcount_seq,options.contig,tcount_contig)
    elif have_XP and contig_tid == None: print "Total Reads:\t%d\tReadsInclDups:\t%d\tBases:\t%d"%(tcount_reads,tcounts_inc_dups,tcount_seq)
    elif not have_XP and contig_tid != None: print "Total Reads:\t%d\tBases:\t%d\tReads on %s:\t%d"%(tcount_reads,tcount_seq,options.contig,tcount_contig)
    elif not have_XP and contig_tid == None: print "Total Reads:\t%d\tBases:\t%d"%(tcount_reads,tcount_seq)
