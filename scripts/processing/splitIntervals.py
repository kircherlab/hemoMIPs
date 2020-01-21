#!/bin/env python

import sys

for line in sys.stdin:
  chrom,posrange = line.strip().split(":")
  start,end = map(int,posrange.split("-"))
  
  if end-start > 300:
    csize = end-start
    number = (csize/200)+1
    nsize = csize/number
    
    hstart,hend = start-50,start+nsize+50
    print "%s:%d-%d"%(chrom,hstart,hend)
    while (hend < end):
      hstart,hend = hend-100,hend-50+nsize
      print "%s:%d-%d"%(chrom,hstart,hend)
  else:
    print "%s:%d-%d"%(chrom,start,end)
