#!/usr/bin/env python
# -*- coding: ASCII -*-

import sys,os

count=0
cat=None
for line in sys.stdin:
  fields = line.split("\t")
  if len(fields) == 1: 
    fields = line.split()
  if len(fields) >= 2:
    ncat = "\t".join(fields[:-1])
    if ncat != cat:
      if count != 0:
        print "%s\t%d"%(cat,count)
      cat = ncat
      count = 0
    count+=int(fields[-1])
if count != 0:
  print "%s\t%d"%(cat,count)
