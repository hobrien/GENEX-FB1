#!/usr/bin/env python

import sys
import subprocess
import fileinput
import warnings
from string import maketrans

"""
Parse GTF from Tophat and add to DB
""" 

def main(): 
    print '\t'.join(['seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'gene_type', 'gene_status', 'gene_name'])       
    for line in fileinput.input([]):
       line = line.strip()
       try:
           parsed = parse_GTF(line.split('\t'))
       except IndexError:
           warnings.warn("%s not a correclty formatted GTF line" % line)
           continue
       if parsed['feature'] == 'gene':
           print '\t'.join([parsed['seqid'], 
                           parsed['source'], 
                           parsed['feature'], 
                           str(parsed['start']), 
                           str(parsed['end']), 
                           str(parsed['score']), 
                           str(parsed['strand']),
                           str(parsed['frame']),
                           parsed['gene_id'],
                           parsed['gene_type'],
                           parsed['gene_status'],
                           parsed['gene_name']])
                            
    
    
    
def parse_GTF (fields):
  tags = {}
  for attribute in fields[8].split(";")[:-1]:
    attribute = attribute.strip()
    tags[attribute.split(" ")[0]] = " ".join(attribute.split(" ")[1:]).replace('"','')
  try:
    tags['frame'] = int(fields[7])
  except ValueError:
    if fields[7] == '.':
      tags['frame'] = ''
    else:
      sys.exit("frame %s not recognized. Must be 1, 2, 3 or ." % fields[5]) 
  if fields[6] == '-' or fields[6] == 0 or fields[6] == -1:
    tags['strand'] = 0
  elif fields[6] == '+' or fields[6] == 1:
    tags['strand'] = 1
  elif fields[6] == '.':
    tags['strand'] = ''
  else:  
    sys.exit("strand %s not recognized. Must one of +, -, 1, -1 or 0" % fields[6])      
  try:
    tags['score'] = float(fields[5])
  except ValueError:
    if fields[5] == '.':
      tags['score'] = ''
    else:
      sys.exit("score %s not recognized. Must be a number" % fields[5])
  try:
    tags['end'] = int(fields[4])
  except ValueError:
    sys.exit("score %s not recognized. Must be a positive integer" % fields[4])
  try:
    tags['start'] = int(fields[3])
  except ValueError:
    sys.exit("score %s not recognized. Must be a positive integer" % fields[3])
  tags['feature'] = fields[2]
  tags['source'] = fields[1]
  tags['seqid'] = fields[0]
  return tags
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)
           
if __name__ == "__main__":
    warnings.formatwarning = warning_on_one_line
    main()
