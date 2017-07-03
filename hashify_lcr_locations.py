#!/usr/bin/python

import sys, os
import cPickle
import pdb

# This program creates a hash from a tab-delimited file of low complexity regions.
# The structure of the input file is chrom, start stop.
# The input file is not expected to have a header line.
# The data in the file MUST be in ascending genomic order, with X, then Y at the end of the file.

# This program also skips as many header lines as specified.
# The hash is then stored in pickle format in the file specified.
# Values under keys are stored as a list, to allow for duplicate-key records to be stored without over-writing data.
# Arguments
#   -f | --file : input file
#   -o | --output : pickle file to create/overwrite

class HashifyLCR:

  filepath = None
  pickleFile = None
  dataHash = {}
  
  def __init__(self):
    self.parse_cmdline_parameters()
    self.verify_cmdline_parameters()
    return

  def parse_cmdline_parameters(self):
    if 1 == len(sys.argv):
      self.print_short_help()
    index = 1
    while index < len(sys.argv):
      try:
        if '-f' == sys.argv[index] or '--file' == sys.argv[index]:
          self.filepath = sys.argv[index+1]
          index += 2
        elif '-o' == sys.argv[index] or '--output' == sys.argv[index]:
          self.pickleFile = sys.argv[index+1]
          index += 2
        else:
          sys.stdout.write('ERROR - illegal parameter specified: {0}\n'.format(sys.argv[index]))
          self.print_short_help()
      except IndexError as exc:
        sys.stdout.write('ERROR - incorrect number of parameters specified\n')
        self.print_short_help()
    return

  def print_short_help(self):
    sys.stderr.write('usage: {0} parameters\n'.format(sys.argv[0]))
    sys.stderr.write('Parameters:\n')
    sys.stderr.write('\t-f | --file\t\tfile containing tab-separated values\n')
    sys.stderr.write('\t-o | --output\t\tname of pickle file to write to\n')
    sys.stderr.flush()
    sys.exit(1)
    return

  def verify_cmdline_parameters(self):
    if None ==  self.filepath or None == self.pickleFile:
      sys.stderr.write('ERROR - you must specify all parameters\n')
      self.print_short_help()
    if not os.path.isfile(self.filepath):
      sys.stderr.write('ERROR - unable to acces data file\n'.)
      self.print_short_help()
    dirName = os.path.dirname(self.pickleFile)
    if 0 < len(dirName):
      if not os.path.isdir(os.path.dirname(dirName)):
        sys.stderr.write('ERROR - unable to access directory containing \"output file\"\n'.)
        self.print_short_help()
    return

  def read_data(self):
    fh = open(self.filepath, 'r')
    line = fh.readline().rstrip()
    while line:
      (chrom, start, end) = line.split('\t')
      start = int(start)
      end = int(end)
      if chrom in self.dataHash.keys():
        self.dataHash[chrom].append([start, end])
      else:
        self.dataHash[chrom] = [[start, end]]
      line = fh.readline().rstrip()
    fh.close()
    return

  def write_pickle(self):
    fh = open(self.pickleFile, 'w')
    cPickle.dump(self.dataHash, fh)
    fh.close()
    return

if __name__ == '__main__':
  x = HashifyLCR()
  x.read_data()
  x.write_pickle()
