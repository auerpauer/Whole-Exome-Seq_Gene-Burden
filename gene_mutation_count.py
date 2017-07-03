#!/usr/bin/python

import sys, os
import pdb

# This program reads in two arguments. Each one is a comma-separated list.
# These lists are not allowed to have spaces due to bash command line processing.
#
# The first argument is a list of filepaths to be read in. Each of these files
# must be a variant table created by the python program "exome_burden_script.py".
#
# The second arguemnt is a list of headers that are to define the data pulled
# from the files. The order of the files and the headers must be in the same
# order. These headers appear at the tope of the tab-delimited file this
# program generates.
#
# For each file, the number of variants is counted for each gene.
# From this data, the program builds a tab-delimited file that gives the gene
# name and the variant count for each file.
#
# Output is written to stdout.

class GeneMutationCount:

  fileList = None
  annoList = None
  plDiffIndex = None
  geneIndices = None
  geneMutCounts = {}

  def __init__(self):
    self.parse_cmdline_parameters()
    self.verify_cmdline_parameters()
    return

  ### parse_cmdline_parameters ###
  def parse_cmdline_parameters(self):
    if 1 == len(sys.argv):
      self.print_short_help()
    index = 1
    while index < len(sys.argv):
      if '-f' == sys.argv[index] or '--file-list' == sys.argv[index]:
        self.fileList = sys.argv[index+1]
        index +=2
      elif '-a' == sys.argv[index] or '--annotation' == sys.argv[index]:
        self.annoList = sys.argv[index+1]
        index += 2
      else:
        sys.stdout.write('ERROR - illegal parameter specified\n')
        self.print_short_help()
    return

  ### print_short_help ###
  def print_short_help(self):
    sys.stderr.write('usage: {0} parameters\n'.format(sys.argv[0]))
    sys.stderr.write('Parameters:\n')
    sys.stderr.write('\t-f | --file-list\tcomma-separated list of variant files to read\n')
    sys.stderr.write('\t-a | --annotation\tcomma-separated list of names that appear at top of table\n')
    sys.stderr.write('\t\t\t\t names must match the order of files listed for \"file-list\"\n')
    sys.stderr.flush()
    sys.exit(1)
    return

  ### verify_cmdline_parameters ###
  def verify_cmdline_parameters(self):
    if None == self.fileList or None == self.annoList:
      sys.stderr.write('ERROR - you must specify all parameters\n')
      self.print_short_help()
    self.fileList = self.fileList.split(',')
    for elt in self.fileList:
      if not os.path.isfile(elt):
        sys.stderr.write('ERROR - unable to acces file: {0}\n'.format(elt))
        self.print_short_help()
    self.annoList = self.annoList.split(',')
    if len(self.annoList) != len(self.fileList):
      sys.stderr.write('ERROR - you must specify an annotation for each variant file\n')
      self.print_short_help()
    return

  ### process_files ###
  def process_files(self):
    index = 0
    while index < len(self.fileList):
      fh = open(self.fileList[index], 'r')
      line = fh.readline().rstrip() # read header line
      self.find_indices(line)
      self.read_variant_file(fh, index)
      fh.close()
      index += 1
    return

  ### find_indices ###
  def find_indices(self, line):
    headerColumns = line.split('\t')
    self.plDiffIndex = headerColumns.index('PLDIFF/DEPTH')
    # There are two columns for gene name, one for each allele
    self.geneIndices = [i for i, x in enumerate(headerColumns) if x == 'Gene.refGene']
    return

  ### read_variant_file ###
  # Check lines for usable data. When found, save the gene name.
  # Check if this gene name already exists in the hash 'geneMutCounts'.
  # 'fileIndex' is the index of the current file being read from the 'fileList'.
  # 'fileList' holds the paths to the variant tables this program reads to get the gene counts.
  # When a new gene name is found, we add it to 'geneMutCounts' along with a number of
  # empty lists equal to the number of variant files specified by the user
  def read_variant_file(self, fh, fileIndex):
    line = fh.readline().rstrip()
    while line:
      values = line.split('\t')
      plDiff = values[self.plDiffIndex]
      if 'lcr_region' == plDiff: # skip low-complexity regions (LCR)
        line = fh.readline().rstrip()
        continue
      geneName = self.determine_current_gene(values)
      if '.' == geneName:
        sys.stdout.write('ERROR - unable to determine gene name for current line:\n{0}\n'.format(line))
        sys.exit(1)
      if geneName not in self.geneMutCounts.keys():
        # create a new entry for this gene with all columns set to zero
        self.geneMutCounts[geneName] = [0]*len(self.fileList)
      self.geneMutCounts[geneName][fileIndex] += 1 # increment the correct column for this gene
      line = fh.readline().rstrip()
    return

  ### determine_current_gene ###
  def determine_current_gene(self, values):
    if '.' == values[self.geneIndices[0]]:
      return values[self.geneIndices[1]]
    return values[self.geneIndices[0]]

  ### write_results ###
  def write_results(self):
    strList = [0]*len(self.fileList)
    sys.stdout.write( 'Gene\t{0}\n'.format('\t'.join(self.annoList)) )
    geneList = sorted(self.geneMutCounts.keys())
    # in order to join the list of values, they must be strings
    for gene in geneList:
      strList = [str(x) for x in self.geneMutCounts[gene]]
      sys.stdout.write( '{0}\t{1}\n'.format(gene, '\t'.join(strList)) )
    return

if __name__ == '__main__':
  x = GeneMutationCount()
  x.process_files()
  x.write_results()
