#!/usr/bin/python

from __future__ import division
import sys, os
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# UNNECESSARY IMPORTS DUE TO USE OF JUPYTER NOTEBOOK
#import numpy as np
#import matplotlib
#matplotlib.use('Agg') # b/c "display" is not set
#import matplotlib.pyplot as plt
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
import cPickle
import pdb

# This program generates hashes
#   Key: sample name
#   Value: number of times this sample appears in the variant file
# In other words, this program counts the number SNVs per sample.
# This is then stored in a pickle file.
# These pickle files are then read in to create histograms in matplotlib, using a jupyter (ipython) notebook.
# These files are also used by the program wes_spreadsheet_generator.py
#
# The number of times the sample occurs in the file is also the number of variants it contributes to the data.
# This count can be limited by only counting variants in a given set of significant genes.

class CountSampleVariants:

  plDiffCutoff = 8 # 7 loose; 8 stringent

  caseFilePath = None
  caseFH = None
  dirName = None
  outPickle = None
  signifGenesPath = None
  signifGenes = None
  mutCount = {}
  
  def __init__(self):
    self.parse_cmdline_parameters()
    self.verify_cmdline_parameters()
    return

  def __del__(self):
    if self.caseFH: self.caseFH.close()
    #if self.ctrlFH: self.ctrlFH.close()

  ### parse_cmdline_parameters ###
  def parse_cmdline_parameters(self):
    if 1 == len(sys.argv):
      self.print_short_help()
    index = 1
    while index < len(sys.argv):
      if '-v' == sys.argv[index] or '--variant-file' == sys.argv[index]:
        self.caseFilePath = sys.argv[index+1]
        index += 2
      elif '-o' == sys.argv[index] or '--output-file' == sys.argv[index]:
        self.outPickle = sys.argv[index+1]
        index += 2
      elif '-s' == sys.argv[index] or '--significant-genes' == sys.argv[index]:
        self.signifGenesPath = sys.argv[index+1]
        index += 2
      else:
        sys.stdout.write('ERROR - illegal parameter specified\n')
        self.print_short_help()
    return

  ### print_short_help ###
  def print_short_help(self):
    sys.stderr.write('usage: {0} parameters\n'.format(sys.argv[0]))
    sys.stderr.write('Parameters:\n')
    sys.stderr.write('\t-v | --variant-file\tvariant table file generated from \"exome_burden_script.py\"\n')
    sys.stderr.write('\t-o | --output-file\tfull path of file to new pickle file\n')
    sys.stderr.write('\t\t\t\t  existing files will be overwritten\n')
    sys.stderr.write('\t-s | --significant-genes\tfile containing HUGO names of significant genes to analyze\n')
    sys.stderr.write('\t\t\t\t\t  lines starting with \"#\" will be ignored\n')
    sys.stderr.write('\t\t\t\t\t  use \"all\" to count all genes\n')
    sys.stderr.flush()
    sys.exit(1)
    return

  ### verify_cmdline_parameters ###
  def verify_cmdline_parameters(self):
    # verify and open variant files
    if None == self.caseFilePath or None == self.outPickle or None == self.signifGenesPath:
      sys.stderr.write('ERROR - you must specify all parameters\n')
      self.print_short_help()
    if not os.path.isfile(self.caseFilePath):
      sys.stderr.write('ERROR - unable to access case file: {0}\n'.format(self.caseFilePath))
      sys.exit(1)
    self.caseFH = open(self.caseFilePath, 'r')

    # verify other files or parent directories of files
    self.dirName = os.path.dirname(self.outPickle)
    if 0 < len(self.dirName):
      pickleDir = os.path.dirname(self.outPickle)
      if 0 < len(pickleDir):
        if not os.path.exists(pickleDir):
          sys.stderr.write('ERROR - unable to access parent directory of output pickle file: {0}\n'.format(pickleDir))
          self.print_short_help()
    # if the output file does not end in ".pickle", add it
    pickleExtension = os.path.basename(self.outPickle).split('.')[0]
    if 'pickle' != pickleExtension.lower():
      pickleExtension = '{0}.pickle'.format(pickleExtension)
    if ( 'all' != self.signifGenesPath ) and  ( not os.path.isfile(self.signifGenesPath) ):
      sys.stderr.write('ERROR - unable to access parent directory of significant genes file: {0}\n'.format(self.signifGenesPath))
      self.print_short_help()
    return

  def read_signif_gene_file(self):
    if 'all' == self.signifGenesPath:
      return
    self.signifGenes = []
    fh = open(self.signifGenesPath, 'r')
    line = fh.readline().rstrip()
    while line:
      if not '#' == line[0]:
        self.signifGenes.append(line)
      line = fh.readline().rstrip()
    fh.close()
    # if there are no genes read in, count hits for all genes
    if 0 == len(self.signifGenes):
      self.signifGenesPath = 'all'
    return
    
  ### read_case_variants ###
  def read_case_variants(self):
    headerData = self.caseFH.readline().rstrip().split('\t')
    plDiffIndex = headerData.index('PLDIFF/DEPTH')
    geneIndices = [i for i,x in enumerate(headerData) if x == 'Gene.refGene']
    line = self.caseFH.readline().rstrip()
    while line:
      values = line.split('\t')
      if not self.plDiff_score_pass(values, plDiffIndex):
        line = self.caseFH.readline().rstrip()
        continue
      geneName = self.determine_gene_name(geneIndices, values)
      # filter out non-significant genes, if any were specified
      if 'all' != self.signifGenesPath and geneName not in self.signifGenes:
        line = self.caseFH.readline().rstrip()
        continue
        
      sampleName = values[0]
      if sampleName in self.mutCount.keys():
        self.mutCount[sampleName] += 1
      else:
        self.mutCount[sampleName] = 1
      line = self.caseFH.readline().rstrip()
    fh = open(self.outPickle, 'w')
    cPickle.dump(self.mutCount, fh)
    self.caseFH.close()
    return

  ### plDiff_score_pass ###
  def plDiff_score_pass(self, values, plDiffIndex):
    plScore = values[plDiffIndex]
    if 'lcr_region' == plScore:
      return False
    if 'hom' != plScore:
      plScore = int(plScore)
      if self.plDiffCutoff > plScore:
        return False
    return True

  def determine_gene_name(self, geneIndices, values):
    name = values[geneIndices[0]]
    if '.' == name:
      name = values[geneIndices[1]]
    return name

  ### plot_data ###
  # NO LONGER NECESSARY DUE TO USE OF JUPYTER NOTEBOOK
  #def plot_data(self):
  #  graphFilename = os.path.join(self.dirName, 'EC_count_histogram.pdf')
  #  countArray = np.array(sorted(self.mutCount.values()))
  #  ind = np.arange(len(countArray))
  #  width = 1/len(countArray) + 0.1
  #  plt.bar(ind, countArray, width, align='center')
  #  plt.title('Mutations per Elite Controller')
  #  plt.ylabel('Count')
  #  plt.savefig(graphFilename)
  #  return


if __name__ == '__main__':
  x = CountSampleVariants()
  x.read_signif_gene_file()
  x.read_case_variants()
  #x.plot_data()
