#!/usr/bin/python

import sys, os
import scipy.stats
import pdb

# Performs Fisher exact test on gene burden counts for cases & controls, for both allele counts and sample counts.
# No decision is made as to the significance of the results. Only p-values are reported for later decision making.
#
# Two files are read in. One is the result of gene burden analysis performed on the cases.
# The other is the result of gene burden analysis performed on ethnic controls for the cases.
# The number of cases and controls must be specified in order to perform the statistical test.
# These files must be formatted in the exact same way. This shouldn't be a problem unless you have altered
# the program "exome_burden_script.py" at some point during your analysis pipeline.
#
# Each file is stored as a hash using gene name (first column) as the key. Only those genes appearing in the
# cases dataset are tested and reported.
#
# Output is written to stdout. These stdout statements are scattered throughout the code, instead of in a
# logical manner. Sorry.
#
#
### The class variable "nonadditiveCounts" requires some explaining. ###
# The original logic of the gene burden program was to add the counts from loss-of-function (lof) to
# deleterious. And then add this value to missense counts.
#
# By specifying "--individual-counts" this logic is reversed. So, that the value of lof is subtracted
# from deleterious and the value of deleterious is subtracted from missense. However, these two subtractions
# are done independently; not one after the other (otherwise, this could result in negative values for
# missense mutations).
#
# I added this logic in order to get the p-value for each mutation type on its own.
#
#
# gene_burden_fisher_exact.py -a <cases_file> -n <num-cases> -o <controls_file> -m <num-controls> [-i]


class FisherBurdenRight:

  casesFilepath = None
  numCases = None
  casesFH = None
  controlsFilepath = None
  numControls = None
  nonadditiveCounts = False
  controlsFH = None
  headerLine = None
  
  def __init__(self):
    self.parse_cmdline_parameters()
    self.verify_cmdline_parameters()
    return

  def __del__(self):
    if self.casesFH: self.casesFH.close()
    if self.controlsFH: self.controlsFH.close()

  def parse_cmdline_parameters(self):
    if 1 == len(sys.argv):
      self.print_short_help()
    index = 1
    while index < len(sys.argv):
      if '-a' == sys.argv[index] or '--cases-file' == sys.argv[index]:
        self.casesFilepath = sys.argv[index+1]
        index += 2
      elif '-n' == sys.argv[index] or '--num-cases' == sys.argv[index]:
        self.numCases = sys.argv[index+1]
        index += 2
      elif '-o' == sys.argv[index] or '--controls-file' == sys.argv[index]:
        self.controlsFilepath = sys.argv[index+1]
        index += 2
      elif '-m' == sys.argv[index] or '--num-controls' == sys.argv[index]:
        self.numControls = sys.argv[index+1]
        index += 2
      elif '-i' == sys.argv[index] or '--individual-counts' == sys.argv[index]:
        self.nonadditiveCounts = True
        index += 1
      else:
        sys.stderr.write('ERROR - invalid command line argument: {0}\n'.format(sys.argv[index]))
        self.print_short_help()
    return

  def print_short_help(self):
    sys.stderr.write('usage: {0} parameters\n'.format(sys.argv[0]))
    sys.stderr.write('Parameters:\n')
    sys.stderr.write('\t-a | --cases-file\tCASES gene burden output \n')
    sys.stderr.write('\t-n | --num-cases\tnumber of cases in gene burden output \n')
    sys.stderr.write('\t-o | --controls-file\tCONTROLS gene burden output \n')
    sys.stderr.write('\t-m | --num-controls\tnumber of controls in gene burden output \n')
    sys.stderr.write('Options:\n')
    sys.stderr.write('\t-i | --individual-counts\t FLAG: independent counts yield independent p-values\n')
    sys.stderr.flush()
    sys.exit(1)
    return

  def verify_cmdline_parameters(self):
    if None == self.casesFilepath or None == self.numCases or None == self.controlsFilepath or None == self.numControls:
      sys.stderr.write('ERROR - you must specify all four paramters\n')
      self.print_short_help()
    
    if not os.path.isfile(self.casesFilepath):
      sys.stderr.write('ERROR - unable to acces cases file: {0}\n'.format(self.casesFilepath))
      sys.exit(1)
    self.casesFH = open(self.casesFilepath, 'r')
    try:
      self.numCases = int(self.numCases)
    except ValueError as exc:
      sys.stderr.write('ERROR - you must supply an integer for the number of cases\n')
      sys.exit(1)
    
    if not os.path.isfile(self.controlsFilepath):
      sys.stderr.write('ERROR - unable to acces controls file: {0}\n'.format(self.controlsFilepath))
      sys.exit(1)
    self.controlsFH = open(self.controlsFilepath, 'r')
    try:
      self.numControls = int(self.numControls)
    except ValueError as exc:
      sys.stderr.write('ERROR - you must supply an integer for the number of controls\n')
      sys.exit(1)
    return

  def read_files(self):
    self.casesHash = self.hashify_file(self.casesFH)
    self.controlsHash = self.hashify_file(self.controlsFH)
    return

  def hashify_file(self, fh):
    fileHash = {}
    line = fh.readline().rstrip()
    if not self.headerLine:
      self.headerLine = line
    line = fh.readline().rstrip()
    while line:
      values = line.split('\t')
      key = values[0]
      data = values[1:]
      fileHash[key] = data
      line = fh.readline().rstrip()
    return fileHash

  def determine_p_value(self):
    self.write_header()
    controlGenes = sorted(self.controlsHash.keys())
    caseGenes = sorted(self.casesHash.keys())
    for geneName in caseGenes:
      sys.stdout.write(geneName)
      casesCounts = self.casesHash[geneName]
      controlsCounts = None
      if geneName in controlGenes:
        controlsCounts = self.controlsHash[geneName]
      else:
        controlsCounts = [0]*len(casesCounts)
      self.compare_arrays(casesCounts, controlsCounts)
    return

  def write_header(self):
    columns = self.headerLine.split('\t')
    sys.stdout.write(columns[0])
    index = 1
    while index < 13:
      sys.stdout.write('\t{0}(case)\t(control)\tp-value'.format(columns[index]))
      index += 1
    # set visual break to separate variant counts from sample counts
    sys.stdout.write('\tX')
    while index < len(columns):
      sys.stdout.write('\t{0}(case)\t(control)\tp-value'.format(columns[index]))
      index += 1
    sys.stdout.write('\n')
    return
    
  def compare_arrays(self, cases, controls):
    index = 0
    # test variant counts (two alleles at each position)
    while index < 12:
      casesMutated = int(cases[index])
      controlsMutated = int(controls[index])
      if self.nonadditiveCounts:
        # subtract LoF count to get individual counts for deleterious and missense variants
        if not (index + 1) % 3 == 0:
          casesMutated -= int(cases[index+1])
          controlsMutated -= int(controls[index+1])
      try:
        (oddsRatio, pValue) = scipy.stats.fisher_exact([[casesMutated, controlsMutated],
                                                        [self.numCases*2,  self.numControls*2]],
                                                       alternative='greater')
      except ValueError as exc:
        sys.stderr.write('ERROR - incorrect arguments sent to scipy.stats.fisher_exact()\n')
        sys.stderr.flush()
        sys.exit(1)
      sys.stdout.write('\t{0}\t{1}\t{2}'.format(cases[index], controls[index], pValue))
      index += 1

    sys.stdout.write('\tX')
      
    # test sample counts
    while index < 24:
      casesMutated = int(cases[index])
      controlsMutated = int(controls[index])
      (oddsRatio, pValue) = scipy.stats.fisher_exact([[casesMutated, controlsMutated],
                                                      [self.numCases,  self.numControls]],
                                                     alternative='greater')
      sys.stdout.write('\t{0}\t{1}\t{2}'.format(cases[index], controls[index], pValue))
      index += 1
    sys.stdout.write('\n')
    sys.stdout.flush()
    return


if __name__ == '__main__':
  x = FisherBurdenRight()
  x.read_files()
  x.determine_p_value()
