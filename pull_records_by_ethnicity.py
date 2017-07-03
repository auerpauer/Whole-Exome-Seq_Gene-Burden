#!/usr/bin/python

#from optparse import OptionParser
import sys, os
import pdb

# This script copies the meta-data, header line and variant records of the input VCF file.
# However, only the columns containing samples specified in the ethnicity file are copied over.
# Once the new record containing only the samples of interest is built, the record is checked
# to ensure at least one of these samples still has a useful genotype.
# This is done by filtering against the class variable 'unwantedGenotypes'.
# The header line is changed to include only the sample names that were copied over.
# The input ethnicity file must list one sample name per line, with no header lines.
#
# All normal output is sent to stdout.
#
# Note, this script can be used to copy any subset of samples to a new file.
#
# This program can also work on a file that it previously created.
# This means you can first separate out samples based on a filter.
# Then you can call this program on that output to further separate out by a more stringent filter.
#
# PLEASE NOTE: The file containing sample names is called ethnicity files due to the need for my projects
# to separate samples by ethnicity. The ethnicity of samples IS NOT CHECKED by this program.

class PullSampleSubset:

  unwantedGenotypes = ['0/0', '0|0', '.\.', '.|.']
  
  vcfFilePath = None
  subsetFilePath = None
  desiredSamples = {}
  missingSamples = {}
  vcfFH = None
  ethnicityFH = None
  
  def __init__(self):
    self.parse_cmdline_parameters()
    self.verify_cmdline_parameters()
    return

  def __del__(self):
    if self.vcfFH: self.vcfFH.close()
    if self.ethnicityFH: self.ethnicityFH.close()
    return

  def parse_cmdline_parameters(self):
    if 1 == len(sys.argv):
      self.print_short_help()
    index = 1
    while index < len(sys.argv):
      try:
        if '-v' == sys.argv[index] or '--vcf-file' == sys.argv[index]:
          self.vcfFilePath = sys.argv[index+1]
          index += 2
        elif '-e' == sys.argv[index] or '--ethnicity-file' == sys.argv[index]:
          self.subsetFilePath = sys.argv[index+1]
          index += 2
        else:
          sys.stderr.write('ERROR - illegal parameter specified\n')
          self.print_short_help()
      except IndexError as exc:
        sys.stderr.write('ERROR - incorrect number of arguments\n')
        self.print_short_help()
    return
  
  def print_short_help(self):
    sys.stderr.write('usage: {0} parameters\n\n'.format(os.path.basename(sys.argv[0])))
    sys.stderr.write('Parameters:\n')
    sys.stderr.write('\t-v | --vcf\t\tVCF file to read in\n')
    sys.stderr.write('\t-e | --ethnicity-file\tfile containing samples to retain\n')
    sys.stderr.flush()
    sys.exit(1)
    return

  def verify_cmdline_parameters(self):
    if None == self.vcfFilePath or None == self.subsetFilePath:
      sys.stdout.write('ERROR - you must specify all parameters\n')
      self.print_short_help()
    if not os.path.isfile(self.vcfFilePath):
      sys.stderr.write('ERROR - VCF file is not accessible: {0}\n'.format(self.vcfFilePath))
      self.print_short_help()
    self.vcfFH = open(self.vcfFilePath, 'r')
    if not os.path.exists(self.subsetFilePath):
      sys.stderr.write('ERROR - ethnicity file is not accessible: {0}\n'.format(self.subsetFilePath))
      self.print_short_help()
    return

  ### read_sample_names ###
  def read_sample_names(self):
    self.ethnicityFH = open(self.subsetFilePath, 'r')
    line = self.ethnicityFH.readline().rstrip()
    while line:
      self.desiredSamples[line.upper()] = None
      line = self.ethnicityFH.readline().rstrip()
    self.ethnicityFH.close()
    self.ethnicityFH = None
    return

  ### create_new_vcf_file ###
  def create_new_vcf_file(self):
    headerLine = self.skip_vcf_meta_data()
    newHeaderLine = self.save_sample_positions(headerLine)
    sys.stdout.write('{0}\n'.format(newHeaderLine))
    line = self.vcfFH.readline().rstrip()
    while line:
      newRecord = self.build_new_record(line)
      if self.verify_nonempty_record(newRecord):
        print '\t'.join(newRecord)
      line = self.vcfFH.readline().rstrip()
    self.vcfFH.close()
    self.vcfFH = None
    return

  ### skip_vcf_meta_data ###
  # called by create_new_vcf_file()
  # Echo all records to stdout until the header line is found.
  # This way, we maintain all the meta-data of the VCF file.
  # This may include data that is no longer relevant, but this is unlikely.
  def skip_vcf_meta_data(self):
    headerLine = None
    line = self.vcfFH.readline().rstrip()
    while '#CHROM' != line[0:6]:
      sys.stdout.write('{0}\n'.format(line))
      line = self.vcfFH.readline().rstrip()
    return line
      
  ### save_sample_positions ###
  # called by create_new_vcf_file()
  # Save the index of each sample name in the header.
  # This index also points to the data for each sample in the records that follow.
  def save_sample_positions(self, headerLine):
    headerColumns = headerLine.split('\t')
    genomicColumns = headerColumns[0:9]
    index = 0
    sampleNames = sorted(self.desiredSamples.keys()) # is it necessary to sort the keys?
    for name in sampleNames:
      self.desiredSamples[name] = headerColumns.index(name)
    genomicColumns.extend(sampleNames)
    newHeader = '\t'.join(genomicColumns)
    return newHeader

  ### build_new_record ###
  # called by create_new_vcf_file()
  def build_new_record(self, line):
    columns = line.split('\t')
    newRecord = columns[0:9]
    for name in sorted(self.desiredSamples.keys()):
      sampleIndex = self.desiredSamples[name]
      if None == sampleIndex:
        sys.stderr.write('ERROR - requested sample does not exist in VCF file: {0}\n'.format(name))
        self.print_short_help()
      newRecord.append(columns[sampleIndex])
    return newRecord

  ### verify_nonempty_record ###
  # After removing all unwanted samples from this record, we need to verify that at least one
  # of these samples contains at least one alternate allele at this genomic locations.
  # Upon seeing the first sample that contains an alterante allele at this position, 'True' is returned
  # Otherwise, if all samples are checked and none of them contain an alternate allele, 'False' is returned.
  def verify_nonempty_record(self, record):
    sampleData = record[9:]
    for sample in sampleData:
      genotype = sample.split(':')[0]
      if not genotype in self.unwantedGenotypes:
        return True
    return False


if __name__ == '__main__':
  x = PullSampleSubset()
  x.read_sample_names()
  x.create_new_vcf_file()
