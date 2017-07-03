#!/usr/bin/python

#from optparse import OptionParser
import sys, os, re
import traceback
import cPickle
import pdb
  
# INTERNAL EXCEPTION CLASS
# This is used to pass error message up the stack to a context that has more info to understand and report the error.
# The "stack" holds the methods that raise and pass this object up the stack.
# New methods are appended to the end, but when reporting the stack, this order is reveresed, to give highest level method first.
class BurdenException(Exception):
  # error codes
  input_error = 1010
  parsing_error = 2010
  filtered_out = 3000 # All error codes larger in value than this one represent a record that is filtered out, but not erroneous
  filter_freq_too_high = 3010 # All 301X error codes represent a variant that occurs too often to be of interests to us
  filter_freq_gnomad_exome_too_high = 3011
  filter_freq_gnomad_genome_too_high = 3012
  filter_freq_exac_too_high = 3013
  filter_on_refgene = 3020 # represents a record that did not pass filters on refGene values
  filter_nonexonic_nonsplicing = 3021
  filter_nonframshift_indel = 3022
  filter_synonymous_unknown = 3023
  filter_lcr_indel = 3024
  #filter_named_snp = 3030 # This variant is recorded in dbSNP and therefore MAY be too common for our interest

  def __init__(self, errCode, message, currentMethod):
    # Call the base class constructor with the parameters it needs
    #       super(BurdenException, self).__init__(message) # causes a deprecation warning
    self.errCode = errCode
    self.message = message
    self.stack = [currentMethod]
    return
  def add_to_stack(self, newMethod):
    self.stack.append(newMethod)
    return
  def add_to_message(self, newMessage):
    self.message = '{0}\n{1}'.format(self.message, newMessage)
    return
  def get_message(self):
    printableStack = ':'.join(reversed(self.stack)) # 'reversed()' returns an iterable, not a list
    msgString = '{0}\nSTACK TRACE->{1}'.format(self.message, printableStack)
    return msgString

#########################
##### MutationCount #####
class MutationCount:
  # This class manages all the mutations (counted by variants and samples) for each gene
  # The reason for using this class is because the gene records do not occur in uninterrupted order in the VCF file
  # Therefore, context switching must be implemented to allow proper counting of compound heterozygous variants
  # Each gene name is a key into the 'countingSets' hash
  # The value matching each key is a list which contains three more lists for missense, deleterious & lof variant types
  # Each of these inner lists contains two more lists which contain names of samples that have heterozygous or homozygous genotypes
  # When the number of genes stored in 'countingSets' passes a threshold, the oldest gene in the hash is processed and written to file
  # This happens without the calling class (ExomeBurden) knowing that a gene is being processed and deleted from this class
  missense = 0
  deleterious = 1
  lof = 2 # Loss of Function
  het = 0
  hom = 1
  geneCountThreshold = 10 # assuming that not more than 10 gene records are interleaved in the VCF file

  def __init__(self, vcfFilename, outputDir, excel):
    self.excel = excel
    self.open_table_files(vcfFilename, outputDir)
    self.addedGeneNames = [] # append new names and pop(0) processed names; FIFO order
    self.countingSets = {}
    return

  def __del__(self):
    if self.tableFH: self.tableFH.close()
    return

    ### open_table_files ###
    # called by MutationCount.init()
  def open_table_files(self, filename, outputDir):
    self.tableFH = open(os.path.join(outputDir, filename), 'w')
    if self.excel:
      self.write_excel_header()
    else:
      self.write_basic_header()
    return

    ### write_excel_header ###
    # This method prints a header that formats better when the resultant file is opened in a spreadsheet program
  def write_excel_header(self):
    self.tableFH.write('\t\t\t\t\t\tCounts By Variant\t\t\t\t\t\t\tx\t\t\t\t\t\t\tCounts By Sample\t\t\t\t\t\n')
    self.tableFH.write('Gene\tTotal_Missense\tTotal_Deleterious\tTotal_LoF\t')
    self.tableFH.write('Het_Missense\tHet_Deleterious\tHet_LoF\t')
    self.tableFH.write('HomComHet_Missense\tHomComHet_Deleterious\tHomComHet_LoF\t')
    self.tableFH.write('Hom_Missense\tHom_Deleterious\tHom_LoF\t')
    self.tableFH.write('X\t')
    self.tableFH.write('Total_Missense\tTotal_Deleterious\tTotal_LoF\t')
    self.tableFH.write('Het_Missense\tHet_Deleterious\tHet_LoF\t')
    self.tableFH.write('HomComHet_Missense\tHomComHet_Deleterious\tHomComHet_LoF\t')
    self.tableFH.write('Hom_Missense\tHom_Deleterious\tHom_LoF\n')
    return

    ### write_basic_header ###
    # This method prints a very basic, one-line header that is usually ignored when the one-tailed test script is called on the resultant file
  def write_basic_header(self):
    self.tableFH.write('Gene\tTotal_Missense\tTotal_Deleterious\tTotal_LoF\t')
    self.tableFH.write('Het_Missense\tHet_Deleterious\tHet_LoF\t')
    self.tableFH.write('HomComHet_Missense\tHomComHet_Deleterious\tHomComHet_LoF\t')
    self.tableFH.write('Hom_Missense\tHom_Deleterious\tHom_LoF\t')
    self.tableFH.write('Total_Missense\tTotal_Deleterious\tTotal_LoF\t')
    self.tableFH.write('Het_Missense\tHet_Deleterious\tHet_LoF\t')
    self.tableFH.write('HomComHet_Missense\tHomComHet_Deleterious\tHomComHet_LoF\t')
    self.tableFH.write('Hom_Missense\tHom_Deleterious\tHom_LoF\n')
    return

    ### add_mutation ###
    # called by external class
  def add_mutation(self, geneName, variantType, genotype, sampleName):
    if geneName not in self.countingSets.keys():
      self.countingSets[geneName] = [([],[]), ([],[]), ([],[])] # missense:(het, hom); deleterious:(het,hom); lof:(het,hom)
      self.addedGeneNames.append(geneName)
    self.countingSets[geneName][variantType][genotype].append(sampleName)
    if len(self.countingSets.keys()) > self.geneCountThreshold:
      self.process_oldest_gene()
    return

    ### process_oldest_gene ###
    # called by MutationCount.add_mutation() and MutationCount.flush_genes()
    # This consolidates all the mutations for a specific gene and writes them to their respective table files.
    # IDEALLY, this method would be called as a separate thread, but I do not have time to design & test that implementation.
    # The main logic here is to count the variants and the samples and figure out compound heterozygotes from single heterozygotes
    # for a given sample name. If a sample contributes only one heterozygote, then we count that single variant as 'het'.
    # However, if a sample contributes more than one heterozygote to a gene, it is considered a compound heterozygote.
  def process_oldest_gene(self):
    oldestGene = self.addedGeneNames.pop(0) # pop the first name off the list
    #if oldestGene == 'LFNG': pdb.set_trace() # FDFT1, PABPC3, MPRIP ABCA12
    variantCount = ([0,0,0], [0,0,0], [0,0,0]) # varaintCount = (het:[missense, deleterious, lof], cmpdHet[..], hom[..])
    sampleCount = ([0,0,0], [0,0,0], [0,0,0]) # same structure as 'variantCount'
    varTypeIndex = 0
    # walk through the list of variant types holding sample names for heterozygous and homozygous mutations in 'self.countingSets'
    while varTypeIndex < len(self.countingSets[oldestGene]):
      hetNamesCounted = []
      cmpdNamesCounted = []
      homNamesCounted = []
      hetSamples = self.countingSets[oldestGene][varTypeIndex][self.het] # all sample names
      uniqSampleNames = list(set(hetSamples)) # unique sample names
      for name in uniqSampleNames:
        # count the number of times this name occurs in 'hetSamples'
        occurrences = hetSamples.count(name)
        if 1 < occurrences:
          # more than one het mutation means this is a "compound heterozygote"
          variantCount[1][varTypeIndex] += occurrences * 2
          if name not in cmpdNamesCounted:
            sampleCount[1][varTypeIndex] += 1
            cmpdNamesCounted.append(name)
        # increment heterozygote count for variant and sample
        variantCount[0][varTypeIndex] += occurrences
        if name not in hetNamesCounted:
          sampleCount[0][varTypeIndex] += 1
          hetNamesCounted.append(name)
      
      homSamples = self.countingSets[oldestGene][varTypeIndex][self.hom]
      uniqSampleNames = list(set(homSamples))
      variantCount[2][varTypeIndex] = len(homSamples) * 2
      for name in uniqSampleNames:
        if name not in homNamesCounted:
          sampleCount[2][varTypeIndex] += 1
          homNamesCounted.append(name)
      varTypeIndex += 1
    if sum(sampleCount[0]) > 0 or sum(sampleCount[1]) > 0 or sum(sampleCount[2]) > 0: 
      self.write_to_table(oldestGene, variantCount, sampleCount) # write every gene
    self.countingSets[oldestGene] = None # remove data stored at key:geneName
    self.countingSets.pop(oldestGene, None) # remove key:geneName
    return

  ### write_to_table ###
  # called by process_oldest_gene()
  # This method is insanely ugly, I did not have the time to write it in a modular fashion
  # varCount = het(missense, deleterious, lof), cmpdHet(..), hom(..)
  # smpCount has the same layout as varCount
  # Mis -> missense; Del -> deleterious; LoF -> Loss-of-Function
  def write_to_table(self, geneName, varCount, smpCount):
    #if 'MPRIP' == geneName: pdb.set_trace() # PABPC3 # debug code for interlaced genes in the genome

    # compute and write variant allele counts
    hetLoF = varCount[0][2]
    hetDel = hetLoF + varCount[0][1]
    hetMis = hetDel + varCount[0][0]

    homPlusCmpdHetLoF = varCount[1][2] + varCount[2][2]
    homPlusCmpdHetDel = homPlusCmpdHetLoF + varCount[1][1] + varCount[2][1]
    homPlusCmpdHetMis = homPlusCmpdHetDel + varCount[1][0] + varCount[2][0]

    homLoF = varCount[2][2]
    homDel = homLoF + varCount[2][1]
    homMis = homDel + varCount[2][0]

    totLoF = hetLoF + homPlusCmpdHetLoF
    totDel = hetDel + homPlusCmpdHetDel
    totMis = hetMis + homPlusCmpdHetMis

    if self.excel:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\tx\t'.\
                         format(geneName,
                                totMis, totDel, totLoF,
                                hetMis, hetDel, hetLoF,
                                homPlusCmpdHetMis, homPlusCmpdHetDel, homPlusCmpdHetLoF,
                                homMis, homDel, homLoF))
    else:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t'.\
                         format(geneName,
                                totMis, totDel, totLoF,
                                hetMis, hetDel, hetLoF,
                                homPlusCmpdHetMis, homPlusCmpdHetDel, homPlusCmpdHetLoF,
                                homMis, homDel, homLoF))

    # compute and write sample variant counts
    hetLoF = smpCount[0][2]
    hetDel = hetLoF + smpCount[0][1]
    hetMis = hetDel + smpCount[0][0]

    homPlusCmpdHetLoF = smpCount[1][2] + smpCount[2][2]
    homPlusCmpdHetDel = homPlusCmpdHetLoF + smpCount[1][1] + smpCount[2][1]
    homPlusCmpdHetMis = homPlusCmpdHetDel + smpCount[1][0] + smpCount[2][0]

    homLoF = smpCount[2][2]
    homDel = homLoF + smpCount[2][1]
    homMis = homDel + smpCount[2][0]

    totLoF = hetLoF + homPlusCmpdHetLoF
    totDel = hetDel + homPlusCmpdHetDel
    totMis = hetMis + homPlusCmpdHetMis

    if self.excel:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.\
                         format(totMis, totDel, totLoF,
                                hetMis, hetDel, hetLoF,
                                homPlusCmpdHetMis, homPlusCmpdHetDel, homPlusCmpdHetLoF,
                                homMis, homDel, homLoF))
    else:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.\
                         format(totMis, totDel, totLoF,
                                hetMis, hetDel, hetLoF,
                                homPlusCmpdHetMis, homPlusCmpdHetDel, homPlusCmpdHetLoF,
                                homMis, homDel, homLoF))
      
    self.tableFH.flush()
    return

  ### flush_genes ###
  # called by external class
  def flush_genes(self):
    while 0 < len(self.addedGeneNames):
      self.process_oldest_gene()
    return

#######################
##### ExomeBurden #####
# Each line in the VCF file is processed independently. The only exception to this is counting variations that are
# heterozygous to track compound heterozygous variants.
#
# Most of the information we care about resides in the INFO field of the annotation. This data is packed down a few levels,
# and is extracted with the "<string>.split(<delimiter>)" method. This yields 'key=value' pairs, which are also split.
# The 'key=value' pairs are stored in the hash "infoHash". The list "screenedInfoField" is a list of 'key' values that are
# used for filtering the variant. Thus, we store all values in the INFO field, but only use a few to determine if the variant is
# worth further consideration.
#
# Because this is bioinformatics, there are exceptions to how the data is stored or labeled in the VCF file that require
# extra code to deal with. The most egregious is that there are two different names for ExAC data, depending on the annotation
# applied to the the VCF file. In other words, some VCF files use one string for the ExAC key while other VCFs use a different string.
# To allow for this, "screenedInfoFields" is allowed to have a list instead of a string in one of it's elements. If a list is detected,
# the elements in the list are matched against values in the VCF file. One a match is found, that value is used to overwrite
# the list. The code in method "screened_info_field()" performs this.
#
## OBSOLETE WITH THE USE OF gnomad_ExAC DATA
## vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# Another example of exceptional data handling is the version strings appended to the other two SNV population databases.
# To deal with this, "screenedInfoField" has a string for each of these databases without any version information. These values
# are matched up against the value in the VCF file. When a match is found THE VALUE IN "screenedInfoFields" IS USED to allow
# for smooth processing of the data values for these databases.
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
## OBSOLETE WITH THE USE OF gnomad_ExAC DATA
#
# Another wrinkle to the filtering of variants is that those variants that pass all filters are to have their records "melted".
# This means that, for each sample, all data is written out in tab-delimited format. This creates a very great redundancy in the
# resulting data table. This is why ALL VALUES in the INFO field are saved in "infoHash". This allows all these fields to be
# written to the melted file without their values being changed. Despite the apparently wastefulness of this file, it is very
# useful for post-processing beyond just achieving p-values from the Fisher exact test.
#
# Besides the melted VCF file, the end result of this program is to generate a file of variant counts. This file should be generated
# for both your cases as well as an ethnically comparable background. The resultant variant counts tables are compared with a
# Fisher exact test to determine if the number of significant variants contributed to each gene is statistically significant.

class ExomeBurden:
  allowedFuncDelFilters = ('CADD', 'metaSVM', 'radialSVM', 'MPC') # functional deleteriousness filter
  #allowedPopFilterNames = ('esp6500', '1000g', 'exac') # population frequency filter, OBSOLETE WITH gnomad_ExAC DATA
  desiredFormatFields = ['GT', 'AD', 'DP', 'GQ', 'PL'] # melted table output formatting of FORMAT field
  lcrPickleFile = '/home/user/data/lcr_positions_nov2015_v1.pickle'  # low complexity regions
  lcrHash = None

  # These hold the command line paramters specified by the user at runtime
  desiredPopFreq = None
  vcfFH = None # variant call format - file handle
  funcDel = None
  plDiffCutoff = None
  excel = False
  outputDir = None
  removeNonFrameShiftInDel = False
  doNotFilter = []
  requiredDepth = 8 # default required read depth value
  africaSpecial = False # this was a one-off special case, it can safely be ignored

  # Variables that hold parsed values from the file
  infoHash = None
  alleleList = None

  # Values used to melt VCF records
  meltFH = None
  prefixFieldNames = []
  infoFieldNames = []
  formatFieldNames = []
  
  # Please note that all arrays that count variants are in the same order as this one
  allowedVariantTypes = ('missense', 'deleterious', 'lof')
  # hardcoded indexes into 'mutationCount'
  het = 0
  hom = 1
  # genotypes from sample data that we do not want to process (they don't contain any trackable variants)
  unwantedGenotypes = ('0/0', '0|0', './.', '.|.')

  ### init ###
  def __init__(self):
    self.read_lcr_pickle_file()
    vcfFileName = self.parse_cmdline_arguments()
    
    # create output file names
    filename = os.path.basename(vcfFileName)
    filename = os.path.splitext(filename)[0]
    countTableName = '{0}_{1}_{2}_plDiff{3}_counts_table.tsv'.format(filename, self.desiredPopFreq, self.funcDel, self.plDiffCutoff)
    self.mutCount = MutationCount(countTableName, self.outputDir, self.excel)
    meltedDataFile = '{0}_{1}_{2}_plDiff{3}_variant.table'.format(filename, self.desiredPopFreq, self.funcDel, self.plDiffCutoff)
    meltedDataFile = os.path.join(self.outputDir, meltedDataFile)
    self.meltFH = open(meltedDataFile, 'w')
    return

  def __del__(self):
    if self.vcfFH: self.vcfFH.close()
    if self.meltFH: self.meltFH.close()
    return

  ### read_lcr_pickle_file ###
  def read_lcr_pickle_file(self):
    if not os.path.isfile(self.lcrPickleFile):
      sys.stderr.write('ERROR - unable to access low complexity region pickle file: {0}\n'.format(self.lcrPickleFile))
      sys.stderr.write('\tYou need to edit this program and change the value of \"lcrPickleFile\"\n')
      sys.exit(1)
    fh = open(self.lcrPickleFile, 'r')
    self.lcrHash = cPickle.load(fh)
    fh.close()
    return
    
  ### parse_cmdline_arguments ###
  # see "print_short_help()" for an explanation of the command line options
  def parse_cmdline_arguments(self):
    if 1 == len(sys.argv):
      self.print_short_help()

    argIndex = 1
    while argIndex < len(sys.argv):
      if '-f' == sys.argv[argIndex] or '--frequency' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"frequency\" argument\n')
          self.print_short_help()
        self.parse_frequency(sys.argv[argIndex + 1])
        argIndex += 2
      elif '-v' == sys.argv[argIndex] or '--vcf' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"vcf\" argument\n')
          self.print_short_help()
        self.open_vcf_file(sys.argv[argIndex + 1])
        vcfFileName = sys.argv[argIndex + 1]
        argIndex += 2
      elif '-d' == sys.argv[argIndex] or '--funcDel' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"funcDel\" argument\n')
          self.print_short_help()
        self.parse_functional_deleteriousness_filter(sys.argv[argIndex + 1])
        argIndex += 2
      elif '-p' == sys.argv[argIndex] or '--pldiff-cutoff' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"plDiff-cutoff\" argument\n')
          self.print_short_help()
        try:
          self.plDiffCutoff = int(sys.argv[argIndex+1])
        except ValueError as exc:
          sys.stderr.write('ERROR - \"pldiff-cutoff\" must be an integer\n')
          self.print_short_help()
        argIndex += 2
      elif '-x' == sys.argv[argIndex] or '--excel' == sys.argv[argIndex]:
        self.excel = True
        argIndex += 1
      elif '-o' == sys.argv[argIndex] or '--outputDir' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"outputDir\" argument\n')
          self.print_short_help()
        self.parse_output_directory(sys.argv[argIndex + 1])
        argIndex += 2
      elif '-s' == sys.argv[argIndex] or '--shift' == sys.argv[argIndex]:
        self.removeNonFrameShiftInDel = True
        argIndex += 1
      elif '-c' == sys.argv[argIndex] or '--coverage' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"coverage\" argument\n')
          self.print_short_help()
        self.requiredDepth = int(sys.argv[argIndex+1])
        argIndex += 2
      elif '-a' == sys.argv[argIndex] or '--africa-special' == sys.argv[argIndex]:
        self.africaSpecial = True
        argIndex += 2
      else:
        sys.stderr.write('ERROR - illegal format for command line parameters\n')
        self.print_short_help()

    if None == self.desiredPopFreq or None == self.vcfFH or None == self.funcDel or None == self.plDiffCutoff:
      sys.stderr.write('ERROR - \"frequency\", \"VCF file\", \"functional deleteriousness\" and \"plDiff cutoff\" must be specified\n')
      self.print_short_help()

    if None == self.outputDir:
      self.outputDir = '.'

    return vcfFileName

  ### parse_frequency ###
  def parse_frequency(self, freqArg):
    try:
      freq = float(freqArg)
    except ValueError as exc:
      sys.stderr.write('ERROR - \"frequency\" value must contain a decimal point (floats only): {0}\n')
      self.print_short_help()
    if freq > 1 or freq <= 0:
      sys.stderr.write('ERROR - incorrect range for frequency; must be between [0 and 1): {0}\n')
      self.print_short_help()
    self.desiredPopFreq = freq
    return

  ### open_vcf_file ###
  def open_vcf_file(self, filename):
    try:
      self.vcfFH = open(filename, 'r')
    except IOError as exc:
      sys.stderr.write('ERROR - unable to open VCF file: {0}\n'.format(filename))
      self.print_short_help()
    return

  ### parse_functional_deleteriousness ###
  def parse_functional_deleteriousness_filter(self, filterName):
    index = 0
    while index < len(self.allowedFuncDelFilters):
      if self.allowedFuncDelFilters[index].lower() == filterName.lower():
        self.funcDel = self.allowedFuncDelFilters[index]
        return
      index += 1
    sys.stderr.write('ERROR - illegal value for functional deleteriousness: {0}\n'.format(filterName))
    self.print_short_help()
    return

  ### parse_output_directory ###
  def parse_output_directory(self, dirName):
    if not os.path.exists(dirName):
      sys.stdout.write('ERROR - output directory does not exist: {0}\n'.format(dirName))
      self.print_short_help()
    if not os.path.isdir(dirName):
      sys.stderr.write('ERROR - output directory is not a directory: {0}\n'.format(dirName))
      self.print_short_help()
    self.outputDir = dirName
    return

  ### print_short_help ###
  def print_short_help(self):
    sys.stderr.write('usage: {0} parameters [options] \n\n'.format(os.path.basename(sys.argv[0])))
    sys.stderr.write('Parameters:\n')
    sys.stderr.write('\t-f | --frequency\tmax desired frequency of variants in population databases\n')
    sys.stderr.write('\t-v | --vcf\t\tVCF file containing variants to be analyzed\n')
    sys.stderr.write('\t-d | --funcDel\t\tmethod to use for filtering functional deleteriousness: {0}\n'.format(', '.join(self.allowedFuncDelFilters)))
    sys.stderr.write('\t-p | --pldiff-cutoff\tminimum required value of plDiff for a sample to be included in further processing\n')
    sys.stderr.write('\t\t 7:low stringency 8:high stringency\n')
    sys.stderr.write('Options\n')
    sys.stderr.write('\t-o | --outputDir\toutput directory to write files; default: current directory\n')
    sys.stderr.write('\t-x | --excel\t\tFLAG: use to cause more header information for display in spreadsheet\n')
    sys.stderr.write('\t-s | --shift\t\tFLAG: use to specify that nonframeshift indels should be filtered out\n')
    sys.stderr.write('\t-c | --coverage\t\trequired read depth for a variant to be accepted; default: 8\n')
    sys.stderr.write('\t-a | --africa-special\tFLAG: filter specific info fields\n')
    sys.stderr.flush()
    sys.exit(1)
    return
    
  ### read_vcf_file ###
  # called by instantiated object (see __main__ below)
  def read_vcf_file(self):
    headerLine = self.parse_vcf_meta_data()

    # The 'gnomAD_*_all' and various 'exac' entries are used as an index into 'self.infoHash'. Please do not change their position in this list.
    self.screenedInfoFields = ['ExonicFunc.refGene', 'Gene.refGene', 'Func.refGene', ['ExAC_ALL','exac03', 'ExAC_nontcga_ALL']]
    if 'gnomAD_exome_ALL' in self.infoFieldNames:
      self.screenedInfoFields.append('gnomAD_exome_ALL')
      self.screenedInfoFields.append('gnomAD_genome_ALL')
    if self.africaSpecial:
      self.screenedInfoFields.append('ExAC_nontcga_AFR')
      self.screenedInfoFields.append('gnomAD_exome_AFR')
      self.screenedInfoFields.append('gnomAD_exome_AFR')
    
    if 'CADD' == self.funcDel:
      self.screenedInfoFields.append('CADD_phred')
    elif 'metaSVM' == self.funcDel:
      self.screenedInfoFields.append('MetaSVM_pred')
    elif 'radialSVM' == self.funcDel:
      self.screenedInfoFields.append('RadialSVM_pred')
    elif 'MPC' == self.funcDel:
      self.screenedInfoFields.append('mpc')
    
    if 0 == len(self.infoFieldNames) or 0 == len(self.formatFieldNames):
      errMsg = 'ERROR - unable to find INFO or FORMAT meta-data...does the VCF file have a meta-data header?'
      raise BurdenException(BurdenException.parsing_error, errMsg, 'read_vcf_file')
    if not headerLine:
      errMsg = 'ERROR - unable to find header line in VCF file\n'
      raise BurdenException(BurdenException.parsing_error, errMsg, 'read_vcf_file')
    self.sampleNames = headerLine.split('\t')[9:]
    self.write_melted_header_line()
    try:
      self.read_variants()
    except BurdenException as exc:
      self.vcfFH.close()
      exc.add_to_stack('read_vcf_file')
      raise exc
    self.mutCount.flush_genes()
    self.vcfFH.close()
    return

  ### parse_vcf_meta_data ###
  # Populate "infoFieldNames" and "formatFieldNames".
  # Return the header line, beginning with '#CHROM'.
  def parse_vcf_meta_data(self):
    annovarDate = False
    line = self.vcfFH.readline().rstrip()
    while line:
      if '#CHROM' == line[0:6]:
        return line # return header line
      if '##INFO' == line[0:6]:
        identifier = re.search('##INFO=<ID=(.*?),Number', line)
        if not identifier:
          errMsg = 'ERROR - Unable to parse info name from meta data line\n{0}\n'.format(line)
          raise BurdenException(BurdenException.parsing_error, errMsg, 'parse_vcf_meta_data')
        if 'ALLELE_END' == identifier.group(1) or\
           'NEGATIVE_TRAIN_SITE' == identifier.group(1) or\
           'POSITIVE_TRAIN_SITE' == identifier.group(1) or\
           'CSQ' == identifier.group(1):
          line = self.vcfFH.readline().rstrip()
          continue
        if 'ANNOVAR_DATE' == identifier.group(1):
          annovarDate = True
        if annovarDate:
          self.infoFieldNames.append(identifier.group(1))
        else:
          self.prefixFieldNames.append(identifier.group(1))
      elif '##FORMAT' == line[0:8]:
        identifier = re.search('##FORMAT=<ID=(.*?),Number', line)
        if not identifier:
          errMsg = 'ERROR - Unable to parse format name from meta data line\n{0}\n'.format(line)
          raise BurdenException(BurdenException.parsing_error, errMsg, 'parse_vcf_meta_data')
        self.formatFieldNames.append(identifier.group(1))
      line = self.vcfFH.readline().rstrip()
    return

  ### write_melted_header_line ###
  def write_melted_header_line(self):
    self.meltFH.write( 'SAMPLE\tEFFECTIVE VARIANT TYPE\t{0}\tCHROM\tPOSITION\tREFERENCE ALLELE\tALTERNATE ALLELES\tQUALITY\t'
                       .format('\t'.join(self.desiredFormatFields)) )
    self.meltFH.write( 'FILTER\tPLDIFF/DEPTH\t{0}\tALLELE\tVARIANT TYPE\t{1}\tALLELE\tVARIANT TYPE\t{1}\n'
                       .format('\t'.join(self.prefixFieldNames), '\t'.join(self.infoFieldNames)) )
    return
    
  ### read_variants ###
  # called by read_vcf_file()
  # This is the high-level method that performs most of the work on each record and variant.
  # Due to the size of most VCF files, we need to read line-by-line, as opposed to storing the file contents in memory.
  # For each record:
  # First, determine the position of the GQ, DP and PL values in the FORMAT field.
  # Second, isolate data for each alternate allele.
  # Third, the data for each alternate allele is passed through filters and marked "False" if it fails any of them.
  # Fourth, the number of samples in this record that have passed the filters are counted.
  # Fifth, the entire record is written to a variant table; each sample is melted out to create a large table of all variant data.
  #
  # A word about "alleleList"
  # Each record has a reference allele and at least one alternate allele. "alleleList" holds information for all of these alleles.
  # Each list element consists of another list of three elements:
  # 0) The actual allelele ie: "A", "GT"
  # 1) A hash containing the "Key=Value" pairs from the INFO field for this allele.
  # 2) This value represents the type of mutation this variant would introduce. It's initialized to 'pass' for debugging.
  #    If this variant fails any filter, this value is set to 'fail'. Otherwise, it's set to the integer value of one of
  #    the variables in MutationCount: lof:2, deleterious:1, missense:0
  # 3) Boolean indicating if this allele is an "indel" AND the location is in a Low Complexity Region (lcr).
  def read_variants(self):
    line = self.vcfFH.readline().rstrip()
    while line:
      formatFields = line.split('\t')[8].split(':')
      if 'GQ' not in formatFields:
        errMsg = 'ERROR - genotype quality (GQ) not present in format field: {0}\n'.format(line.split('\t')[8])
        raise BurdenException(BurdenException.parsing_error, errMsg, 'read_variants')
      if 'DP' not in formatFields:
        errMsg = 'ERROR - read depth (DP) not present in format field: {0}\n'.format(line.split('\t')[8])
        raise BurdenException(BurdenException.parsing_error, errMsg, 'read_variants')
      self.GQindex = formatFields.index('GQ')
      self.DPindex = formatFields.index('DP')
      if 'PL' in formatFields:
        self.PLindex = formatFields.index('PL')
      else:
        self.PLindex = '.'
      plDiffScores = None # tracks specific filter failures OR the quality score for heterozygous samples, if filters are passed
      columns = line.split('\t')
      #CHROM  POS       ID  REF        ALT         QUAL     FILTER       INFO       FORMAT
      (chrom, position, ID, refAllele, altAlleles, quality, filterField, infoField, formatField) = columns[0:9]
      samples = columns[9:]
      self.alleleList = [[refAllele, None, 'fail', False]] # place reference allele at head of list
      prefixHash = self.populate_allele_list(altAlleles.split(','), infoField)

      ### DEBUG
      #if chrom == '6' and position == '32485520': pdb.set_trace()
      ### END DEBUG

      index = 1 # skip reference allele
      while index < len(self.alleleList):
        # self.alleleList[index][1] points to the hash of INFO field values that are currently being filtered
        if '.' == self.alleleList[index][1]['Gene.refGene']: # skip records that do not have a gene name
          self.alleleList[index][2] = 'fail'
          index += 1
          continue
        try:
          self.filter_on_database_population_frequencies(index)
          self.filter_on_refgene_values(index, chrom, position)
        except BurdenException as exc:
          if exc.errCode > BurdenException.filtered_out:
            index += 1
            continue
          else:
            if BurdenException.parsing_error == exc.errCode:
              exc.add_to_stack('read_variants')
              exc.add_to_message('in record with Chrom:{0} and Position:{1}'.format(chrom, position))
            raise exc
        self.classify_variant(index)
        index += 1 # end "while index > len(self.alleleList):" loop

      # If all alleles fail the filters, do not bother looking at each individual sample
      allFail = True
      for allele in self.alleleList:
        if 'fail' != allele[3]:
          allFail = False
      if allFail:
        plDiffScores = ['fail']*len(samples) # mark all samples as "fail"
      else:
        plDiffScores = self.count_genotypes(self.alleleList[1][1]['Gene.refGene'], samples)
      self.melt_vcf_record(line, plDiffScores, prefixHash)
      line = self.vcfFH.readline().rstrip() # end "while line:" loop
    return

  ### populate_allele_list ###
  # called by read_variants()
  # Parse out "KEY=VALUE" pairs for each alternate allele.
  # These pairs are stored in a hash named 'infoHash'.
  # This hash is then stored in 'alleleList'..
  # Each alternate allele gets it's own infoHash.
  def populate_allele_list(self, altAlleles, infoField):
    try:
      alleleData = self.parse_allele_info(infoField)
    except BurdenException as exc:
      exc.add_to_stack('populate_allele_list')
      raise exc
    annovarDates = alleleData[1::2] # [start:stop:stride]
    del alleleData[1::2]
    prefixHash = self.parse_prefix_data(alleleData[0])
      
    index = 1
    while index < len(alleleData):
      infoHash = self.parse_info_field(alleleData[index], annovarDates.pop(0))
      self.alleleList.append([altAlleles[index-1], infoHash, 'pass', False])
      index += 1
    return prefixHash

  ### parse_allele_info ###
  # The first line, 're.split', returns a list of all data surrounded by the regex, including
  # the data captured by the parentheses. ie:(.*?)
  # We then need to remove the first and last elements of that list, as they do not contain
  # data that we want.
  # We finally need to remove any empty elements with the call to filter.
  # Final format of "match":
  # a) prefix data
  # b) annovar date value
  # c) allele specific data
  # (b) and (c) repeat for each alternate allele.
  def parse_allele_info(self, infoField):
    match = re.split('ANNOVAR_DATE=(\d{4}-\d\d-\d\d);(.*?)ALLELE_END;', infoField)
    if 1 == len(match):
      match = re.split('ANNOVAR_DATE=(\d{4}-\d\d-\d\d);(.*?)ALLELE_END', infoField) # search w/o semicolon at end
      if 1 == len(match):
        errString = 'ERROR - unable to parse info field with use of \"re.split()\" function\n'
        raise BurdenException(BurdenException.parsing_error, errString, 'parse_allele_info')
    match = match[0:-1] # remove last element
    match = filter(None, match) # remove all empty elements
    return match

  ### parse_prefix_data ###
  def parse_prefix_data(self, prefixData):
    prefixHash = {}
    keyValuePairs = prefixData.split(';')
    index = 0
    while index < len(keyValuePairs):
      if '=' in keyValuePairs[index]:
        (key, value) = keyValuePairs[index].split('=')
        prefixHash[key] = value
      index += 1
    return prefixHash

  ### parse_info_field ###
  def parse_info_field(self, infoField, annovarDate):
    infoHash = {}
    keyValuePairs = infoField.split(';')
    index = 0
    infoHash['ANNOVAR_DATE'] = annovarDate
    while index < len(keyValuePairs):
      if '=' in keyValuePairs[index]:
        (name, value) = keyValuePairs[index].split('=')
        key = self.screened_info_field(name)
        if key:
          infoHash[key] = value
        else:
          infoHash[name] = value
      index += 1
    try:
      self.verify_info_fields(infoHash)
    except BurdenException as exc:
      exc.add_to_stack('populate_allele_list')
      raise exc
    return infoHash
    
  ### screened_info_field ###
  # called by parse_info_field()
  # Allow more general matching of INFO field names to allow matching of different versions of databases.
  # If one of the screened field names is contained in 'currentField', return the value that matched against 'currentField'.
  # This way, we dictate what keys are used in 'self.infoHash', when matched later in 'filter_on_database_population_frequencies()'.
  # This also frees us from using any appended strings to indicate database version.
  #
  # There are also elements in "screenedInfoFields" that are actually a list of values. This represents the fact that different
  # VCF files have different strings identifying the same data. Once the first record is read, the value used by the current
  # VCF file is used to replace the element containing the list.
  def screened_info_field(self, currentField):
    index = 0
    while index < len(self.screenedInfoFields):
      # only process strings in this way; entries that are lists are dealt with below
      if isinstance(self.screenedInfoFields[index], basestring):
        fieldName = self.screenedInfoFields[index]
        matchObj = re.search(fieldName, currentField) # does "fieldName" occur in "currentField"?
        if matchObj:
          return fieldName
      else:
        # This chunk of code only needs to run once for each list in "screenedInfoFields".
        # After that, the field in "screenedInfoFields" has been changed to a string.
        degenerateFieldNames = self.screenedInfoFields[index]
        for fieldName in degenerateFieldNames:
          matchObj = re.search(fieldName, currentField)
          if matchObj:
            self.screenedInfoFields[index] = fieldName # overwrite the list with the string matching the current VCF data
            return fieldName
      index += 1
    return False

  ### verify_info_fields ###
  # called by parse_info_field()
  # Check that all the fields we want to filter on are present
  # If you want to filter on other fields, you need to edit screenedInfoFields at the top of this class definition
  def verify_info_fields(self, infoHash):
    missingFields = []
    index = 0
    while index < len(self.screenedInfoFields):
      fieldName = self.screenedInfoFields[index]
      if fieldName not in infoHash.keys():
        if not isinstance(fieldName, basestring):
          errString = 'ERROR - Unable to find one of the matching field names in : {0}\n'.format(','.join(fieldName))
          raise BurdenException(BurdenException.parsing_error, errString, 'verify_info_fields')
        missingFields.append(fieldName)
      index += 1
    if missingFields:
      errString = 'ERROR - The following fields were not present in the INFO field:\n{0}\n'.format(','.join(missingFields))
      raise BurdenException(BurdenException.parsing_error, errString, 'verify_info_fields')
    return

  ### filter_on_database_population_frequencies ###
  # called by read_variants()
  def filter_on_database_population_frequencies(self, index):
    infoHash = self.alleleList[index][1]

    if 'exac' in self.screenedInfoFields[3].lower():
      if '.' == infoHash[self.screenedInfoFields[3]]:
        exacValue = 0
      else:
        exacValue = float(infoHash[self.screenedInfoFields[3]])
      if exacValue > self.desiredPopFreq:
        self.alleleList[index][2] = 'fail'
        raise BurdenException(BurdenException.filter_freq_exac_too_high, '-', 'filter_on_database_population_frequencies')
    else:
      errMsg = 'unable to parse ExAC population frequency'
      raise BurdenException(BurdenException.parsing_error, errMsg, 'filter_on_database_population_frequencies')

    if 'gnomAD_exome_ALL' in self.screenedInfoFields:
      if '.' == infoHash['gnomAD_exome_ALL']:
        gnomadExome = 0
      else:
        gnomadExome = float(infoHash['gnomAD_exome_ALL'])
      if gnomadExome > self.desiredPopFreq:
        self.alleleList[index][2] = 'fail'
        raise BurdenException(BurdenException.filter_freq_gnomad_exome_too_high, '-', 'filter_on_database_population_frequencies')
      
      if '.' == infoHash['gnomAD_genome_ALL']:
        gnomadGenome = 0
      else:
        gnomadGenome = float(infoHash['gnomAD_genome_ALL'])
      if gnomadGenome > self.desiredPopFreq:
        self.alleleList[index][2] = 'fail'
        raise BurdenException(BurdenException.filter_freq_gnomad_genome_too_high, '-', 'filter_on_database_population_frequencies')

    if self.africaSpecial:
      # ExAC_nontcga_AFR        gnomAD_exome_AFR        gnomAD_genome_AFR
      africaCutoff = 5e-04
      if '.' == infoHash['ExAC_nontcga_AFR']:
        exac = 0
      else:
        exac = float(infoHash['ExAC_nontcga_AFR'])
      if africaCutoff < exac:
        self.alleleList[index][2] = 'fail'
        raise BurdenException(BurdenException.filter_freq_exac_too_high, '-', 'filter_on_database_population_frequencies')

      if '.' == infoHash['gnomAD_exome_AFR']:
        exome = 0
      else:
        exome = float(infoHash['gnomAD_exome_AFR'])
      if '.' == infoHash['gnomAD_genome_AFR']:
        genome = 0
      else:
        genome = float(infoHash['gnomAD_genome_AFR'])
      if africaCutoff < exome or africaCutoff < genome :
        self.alleleList[index][2] = 'fail'
        raise BurdenException(BurdenException.filter_freq_gnomad_genome_too_high, '-', 'filter_on_database_population_frequencies')
    return

  ### classify_variant ###
  # called by read_variants()
  def classify_variant(self, index):
    infoHash = self.alleleList[index][1]
    dScore = self.determine_deleteriousness(infoHash)
    if ('splicing' in infoHash['Func.refGene'] or
        'frameshift' == infoHash['ExonicFunc.refGene'][0:10] or
        'stopgain' in infoHash['ExonicFunc.refGene']):
      variantType = self.mutCount.lof # 2 LoF
    elif (('nonsynonymous_SNV' == infoHash['ExonicFunc.refGene'] and 'deleterious' == dScore) or
          'nonframeshift' in infoHash['ExonicFunc.refGene'] or 'stoploss' in infoHash['ExonicFunc.refGene']):
      variantType = self.mutCount.deleterious # 1
    else:
      variantType = self.mutCount.missense # 0
    self.alleleList[index][2] = variantType
    return

  ### determine_deleteriousness ###
  # called by classify_variant()
  def determine_deleteriousness(self, infoHash):
    if 'metaSVM' == self.funcDel:
      if '.' == infoHash['MetaSVM_pred']:
        #infoHash['MetaSVM_pred'] = 'T'
        return 'tolerated' # assuming no value means score of "tolerated"
      if 'T' == infoHash['MetaSVM_pred']:
        return 'tolerated'
      else:
        return 'deleterious'
      
    if 'radialSVM' == self.funcDel:
      if '.' == infoHash['RadialSVM_pred']:
        #infoHash['RadialSVM_pred'] = 'T'
        return 'tolerated'
      if 'T' == infoHash['RadialSVM_pred']:
        return 'tolerated'
      else:
        return 'deleterious'
      
    if 'CADD' == self.funcDel:
      if '.' == infoHash['CADD_phred']:
        #infoHash['CADD_phred'] = 0
        return 'tolerated'
      caddPhred = float(infoHash['CADD_phred'])
      #if 15 > caddPhred:
      if 20 > caddPhred:
        return 'tolerated'
      else:
        return 'deleterious'

    if 'MPC' == self.funcDel:
      if '.' == infoHash['mpc']:
        return 'tolerated'
      else:
        mpcValue = float(infoHash['mpc'])
      if mpcValue >= 2:
        return 'deleterious'
      else:
        return 'tolerated'
    return None
                                              
  ### filter_on_refgene_values ###
  # called by read_variants()
  def filter_on_refgene_values(self, index, chrom, position):
    infoHash = self.alleleList[index][1]
    if 'exonic' != infoHash['Func.refGene'] and 'splicing' != infoHash['Func.refGene'] and \
       'exonic\x3bsplicing' != infoHash['Func.refGene']:
      self.alleleList[index][2] = 'fail'
      raise BurdenException(BurdenException.filter_nonexonic_nonsplicing, '-', 'filter_on_refgene_values')
    if self.removeNonFrameShiftInDel:
      if 'nonframeshift_insertion' == infoHash['ExonicFunc.refGene'] or \
         'nonframeshift_deletion'  == infoHash['ExonicFunc.refGene']:
        self.alleleList[index][2] = 'fail'
        raise BurdenException(BurdenException.filter_nonframshift_indel, '-', 'filter_on_refgene_values')
    if 'synonymous_SNV' == infoHash['ExonicFunc.refGene'] or 'unknown' == infoHash['ExonicFunc.refGene']:
      self.alleleList[index][2] = 'fail'
      raise BurdenException(BurdenException.filter_synonymous_unknown, '-', 'filter_on_refgene_values')
    if 'insertion' in infoHash['ExonicFunc.refGene'] or 'deletion' in infoHash['ExonicFunc.refGene']:
      if self.inside_low_complexity_region(chrom, position):
        self.alleleList[index][2] = 'fail' # not a valid allele
        self.alleleList[index][3] = True  # "indel" variant inside a Low Complexity Region
        raise BurdenException(BurdenException.filter_lcr_indel, '-', 'filter_on_refgene_values')
    return

  ### inside_low_complexity_region ###
  # This method determines if 'position' is within any regions specified in 'self.lcrHash'.
  # 'self.lcrHash' is keyed on chromosome number or letter in character format (ie: '1', '15', 'X')
  # The value pointed to by each chromosome is a sorted list of lists.
  # The secondary level of lists is the start and end of a low-complexity-region on that chromosome.
  # ie: '15' -> [ [100252703, 100252743], [100693010, 100693043], [101019679, 101019704], ... ]
  # If 'position' is out of range of all given locations, return False.
  # Otherwise, perform a binary search to find a region in which 'position' occurs.
  def inside_low_complexity_region(self, chrom, position):
    position = int(position)
    if chrom not in self.lcrHash.keys():
      errString = 'ERROR - variant chromosome not found in lcrHash keys\n'
      raise BurdenException(BurdenException.parsing_error, errString, 'inside_low_complexity_region')
    lcrLocations = self.lcrHash[chrom]
    if position < lcrLocations[0][0] or lcrLocations[-1][1] < position:
      return False
    result = self.destructive_binary_search(position, lcrLocations)
    return result

  ### destructive_binary_search ###
  # This method recurses on itself by shortening the list to be searched in a binary method
  # The 'start' positions are used to search for the [start, end] pair in which
  # 'start' is less than or equal to 'value' and the 'start' position
  # of the next [start, end] list is larger than 'value'.
  # At this point, the values of 'start' and 'end' are checked
  # to see if they contain 'value'.
  def destructive_binary_search(self, value, searchList):
    if 1 == len(searchList):
      return False
    if 2 == len(searchList):
      if searchList[0][0] <= value < searchList[0][1]:
        return True
      elif searchList[1][0] <= value < searchList[1][1]:
        return True
      else:
        return False
    index = (len(searchList) - 1) / 2
    if searchList[index][0] <= value < searchList[index+1][0]:
      if searchList[index][0] <= value <= searchList[index][1]:
        return True
      return False
    elif searchList[index][0] < value:
      return self.destructive_binary_search(value, searchList[index:])
    elif searchList[index][0] > value:
      return self.destructive_binary_search(value, searchList[0:index])
    return False

  ### count_genotypes ###
  # called by read_variants()
  # Samples are filtered in this method, as opposed to alternate alleles.
  # Essentiall, at least one of the alternate alleles has passed all filters
  # and we now need to know how many samples containing these alleles are valid.
  # For each sample:
  # Parse out the genotype info.
  # Check that it is not all reference alleles or unknown by seeing if it exists in self.unwantedGenotypes.
  # Check that the genotype quality (GQ) is greater than 20.
  # PLEASE NOTE: Genotypes of './.' or '.|.' are too short to index for GQ.
  #              These are already filtered out, but beware, in case the logic changes.
  # Isolate the two alleles by splitting on the character '/' or '|'.
  # Compare the two alleles (strand names are entirely arbitrary).
  # If they are the same, increment the homozygous count for this gene.
  # If they are not the same increment the heterozygous count for this sample.
  # This algorithm should deal with any number of alternate alleles.
  def count_genotypes(self, geneName, samples):
    #if 'PABPC3' == geneName: pdb.set_trace()
    plDiffScores = ['hom']*len(samples) # default to homozygous
    index = 0
    while index < len(samples):
      currentSample = samples[index]
      # filter out samples that have no information for this position
      currentGenotype = currentSample.split(':')[0]
      if currentGenotype in self.unwantedGenotypes:
        plDiffScores[index] = 'fail_genotype'
        index += 1
        continue
      # filter out samples with low genome quality for this position
      if 20 >= int(currentSample.split(':')[self.GQindex]):
        plDiffScores[index] = 'fail_GQ'
        index += 1
        continue
      # filter out low read depth samples
      currentDP = currentSample.split(':')[self.DPindex]
      if '.' == currentDP:
        currentDP = 0
      else:
        currentDP = int(currentDP)
      if currentDP < self.requiredDepth:
        plDiffScores[index] = 'fail_DP'
        index += 1
        continue

      sampleName = self.sampleNames[index]
      if '/' in currentGenotype:
        (plusStrand, minusStrand) = currentGenotype.split('/')
      elif '|' in currentGenotype:
        (plusStrand, minusStrand) = currentGenotype.split('|')
      else:
        errMsg = 'ERROR - Unable to parse genotype from sample data: {0}'.format(currentSample)
        raise BurdenException(BurdenException.parsing_error, errMsg, 'count_genotypes')
      sampleAlleles = [self.alleleList[int(plusStrand)], self.alleleList[int(minusStrand)]]
            
      if 'fail' == sampleAlleles[0][2] and 'fail' == sampleAlleles[1][2]:
        plDiffScores[index] = '.'
        index += 1
        continue
      if plusStrand == minusStrand:
        variantType = sampleAlleles[0][2]
        self.mutCount.add_mutation(geneName, variantType, self.mutCount.hom, sampleName)
      else:
        if '.' == self.PLindex:
          errMsg = 'ERROR - no PL value in format string to check plDiff on heterozygous sample'
          raise BurdenException(BurdenException.parsing_error, errMsg, 'count_genotypes')
        currentPL = currentSample.split(':')[self.PLindex]
        currentDP = currentSample.split(':')[self.DPindex]
        plDiffScores[index] = self.determine_plDiff_score(currentPL, currentDP)
        if plDiffScores[index] < self.plDiffCutoff:
          index += 1
          continue
        variantType = self.assign_variant(sampleAlleles)
        self.mutCount.add_mutation(geneName, variantType, self.mutCount.het, sampleName)
      index += 1
    return plDiffScores

  ### determine_plDiff_score ###
  # This method is a surrogate for visual inspection of variant pileups.
  # ASSUMPTIONS: The PL field contains at least two comma-separated values
  def determine_plDiff_score(self, plString, dpString):
    plNumbers = []
    for elt in plString.split(','):
      plNumbers.append(int(elt))
    plNumbers = sorted(plNumbers)
    plDiffScore = (plNumbers[1] - plNumbers[0]) / int(dpString)
    return plDiffScore

  ### assign_variant ###
  # Finds the variant type that is most deleterious.
  # We already know that both values are NOT "fail".
  def assign_variant(self, sampleAlleles):
    if 'fail' == sampleAlleles[0][2]:
      return sampleAlleles[1][2]
    elif 'fail' == sampleAlleles[1][2]:
      return sampleAlleles[0][2]
    else:
      return max(sampleAlleles[0][2], sampleAlleles[1][2])

  ### melt_vcf_record ###
  def melt_vcf_record(self, line, plDiffScores, prefixHash):
    columns = line.split('\t')
    #CHROM  POS       ID      REF       ALT       QUAL    FILTER       INFO       FORMAT
    (chrom, position, ID, refAllele, altAlleles, quality, filterField, infoField, formatField) = columns[0:9]
    samples = columns[9:]
    genomicLocationData = [chrom, position, refAllele, altAlleles, quality, filterField]
    prefixString = self.build_prefix_string(prefixHash)
    index = 0
    while index < len(self.sampleNames):
      self.write_melted_record(genomicLocationData, prefixString, formatField,
                               self.sampleNames[index], samples[index], plDiffScores[index])
      index += 1
    return

  ### build_prefix_string ###
  def build_prefix_string(self, prefixHash):
    prefixValues = []
    for name in self.prefixFieldNames:
      if not name in prefixHash.keys():
        prefixValues.append('.')
      else:
        prefixValues.append(prefixHash[name])
    return '\t'.join(prefixValues)

  ### write_melted_record ###
  # called by melt_vcf_record()
  # genomicLocationData: [chrom, position, refAllele, altAlleles, quality, filter]
  # formatField: Description of data in each sample for this line. ie: GT:AD:DP:GQ:PL or GT:AD:DP:GQ:PGT:PID:PL
  # sampleName: Name of the sample as specified in the header line starting with "#CHROM".
  # sample: Values as specified in "formatField". ie: 0/1:29,7:36:99:207,0,1389 or 0/1:23,4:27:94:0|1:17365_C_G:94,0,1581
  # plDiffScore: Quality score for heterozygous genotypes. May also hold reason for filter failure in "count_genotypes()".
  # PLEASE NOTE that both "formatField" and "sample" are delimited by colon (:).
  def write_melted_record(self, genomicLocationData, prefixString, formatField, sampleName, sample, plDiffScore):
    sampleData = sample.split(':')
    formatHeader = formatField.split(':')
    genotype = sampleData[formatHeader.index('GT')]
    if genotype in self.unwantedGenotypes:
      return
    if isinstance(plDiffScore, basestring) and 'fail' in plDiffScore:
      return # do not print data for samples that have failed tests in "count_genotypes()"

    ### DEBUG ###
    #if '1684347' == genomicLocationData[1] or '3753056' == genomicLocationData[1] or '6529183' == genomicLocationData[1]:
    #  pdb.set_trace()
    ### END DEBUG ###
    
    if '/' in genotype:
      (plusStrand, minusStrand) = genotype.split('/')
    else:
      (plusStrand, minusStrand) = genotype.split('|')
    sampleAlleles = [self.alleleList[int(plusStrand)], self.alleleList[int(minusStrand)]]
    # Check for 'fail' data value and whether this allele is an indel in a low complexity region
    if  ('fail' == sampleAlleles[0][2] and not sampleAlleles[1][3]) \
    and ('fail' == sampleAlleles[1][2] and not sampleAlleles[0][3]):
      return
    effectiveVariantType = self.assign_variant(sampleAlleles)
    effectiveVariantType = self.convert_variant_to_string(effectiveVariantType)
    try:
      var1String = self.convert_variant_to_string(sampleAlleles[0][2])
      var2String = self.convert_variant_to_string(sampleAlleles[1][2])
    except BurdenException as exc:
      exc.add_to_stack('write_melted_record')
    allele1Data = [var1String]
    allele2Data = [var2String]
    allele1Data.extend(self.melt_allele_data(sampleAlleles[0][1]))
    allele2Data.extend(self.melt_allele_data(sampleAlleles[1][1]))
    
    meltedSample = []
    formatIndex = 0
    while formatIndex < len(self.desiredFormatFields):
      try:
        meltedSample.append( sampleData[formatHeader.index(self.desiredFormatFields[formatIndex])] )
      except ValueError as exc:
        if 'PL' == self.desiredFormatFields[formatIndex]:
          meltedSample.append('.')
          pass
        else:
          errString = 'ERROR - format field {0} does not appear in the FORMAT field {0}\n'.format(formatField)
          raise BurdenException(BurdenException.parsing_error, errString, 'write_melted_record')
      formatIndex += 1

    lcrFlag1 = sampleAlleles[0][3]
    lcrFlag2 = sampleAlleles[1][3]
    #if lcrFlag1 and lcrFlag2:
    #  plDiffScore = 'lcr_region'
    if (lcrFlag1 and 'fail' == sampleAlleles[1][2]) \
    or (lcrFlag2 and 'fail' == sampleAlleles[0][2]):
      plDiffScore = 'lcr_region'

    #     0                    1                2              3                   4
    # sampleName, effective variant type,  meltedSample, genomicLocationData, plDiffScore,
    #     5          6         7           8         9
    # prefixData, allele1, allele1Data, allele2, allele2Data
    newRecord = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(sampleName,
                                                                            effectiveVariantType,
                                                                            '\t'.join(meltedSample),
                                                                            '\t'.join(genomicLocationData),
                                                                            plDiffScore,
                                                                            prefixString,
                                                                            sampleAlleles[0][0],
                                                                            '\t'.join(allele1Data),
                                                                            sampleAlleles[1][0],
                                                                            '\t'.join(allele2Data))
    self.meltFH.write(newRecord)
    self.meltFH.flush()
    return

  def convert_variant_to_string(self, variantType):
    if variantType == self.mutCount.missense:
      return 'missense'
    if variantType == self.mutCount.deleterious:
      return 'deleterious'
    if variantType == self.mutCount.lof:
      return 'lof'
    if 'fail' == variantType:
      return '-'
    errMsg = 'ERROR - variant type does not contain expected value\n'
    raise BurdenException(BurdenException.parsing_error, errMsg, 'convert_variant_to_string')

    
  ### melt_allele_data ###
  def melt_allele_data(self, infoHash):
    if None == infoHash: # reference allele
      meltedInfo = ['.']*len(self.infoFieldNames)
      return meltedInfo
    meltedInfo = []
    index = 0
    while index < len(self.infoFieldNames):
      try:
        meltedInfo.append(infoHash[ self.infoFieldNames[index] ])
      except KeyError as exc:
        # This allele does not have a value for the current value in "infoFieldNames"
        meltedInfo.append('.')
      index += 1
    return meltedInfo

    
if __name__ == '__main__':
  try:
    x = ExomeBurden()
    x.read_vcf_file()
  except BurdenException as exc:
    print exc.get_message()
