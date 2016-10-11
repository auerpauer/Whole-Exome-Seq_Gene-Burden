#!/usr/bin/python

from optparse import OptionParser
import sys, os, re
import traceback
#import pdb

# INTERNAL EXCEPTION CLASS
# This is used to pass error message up the stack to a context that has more info to understand and report the error
# The "stack" holds the methods that raise and pass this object up the stack
# New methods are appended to the end, but when reporting the stack, this order is reveresed, to give highest level method first
class BurdenException(Exception):
  # error codes
  input_error = 1010
  parsing_error = 2010
  filtered_out = 3000 # All error codes larger in value than this one represent a record that is filtered out, but not erroneous
  filter_freq_too_high = 3010 # All 301X error codes represent a variant that occurs too often to be of interests to us
  filter_freq_esp6500_too_high = 3011
  filter_freq_1000g_too_high = 3012
  filter_freq_exac_too_high = 3013
  filter_on_refgene = 3020 # represents a record that did not pass filters on refGene values
  filter_nonexonic_nonsplicing = 3021
  filter_nonframshift_indel = 3022
  filter_synonymous_unknown = 3023
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
  # The value matching each key is a list which contains three more lists for missense, missense_damaging & damaging variant types
  # Each of these inner lists contains two more lists which contain names of samples that have heterozygous or homozygous genotypes
  # When the number of genes stored in 'countingSets' passes a threshold, the oldest gene in the hash is processed and written to file
  # This happens without the calling class (ExomeBurden) knowing that a gene is being processed and deleted from this class
  missense = 0
  missense_damaging = 1
  damaging = 2
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
    self.tableFH.write('x\t')
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
    self.tableFH.write('Hom_Missense\tHom_Deleterious\tHom_LoF\n')
    return

    ### add_mutation ###
    # called by external class
  def add_mutation(self, geneName, variantType, genotype, sampleName):
    if geneName not in self.countingSets.keys():
      self.countingSets[geneName] = [([],[]), ([],[]), ([],[])] # missense:(het, hom); missense_damaging:(het,hom); damaging:(het,hom)
      self.addedGeneNames.append(geneName)
    self.countingSets[geneName][variantType][genotype].append(sampleName)
    if len(self.countingSets.keys()) > self.geneCountThreshold:
      self.process_oldest_gene()
    return

    ### process_oldest_gene ###
    # called by MutationCount.add_mutation() and MutationCount.flush_genes()
    # This consolidates all the mutations for a specific gene and writes them to their respective table files
    # IDEALLY, this method would be called as a separate thread, but I do not have time to design & test that implementation
    # The main logic here is to count the variants and the samples and figure out compound heterozygotes from single heterozygotes
    #   for a given sample name
    # If a sample contributes only one heterozygote, then we count that single variant as 'het'
    # However, if a sample contributes more than one heterozygote to a gene, it is considered a compound heterozygote
  def process_oldest_gene(self):
    oldestGene = self.addedGeneNames.pop(0) # pop the first name off the list
    #if oldestGene == 'PRAMEF18,PRAMEF22': pdb.set_trace() # FDFT1, PABPC3, MPRIP
    variantCount = ([0,0,0], [0,0,0], [0,0,0]) # het:(missense, missense-damaging, damaging), cmpdHet(..), hom(..)
    sampleCount = ([0,0,0], [0,0,0], [0,0,0]) # same structure as 'variantCount'
    varTypeIndex = 0
    sampleNamesCounted = []
    while varTypeIndex < len(self.countingSets[oldestGene]):
      hetSamples = self.countingSets[oldestGene][varTypeIndex][self.het]
      uniqSampleNames = list(set(hetSamples))
      for name in uniqSampleNames:
        # count the number of times this name occurs in 'hetSamples'
        occurrences = self.count_name_in_list(name, hetSamples)
        variantCount[0][varTypeIndex] += occurrences # "normal" heterozygous
        remainder = occurrences % 2
        if remainder > 0:
          occurrences -= 1
        variantCount[1][varTypeIndex] += occurrences # compound heterozygous
        # if this sample hasn't already been counted, increment the count and the name to the list; otherwise, ignore it
        if name not in sampleNamesCounted:
          sampleCount[0][varTypeIndex] += 1
          sampleNamesCounted.append(name)
      del hetSamples
      homSamples = self.countingSets[oldestGene][varTypeIndex][self.hom]
      uniqSampleNames = list(set(homSamples))
      variantCount[2][varTypeIndex] = len(homSamples) * 2
      for name in uniqSampleNames:
        if name not in sampleNamesCounted:
          sampleCount[2][varTypeIndex] += 1 # homozygous
          sampleNamesCounted.append(name)
      varTypeIndex += 1
    if sum(sampleCount[0]) > 0 or sum(sampleCount[1]) > 0 or sum(sampleCount[2]) > 0: 
      self.write_to_tables(oldestGene, variantCount, sampleCount) # write every gene
    self.countingSets[oldestGene] = None # remove data stored at key:geneName
    self.countingSets.pop(oldestGene, None) # remove key:geneName
    return

    ### count_name_in_list ###
  def count_name_in_list(self, sampleName, nameList):
    count = 0
    for currentName in nameList:
      if sampleName == currentName:
        count += 1
    return count

  ### write_to_tables ###
  # called by process_oldest_gene()
  # This method is insanely ugly, I did not have the time to write it in a modular fashion
  # varCount = het(missense, missense-damaging, damaging), cmpdHet(..), hom(..)
  # smpCount has the same layout as varCount
  # Mis -> missense; MisDmg -> missense-damaging; Dmg -> damaging
  def write_to_tables(self, geneName, varCount, smpCount):
    #if 'MPRIP' == geneName: pdb.set_trace() # PABPC3 # debug code for interlaced genes in the genome

    # compute and write variant allele counts
    hetDmg = varCount[0][2]
    hetMisDmg = hetDmg + varCount[0][1]
    hetMis = hetMisDmg + varCount[0][0]

    homPlusCmpdHetDmg = varCount[1][2] + varCount[2][2]
    homPlusCmpdHetMisDmg = homPlusCmpdHetDmg + varCount[1][1] + varCount[2][1]
    homPlusCmpdHetMis = homPlusCmpdHetMisDmg + varCount[1][0] + varCount[2][0]

    homDmg = varCount[2][2]
    homMisDmg = homDmg + varCount[2][1]
    homMis = homMisDmg + varCount[2][0]

    totDmg = hetDmg + homPlusCmpdHetDmg
    totMisDmg = hetMisDmg + homPlusCmpdHetMisDmg
    totMis = hetMis + homPlusCmpdHetMis

    if self.excel:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\tx\t'.\
                         format(geneName,
                                totMis, totMisDmg, totDmg,
                                hetMis, hetMisDmg, hetDmg,
                                homPlusCmpdHetMis, homPlusCmpdHetMisDmg, homPlusCmpdHetDmg,
                                homMis, homMisDmg, homDmg))
    else:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t'.\
                         format(geneName,
                                totMis, totMisDmg, totDmg,
                                hetMis, hetMisDmg, hetDmg,
                                homPlusCmpdHetMis, homPlusCmpdHetMisDmg, homPlusCmpdHetDmg,
                                homMis, homMisDmg, homDmg))

    # compute and write sample variant counts
    hetDmg = smpCount[0][2]
    hetMisDmg = hetDmg + smpCount[0][1]
    hetMis = hetMisDmg + smpCount[0][0]

    homPlusCmpdHetDmg = smpCount[1][2] + smpCount[2][2]
    homPlusCmpdHetMisDmg = homPlusCmpdHetDmg + smpCount[1][1] + smpCount[2][1]
    homPlusCmpdHetMis = homPlusCmpdHetMisDmg + smpCount[1][0] + smpCount[2][0]

    homDmg = smpCount[2][2]
    homMisDmg = homDmg + smpCount[2][1]
    homMis = homMisDmg + smpCount[2][0]

    totDmg = hetDmg + homPlusCmpdHetDmg
    totMisDmg = hetMisDmg + homPlusCmpdHetMisDmg
    totMis = hetMis + homPlusCmpdHetMis

    if self.excel:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.\
                         format(totMis, totMisDmg, totDmg,
                                hetMis, hetMisDmg, hetDmg,
                                homPlusCmpdHetMis, homPlusCmpdHetMisDmg, homPlusCmpdHetDmg,
                                homMis, homMisDmg, homDmg))
    else:
      self.tableFH.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n'.\
                         format(totMis, totMisDmg, totDmg,
                                hetMis, hetMisDmg, hetDmg,
                                homPlusCmpdHetMis, homPlusCmpdHetMisDmg, homPlusCmpdHetDmg,
                                homMis, homMisDmg, homDmg))
      
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
# Each line in the VCF file is processed independently. The only exception to this is counting variations that heterozygous
# for compound variant tracking.
#
# Most of the information we care about resides in the INFO field of the annotation. This data is packed down by a few levels,
# which is extracted with the "<string>.split(<delimiter>)" method. This yields 'key=value' pairs, which are also split.
# The 'key=value' pairs are stored in the hash "infoHash". The list "screenedInfoField" is a list of 'key' values that are
# used for filtering the variant. Thus, we store all values in the INFO field, but only use a few to determine if the variant is
# worth further consideration.
#
# Because this is bioinformatics, there are exceptions to how the data is stored or labeled in the VCF file that require
# extra code to deal with. The most egregious is that there are two different names for ExAC data, used by different VCF
# files. In other words, some VCF files use one string throughout the entire file, other VCF files use the other string
# throughout the entire file.
# To allow for this, "screenedInfoFields" is allowed to have a list instead of a string in it's elements. If a list is detected,
# the elements in the list are matched against values in the VCF file. One a match is found, that value is used to overwrite
# the list. The code in method "screened_field()" performs this.
#
# Another example of exceptional data handling is the version strings appended to the other two SNV population databases.
# To deal with this, "screenedInfoField" has a string for each of these databases without any version information. These values
# are matched up against the value in the VCF file. When a match is found THE VALUE IN "screenedInfoFields" IS USED to allow
# for smooth processing of the data values for these databases.
#
# Another wrinkle to the filtering of variants is that those variants that pass all filters are to have their records "melted".
# This means that, for each sample, all data is written out in tab-delimited format. This creates a very great redundancy in the
# resulting data table. This is why ALL VALUES in the INFO field are saved to in "infoHash". This allows all these fields to be
# written to the melted file.

class ExomeBurden:
  # Values that users can specify for the functional deleterious filter to be used
  allowedFuncDelFilters = ('cadd', 'metaSVM')
  allowedPopFilterNames = ('esp6500', '1000g', 'exac')
  desiredFormatFields = ['GT', 'AD', 'DP', 'GQ', 'PL'] # melted table output formatting of FORMAT field

  # These hold the command line paramters specified by the user at runtime
  desiredPopFreq = None
  vcfFH = None
  funcDel = None
  excel = False
  outputDir = None
  removeNonFrameShiftInDel = False
  doNotFilter = []
  requiredDepth = 8 # default required read depth value

  # Values used to melt VCF records
  meltFH = None
  infoFieldNames = []
  formatFieldNames = []
  
  # Please note that all arrays that count variants are in the same order as this one
  allowedVariantTypes = ('missense', 'missense_damaging', 'damaging')
  # hardcoded indexes into 'mutationCount'
  het = 0
  hom = 1
  # genotypes from sample data that we do not want to process (they don't contain any trackable variants)
  unwantedGenotypes = ('0/0', '0|0', './.', '.|.')

  ### init ###
  def __init__(self):
    vcfFileName = self.parse_cmdline_arguments()
 
    # YOU CANNOT CHANGE THE ORDER OF THESE STRINGS they are used as keys into 'self.infoHash' when filtering variants
    self.screenedInfoFields = ['ExonicFunc.refGene', 'Gene.refGene', 'Func.refGene', 'esp6500',
                               '1000g', ['ExAC_ALL','exac03']] #, 'CADD_phred']
    if 'metaSVM' == self.funcDel:
      self.screenedInfoFields.append('MetaSVM_pred')
    elif 'cadd' == self.funcDel:
      self.screenedInfoFields.append('CADD_phred')
    
    # create output file name
    filename = os.path.basename(vcfFileName)
    filename = os.path.splitext(filename)[0]
    countTableName = '{0}_{1}_{2}_counts_table.tsv'.format(filename, self.desiredPopFreq, self.funcDel)
    self.mutCount = MutationCount(countTableName, self.outputDir, self.excel)
    meltedDataFile = '{0}_{1}_{2}_variant.table'.format(filename, self.desiredPopFreq, self.funcDel)
    meltedDataFile = os.path.join(self.outputDir, meltedDataFile)
    self.meltFH = open(meltedDataFile, 'w')
    return

  def __del__(self):
    if self.vcfFH: self.vcfFH.close()
    if self.meltFH: self.meltFH.close()
    return

  ### parse_cmdline_arguments ###
  # see "print_short_help()" for an explanation of the command line options
  def parse_cmdline_arguments(self):
    if 1 == len(sys.argv):
      self.print_short_help()
    if 7 > len(sys.argv) or 14 < len(sys.argv):
      sys.stderr.write('WARNING - Incorrect number of arguments\n')
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
      elif '-e' == sys.argv[argIndex] or '--exempt-database' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"exempt-database\" argument\n')
          self.print_short_help()
        self.doNotFilter = sys.argv[argIndex+1].split(',')
        argIndex += 2
      elif '-c' == sys.argv[argIndex] or '--coverage' == sys.argv[argIndex]:
        if argIndex + 1 >= len(sys.argv):
          sys.stderr.write('ERROR - missing value for \"coverage\" argument\n')
          self.print_short_help()
        self.requiredDepth = int(sys.argv[argIndex+1])
        argIndex += 2
      else:
        sys.stderr.write('ERROR - illegal format for command line parameters\n')
        self.print_short_help()

    if None == self.desiredPopFreq or None == self.vcfFH or None == self.funcDel:
      sys.stderr.write('ERROR - frequency, VCF file and functional deleteriousness must be specified\n')
      self.print_short_help()

    for name in self.doNotFilter:
      if name not in self.allowedPopFilterNames:
        sys.stderr.write('ERROR - illegal population database name: {0}\n'.format(name))
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
      fh = open(filename, 'r')
    except IOError as exc:
      sys.stderr.write('ERROR - unable to open VCF file: {0}\n'.format(filename))
      self.print_short_help()
    self.vcfFH = fh
    return

  ### parse_functional_deleteriousness ###
  def parse_functional_deleteriousness_filter(self, filterName):
    for value in self.allowedFuncDelFilters:
      if value.lower() == filterName.lower():
        self.funcDel = filterName
        return 
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
    sys.stderr.write('usage: {0} parameters [options] \n\n'.format(sys.argv[0]))
    sys.stderr.write('Parameters:\n')
    sys.stderr.write('\t-f | --frequency\tmax desired frequency of variants in population databases\n')
    sys.stderr.write('\t-v | --vcf\t\tVCF file to read in containing variants to be analyzed\n')
    sys.stderr.write('\t-d | --funcDel\t\tmethod to use for filtering functional deleteriousness: CADD or metaSVM\n')
    sys.stderr.write('Options\n')
    sys.stderr.write('\t-o | --outputDir\toutput directory to write files; default: current directory\n')
    sys.stderr.write('\t-x | --excel\t\tFLAG: use to cause more header information for display in spreadsheet\n')
    sys.stderr.write('\t-s | --shift\t\tFLAG: use to specify that nonframeshift indels should be filtered out\n')
    sys.stderr.write('\t-c | --coverage\t\trequired read depth for a variant to be accepted; default: 8\n')
    sys.stderr.write('\t-e | --exempt-database\tcomma-delimited list of population databases to NOT use as filters: {0}\n'.\
                     format(','.join(self.allowedPopFilterNames)))
    sys.stderr.flush()
    sys.exit(1)
    return
    
  ### read_vcf_file ###
  # called by instantiated object (see __main__ below)
  def read_vcf_file(self):
    headerLine = self.parse_vcf_meta_data()
    if 0 == len(self.infoFieldNames) or 0 == len(self.formatFieldNames):
      errMsg = 'ERROR - unable to find INFO or FORMAT meta-data...does the VCF file have a meta-data header?'
      raise BurdenException(BurdenException.parsing_error, errMsg, 'read_vcf_file')
    if not headerLine:
      errMsg = 'ERROR - unable to find header line in VCF file\n'
      raise BurdenException(BurdenException.parsing_error, errMsg, 'read_vcf_file')
    self.sampleNames = headerLine.split('\t')[9:]
    self.write_melted_header_line()
    self.sampleNames = headerLine.split('\t')[9:]
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
    line = self.vcfFH.readline().rstrip()
    while line:
      if '#CHROM' == line[0:6]:
        return line # return header line
      if '##INFO' == line[0:6]:
        identifier = re.search('##INFO=<ID=(.*?),Number', line)
        if not identifier:
          errMsg = 'ERROR - Unable to parse info name from meta data line\n{0}\n'.format(line)
          raise BurdenException(BurdenException.parsing_error, errMsg, 'parse_vcf_meta_data')
        self.infoFieldNames.append(identifier.group(1))
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
    self.meltFH.write( 'SAMPLE\t{0}\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER{1}'.
                       format('\t'.join(self.desiredFormatFields), '\t'.join(self.infoFieldNames)) )
    for name in self.formatFieldNames:
      self.meltFH.write('\t{0}'.format(name))
    self.meltFH.write('\n')
    self.meltFH.flush()
    return
    
  ### read_variants ###
  # called by read_vcf_file()
  # This is the high-level method that performs most of the work on each record
  # Due to the size of most VCF files, we need to read line-by-line, as opposed to storing the file contents in memory
  # Determine the position of the GQ value in the FORMAT field
  # First, the variant is passed through filters and skipped if it doesn't pass the criteria
  # Second, the variant is classified into one of three categories: Damaging, Damaging Missense or Missense
  # Third, the number of samples with this variant are counted; this is performed in concert with the class MutationCount
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
      columns = line.split('\t')
      #CHROM  POS       ID      REF       ALT       QUAL    FILTER       INFO       FORMAT
      (chrom, position, ID, refAllele, altAlleles, quality, filterField, infoField, formatField) = columns[0:9]
      samples = columns[9:]
      try:
        self.parse_info_field(infoField)
        if '.' == self.infoHash['Gene.refGene']: # skip records that do not have a gene name
          line = self.vcfFH.readline().rstrip()
          continue
        self.filter_on_database_population_frequencies()
        self.filter_on_refgene_values()
        self.count_genotypes(self.infoHash['Gene.refGene'], samples)
        # PARSE OUT RECORD AND WRITE TO FILE
        self.melt_vcf_record(line)
      except BurdenException as exc:
        if exc.errCode > BurdenException.filtered_out:
          pass # skip this record and read the next one
        else:
          if BurdenException.parsing_error == exc.errCode:
            exc.add_to_stack('read_variants')
            exc.add_to_message('in record with Chrom:{0} and Position:{1}'.format(chrom, position))
          raise exc
      line = self.vcfFH.readline().rstrip()
    return

  ### parse_info_field ###
  # called by read_variants()
  # This method does NOT check if a field is overwritten
  # That is, if the same field occurs twice in a single INFO field, the second value will be used without errors or warnings
  # ja623 May 18, 2016
  # To allow for newer versions of field names (1000g2012apr vs 1000g2014oct)
  # better string matching must be applied to some field names.
  # Currently these are: 1000 Genomes "1000g"; NHLBI Exome Sequencing Project "esp6500"
  def parse_info_field(self, infoField):
    self.infoHash = {}
    keyValuePairs = infoField.split(';')
    index = 0
    while index < len(keyValuePairs):
      if "=" in keyValuePairs[index]:
        (name, value) = keyValuePairs[index].split('=')
        key = self.screened_field(name)
        if key:
          self.infoHash[key] = value
        else:
          self.infoHash[name] = value
      index += 1
    try:
      self.verify_info_fields()
    except BurdenException as exc:
      exc.add_to_stack('parse_info_field')
      raise exc
    return

  ### screened_field ###
  # ja623 May 18, 2016
  # Allow more general matching of INFO field names to allow matching of different versions of databases.
  # If one of the screened field names is contained in 'currentField', return the value that matched against 'currentField'.
  # This way, we dictate what keys are used in 'self.infoHash', when matched later in 'filter_on_database_population_frequencies()'.
  # This also frees us from using any appended strings to indicate database version.
  #
  # There are also elements in "screenedInfoFields" that are actually a list of values. This represents the fact that some
  # VCF files have different strings identifying the same data. Once the first record is read, the value used by the current
  # VCF file is used to replace the element containing the list.
  def screened_field(self, currentField):
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
        #matchObj = None # unneccesary
        degenerateFieldNames = self.screenedInfoFields[index]
        for fieldName in degenerateFieldNames:
          matchObj = re.search(fieldName, currentField)
          if matchObj:
            self.screenedInfoFields[index] = fieldName # alter the value to reflect what the current VCF contains
            if 'exac' in self.doNotFilter and 'exac' in fieldName.lower():
              position = self.doNotFilter.index('exac')
              self.doNotFilter[position] = fieldName
            return fieldName
      index += 1
    return False

  ### verify_info_fields ###
  # called by parse_info_field()
  # Check that all the fields we want to filter on are present
  # If you want to filter on other fields, you need to edit screenedInfoFields at the top of this class definition
  def verify_info_fields(self):
    missingFields = []
    index = 0
    while index < len(self.screenedInfoFields):
      fieldName = self.screenedInfoFields[index]
      if fieldName not in self.infoHash.keys():
        if not isinstance(fieldName, basestring):
          errString = 'ERROR - Unable to find one of the fields names in : {0}\n'.format(','.join(fieldName))
          raise BurdenException(BurdenException.parsing_error, errString, 'verify_info_fields')
        missingFields.append(fieldName)
      index += 1
    if missingFields:
      errString = 'ERROR - The following fields were not present in the INFO field: {0}\n'.format(','.join(missingFields))
      raise BurdenException(BurdenException.parsing_error, errString, 'verify_info_fields')
    return

  ### filter_on_database_population_frequencies ###
  # called by read_variants()
  def filter_on_database_population_frequencies(self):
    #if '.' == self.infoHash['esp6500siv2_all']: self.infoHash['esp6500siv2_all'] = 0
    #if '.' == self.infoHash['1000g2014oct_all']: self.infoHash['1000g2014oct_all'] = 0
    #if '.' == self.infoHash['ExAC_ALL']: self.infoHash['ExAC_ALL'] = 0
    if '.' == self.infoHash[self.screenedInfoFields[3]]:
      self.infoHash[self.screenedInfoFields[3]] = 0
    else:
      self.infoHash[self.screenedInfoFields[3]] = float(self.infoHash[self.screenedInfoFields[3]])
    if '.' == self.infoHash[self.screenedInfoFields[4]]:
      self.infoHash[self.screenedInfoFields[4]] = 0
    else:
      self.infoHash[self.screenedInfoFields[4]] = float(self.infoHash[self.screenedInfoFields[4]])
    if '.' == self.infoHash[self.screenedInfoFields[5]]:
      self.infoHash[self.screenedInfoFields[5]] = 0
    else:
      self.infoHash[self.screenedInfoFields[5]] = float(self.infoHash[self.screenedInfoFields[5]])

    if self.screenedInfoFields[3] not in self.doNotFilter:
      if self.infoHash[self.screenedInfoFields[3]] > self.desiredPopFreq:
        raise BurdenException(BurdenException.filter_freq_esp6500_too_high, '-', 'filter_on_database_population_frequencies')
    if self.screenedInfoFields[4] not in self.doNotFilter:
      if self.infoHash[self.screenedInfoFields[4]] > self.desiredPopFreq:
        raise BurdenException(BurdenException.filter_freq_1000g_too_high, '-', 'filter_on_database_population_frequencies')
    if self.screenedInfoFields[5] not in self.doNotFilter:
      if self.infoHash[self.screenedInfoFields[5]] > self.desiredPopFreq:
        raise BurdenException(BurdenException.filter_freq_exac_too_high, '-', 'filter_on_database_population_frequencies')
    return

  ### filter_on_refgene_values ###
  # called by read_variants()
  def filter_on_refgene_values(self):
    if 'exonic' != self.infoHash['Func.refGene'] and 'splicing' != self.infoHash['Func.refGene'] and \
       'exonic\x3bsplicing' != self.infoHash['Func.refGene']:
      raise BurdenException(BurdenException.filter_nonexonic_nonsplicing, '-', 'filter_on_refgene_values')
    if self.removeNonFrameShiftInDel:
      if 'nonframeshift_insertion' == self.infoHash['ExonicFunc.refGene'] or \
         'nonframeshift_deletion'  == self.infoHash['ExonicFunc.refGene']:
        raise BurdenException(BurdenException.filter_nonframshift_indel, '-', 'filter_on_refgene_values')
    if 'synonymous_SNV' == self.infoHash['ExonicFunc.refGene'] or 'unknown' == self.infoHash['ExonicFunc.refGene']:
      raise BurdenException(BurdenException.filter_synonymous_unknown, '-', 'filter_on_refgene_values')
    return

  ### classify_variant ###
  # called by read_variants()
  def classify_variant(self):
    dScore = self.determine_deleteriousness()
    if ('splicing' in self.infoHash['Func.refGene'] or
        'frameshift' == self.infoHash['ExonicFunc.refGene'][0:10] or
        'stop' in self.infoHash['ExonicFunc.refGene']):
      variantType = self.mutCount.damaging # 2
    elif (('nonsynonymous_SNV' == self.infoHash['ExonicFunc.refGene'] and 'deleterious' == dScore) or
          'nonframeshift' in self.infoHash['ExonicFunc.refGene']):
      variantType = self.mutCount.missense_damaging # 1
    else:
      variantType = self.mutCount.missense # 0
    return variantType

  ### determine_deleteriousness ###
  def determine_deleteriousness(self):
    dScore = None
    if 'metaSVM' == self.funcDel:
      if '.' == self.infoHash['MetaSVM_pred']: self.infoHash['MetaSVM_pred'] = 'T'
      if 'T' == self.infoHash['MetaSVM_pred']:
        dScore = 'tolerated'
      else:
        dScore = 'deleterious'
    if 'cadd' == self.funcDel:
      if '.' == self.infoHash['CADD_phred']: self.infoHash['CADD_phred'] = 0
      caddPhred = round(float(self.infoHash['CADD_phred']))
      if 15 > caddPhred:
        dScore == 'tolerated'
      else:
        dScore = 'deleterious'
    return dScore

  ### count_genotypes ###
  # called by read_variants()
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
    variantType = self.classify_variant()
    index = 0
    while index < len(samples):
      currentSample = samples[index]
       # filter out low quality genotypes
      currentGenotype = currentSample.split(':')[0]
      if currentGenotype in self.unwantedGenotypes:
        index += 1
        continue
      if 20 >= int(currentSample.split(':')[self.GQindex]):
        index += 1
        continue

      # filter out low read depth variants
      currentDP = currentSample.split(':')[self.DPindex]
      if '.' == currentDP:
        currentDP = 0
      else:
        currentDP = int(currentDP)
      if currentDP < self.requiredDepth:
        index += 1
        continue
      
      if '/' in currentGenotype:
        (plusStrand, minusStrand) = currentGenotype.split('/')
      elif '|' in currentGenotype:
        (plusStrand, minusStrand) = currentGenotype.split('|')
      else:
        errMsg = 'ERROR - Unable to parse genotype from sample data: {0}'.format(currentSample)
        raise BurdenException(BurdenException.parsing_error, errMsg, 'count_genotypes')
      sampleName = self.sampleNames[index]
      if plusStrand == minusStrand:
        self.mutCount.add_mutation(geneName, variantType, self.mutCount.hom, sampleName)
      else:
        self.mutCount.add_mutation(geneName, variantType, self.mutCount.het, sampleName)
      index += 1
    return

  ### melt_vcf_record ###
  # "self.infoHash" is already populated!
  def melt_vcf_record(self, line):
    newRecord = None
    columns = line.split('\t')
    #CHROM  POS       ID      REF       ALT       QUAL    FILTER       INFO       FORMAT
    (chrom, position, ID, refAllele, altAlleles, quality, filterField, infoField, formatField) = columns[0:9]
    samples = columns[9:]
    newRecord = [chrom, position, ID, refAllele, altAlleles, quality, filterField] # 7 fields populated
    newRecord.extend( self.melt_info_field(infoField) )
    index = 0
    while index < len(self.sampleNames):
      self.write_melted_record(newRecord, formatField, self.sampleNames[index], samples[index])
      index += 1
    return

  ### melt_info_field ###
  def melt_info_field(self, infoField):
    #meltedInfo = ['.'] * len(self.infoFieldNames)
    meltedInfo = []
    index = 0
    while index < len(self.infoFieldNames):
      try:
        meltedInfo.append(self.infoHash[ self.infoFieldNames[index] ])
      except KeyError as exc:
        # This record does not have a value for the current value in "infoFieldNames"
        meltedInfo.append('.')
      index += 1
    return meltedInfo

  ### write_melted_record ###
  def write_melted_record(self, newRecord, formatField, sampleName, sample):
    sampleData = sample.split(':')
    if sampleData[0] in self.unwantedGenotypes:
    #if './.' == sampleData[0] or '.|.' == sampleData[0]:
      return # do not print data for samples that have no genotype information
    meltedSample = []
    formatHeader = formatField.split(':')
    index = 0
    #pdb.set_trace()
    while index < len(self.desiredFormatFields):
      try:
        meltedSample.append( sampleData[formatHeader.index(self.desiredFormatFields[index])] )
      except ValueError as exc:
        if 'PL' == self.desiredFormatFields[index]:
          meltedSample.append('.')
          pass
        else:
          sys.stderr.write('ERROR - format field {0} does not appear in the FORMAT field {0}\n'.format(formatField))
          sys.stderr.flush()
          sys.exit(1)
      index += 1
    self.meltFH.write( '{0}\t{1}\t{2}\n'.format(sampleName, '\t'.join(meltedSample), '\t'.join('{0}'.format(n) for n in newRecord)) )
    self.meltFH.flush()
    return


if __name__ == '__main__':
  try:
    x = ExomeBurden()
    x.read_vcf_file()
  except BurdenException as exc:
    print exc.get_message()
