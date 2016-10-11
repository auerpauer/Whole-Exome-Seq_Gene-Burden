#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: burden_script_fisher.pl
#
#        USAGE: ./burden_script_fisher.pl  
#
#  DESCRIPTION: :wq
#
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Francesc Lopez (FL), francesc.lopez@gmail.com
# ORGANIZATION: YCGA
#      VERSION: 1.0
#      CREATED: 06/19/15 13:56:08
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use autodie;
use Data::Dumper;
use Test::More;
use Carp;
use Text::NSP::Measures::2D::Fisher::right;

my $tab_rx = qr/\t/xm;
my $USAGE = "\nPlease, provide the counts for cases, number of cases, file and control file, number of controls.\n" 
          . "\tEx: $0 table_counts.cases.0.001.txt 100 table_counts.controls.0.001.txt 2700\n\n";

my $cases_infile    = shift @ARGV or croak $USAGE;
my $num_cases       = shift @ARGV or croak $USAGE;
my $controls_infile = shift @ARGV or croak $USAGE;
my $num_controls    = shift @ARGV or croak $USAGE;

my ($cases_header, $cases_count_for_ref) = get_counts($cases_infile);
my ($controls_header, $controls_count_for_ref) = get_counts($controls_infile);


my %cases_count_for = %$cases_count_for_ref;
my %controls_count_for = %$controls_count_for_ref;

my @header_cols = split $tab_rx, $cases_header;

my $prefix1 = $cases_infile;
$prefix1 =~ s/\.tsv$//;
#my $prefix2 = $controls_infile;
#$prefix2 =~ s/\.tsv$//;


my $outfile = "$prefix1.burden_analysis_fisher_table.txt";
print "Printing output in [$outfile]\n";
open my $OUT, '>', $outfile;


#Printing header
local $" = "\tp-value\t"; # $" is printed between values of an expanded array when printed out
# ja623 28/09/2015 Fixed header to match printed data columns
#print {$OUT} "Gene\t@header_cols[1..$#header_cols]\tp-value\n";
print {$OUT} "Gene";
my $index = 1;
while ($index < @header_cols) {
		print {$OUT} "\t$header_cols[$index](case)\t(control)\tp-value";
		$index += 1;
}
print {$OUT} "\tX"; # Uncomment this and look at the error msg, this is why perl sucks
$index = 1;
while ($index < @header_cols) {
		print {$OUT} "\t$header_cols[$index](case)\t(control)\tp-value";
		$index += 1;
}
print {$OUT} "\n";
# ja623 28/09/2015 end

foreach my $gene (sort keys %cases_count_for) {
    print {$OUT} "$gene";
    my @controls_values = $controls_count_for{$gene} ? @{$controls_count_for{$gene}} 
                        : (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0); # jma623 28/09/2015 Extend array to 24 elements

    #Fisher for alleles section
    foreach my $col (0..11) {
        my $right_value = calculateStatistic(n11=>$cases_count_for{$gene}[$col],
                                            n1p=>$cases_count_for{$gene}[$col]+$controls_values[$col],
                                            np1=>$num_cases*2,
                                            npp=>$num_controls*2);
				if (!(length $right_value)) { $right_value = "NA"; } # ja623 28/09/2015 prevent no value printing for some p-values
        print {$OUT} "\t$cases_count_for{$gene}[$col]\t$controls_values[$col]\t$right_value";
    }

		print {$OUT} "\tX"; # Uncomment this and look at the error msg, this is why perl sucks

    #Fisher for individuals section
    foreach my $col (12..23) {
         my $right_value = calculateStatistic(n11=>$cases_count_for{$gene}[$col],
                                            n1p=>$cases_count_for{$gene}[$col]+$controls_values[$col],
                                            np1=>$num_cases,
                                            npp=>$num_cases+$num_controls);
				if (!(length $right_value)) { $right_value = "NA"; } # ja623 28/09/2015 prevent no value printing for some p-values
        print {$OUT} "\t$cases_count_for{$gene}[$col]\t$controls_values[$col]\t$right_value";

    }
    print {$OUT} "\n";
}


exit;

################################################################################
#
# SUBROUTINES
#
################################################################################

sub get_counts {
    my ($infile) = @_;
    my %counts_for;
    
    print "\nGetting counts for [$infile]\n";
    open my $IN, '<', $infile;
    my $header = readline $IN;
    chomp $header;
    print "\tHeader::: [$header]\n";

    my $total_genes_processed = 0;
    while (my $line = readline $IN) {
        chomp $line;
        ++$total_genes_processed;
        my @cols = split $tab_rx, $line;
        $counts_for{$cols[0]} = [@cols[1..$#cols]]; 
    }

    print "\tTotal genes processed ::: $total_genes_processed\n";
    return $header, \%counts_for;
}
