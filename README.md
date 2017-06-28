# Whole Exome Sequencing: Gene Burden Overview
Which genes have the most harmful variants in your set of whole exome sequencing samples?

### WHAT'S GOING ON HERE?
The general idea of "gene burden" is finding those genes that have a significant number of harmful variants among your samples. One way to do this is to compare the results from your samples against ethnic controls. This way, any genes that are highly mutated in the proper ethnic controls will not trigger a false positive in your result set.


So, what is a harmful variant? This analysis pipeline makes the following assumptions of what a harmful variant is within the constraints of VCF annotation:
1. The variant frequency in the ExAC databases is less than or equal to a given threshold. You decide what the threshold is at run-time. Suggestions: 5e-05 for dominant (homozygous) variants; 1e-04 for recessive (heterozygous) variants
2. RefSeq annotation must annotate the variant as being either exonic or splicing.
3. RefSeq annotation must NOT annotate the variant as synonymous or unknown.
4. Optionally, RefSeq annotation can be used to filter out non-frameshift deletions. You decide at run-time if you want this filtering turned on.

If all of these filters are passed, then the variant is classified according to the following criteria:
1. Loss-of-Function: RefSeq annotates the variant as "stopgain", "splicing" or "frameshift".
2. Deleterious: One of the following three conditions must be true.
*   The variant is a nonframeshift insertion or deletion.
*   The variant causes a "stoploss".
*   The variant is a non-synonymous SNV and either:  
 	 a) CADD score is 15 or greater or  
	 b) The score of the other chosen functional deleteriousness method is "D" (for deleterious).  
  You choose between CADD or another scoring method at runtime.
3. Missense: This is the default classification, if the other two criteria are not met. Remember, we have already filtered out all variants that we are not interested in. All remaining variants must receive a classification.

The variants are further classified based on whether the variant is homozygous, heterozygous, or compound heterozygous. Compound hets occur when the same sample contributes more than one heterozygous variant to the same gene.

|Missense           |Deleterious           |Loss of Function|
|-------------------|----------------------|----------------|
|Total_Missense     |Total_Deleterious     |Total_LoF       |
|Het_Missense       |Het_Deleterious       |Het_LoF         |
|HomComHet_Missense |HomComHet_Deleterious |HomComHet_LoF   |
|Hom_Missense       |Hom_Deleterious       |Hom_LoF         |

Where
* Het: heterozygous mutation
* HomComHet: compound heterozygous mutations
* Hom: homozygous mutation


These variants are counted twice in two separate, but related ways:  
1. The total number of variants in a gene.  
2. The number of samples contributing one or more variants in a gene.  

Just to be clear, if a one or more samples contains more than one variant for a gene, these number are not going to be the same.

The results of these filters and classifications are saved to a tab-delimited file that has the counts for all the above variant types.
You must also run these filters and classifications on a set of controls that share the same ethnic background as your samples.
For each gene, the variant counts of your samples and ethnic controls are then compared using a right-tailed Fisher's exact test.
This generates a p-value for each variant classification.
These p-values are written to a tab-delimited file, which can be imported into a spreadsheet program.
Then you can sort by p-value to find the most significantly mutated genes under the variant category of your choosing.
The most commonly used categories are Total_LoF, Total_Deleterious and Total_Missense.

The annotations used for filtering and mutation classification are better explained at: http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#output-file-1-refseq-gene-annotation

In addition to using the output file that holds the counts of mutations, another file is produced that describes information about each significant variant. This file is intended for further programmatic analysis, but can be viewed in a spreadsheet program, if you wish. This second file is very large and holds redundant data, compared to the VCF file. However, we deemed it necessary to create this file for further analysis.


# HOW DO I USE THESE PROGRAMS?
You must first have an annotated VCF file. The annotations must include ExAC frequency, RefSeq annotation and the deleteriousness score you choose to use (such as CADD).
The VCF can have all the possible deleteriousness scores, if you wish.
I have only used VCFs generated by the GATK pipeline and annotated with ANNOVAR.
Your VCF file MUST include the meta-data header generated by ANNOVAR.

You need to run two programs in order: exome_burden_script.py (once on your samples and once on your controls). You then need to run gene_burden_fisher_exact.py, which takes the mutation counts for both your samples and their ethnic control.
You must be running python version 2.7.5 or higher. These programs have NOT been written for Python 3.

Assumptions:
1. The programs in this repository are located on your machine at /home/user/programs
2. Your VCF files are located on your machine at /home/user/data
3. Results of these analysis will be stored at /home/user/results

## Filter and Classify Mutations: exome_burden_script.py
```text
$ /home/user/programs/exome_burden_script.py
usage: exome_burden_script.py parameters [options]

Parameters:
        -f | --frequency        max desired frequency of variants in population databases
        -v | --vcf              VCF file to read in containing variants to be analyzed
        -d | --funcDel          method to use for filtering functional deleteriousness: CADD, metaSVM, radialSVM, MPC
        -p | --pldiff-cutoff    minimum required value of plDiff for a sample to be included in further processing
                 7:low stringency 8:high stringency
Options:
        -o | --outputDir        output directory to write files; default: current directory
        -x | --excel            FLAG: use to cause more header information for display in spreadsheet
        -s | --shift            FLAG: use to specify that nonframeshift indels should be filtered out
        -c | --coverage         required read depth for a variant to be accepted; default: 8
        -a | --africa-special   FLAG: filter specific info fields
```

Here is an example with commonly used options:

```bash
$ /home/user/programs/exome_burden_script.py \
 --frequency 5e-05 \
 --vcf my_annotated_samples.vcf \
 --funcDel CADD \
 --outputDir /home/user/results \
 --shift \
 --coverage 10
```

ALL PARAMETERS MUST BE SPECIFIED FOR THE PROGRAM TO RUN. Also, do NOT use the "--excel" option when following this work flow. It is intended as a last step in a different work flow not described here.

Regardless of the options you use, this program creates two files in the specified output directory.
1. my_annotated_samples_5e-05_cadd_counts_table.tsv
2. my_annotated_samples_5e-05_cadd_variant.table

The first one is a tab-delimited file with the Fisher exact p-values of all the various mutation classifications explained above. The second file is a "melted VCF". That is all data is turned "sideways" for further analysis.

The general formula for output file names is:
<input_vcf_file_name>_<frequency>_<functional_deleteriousness_method>_counts_table.tsv
This way, you can keep track of which output is from a specific sample and which parameter/option values you used.

You now need to run the programs on your controls WITH THE SAME OPTIONS USED BEFORE. Actually, you can use a different output directory, if you wish.

```bash
$ /home/user/programs/exome_burden_script.py \
 --frequency 5e-05 \
 --vcf my_annotated_controls.vcf \
 --funcDel CADD \
 --outputDir /home/user/results \
 --shift \
 --coverage 10
```

## Determine p-values: gene_burden_fisher_exact.py
This program requires the the two files just previously generated and the number of samples in each.

$ ./gene_burden_fisher_exact.py
usage: ./gene_burden_fisher_exact.py parameters
Parameters:
        -a | --cases-file       CASES gene burden output
        -n | --num-cases        number of cases in gene burden output
        -o | --controls-file    CONTROLS gene burden output
        -m | --num-controls     number of controls in gene burden output
Options:
        -i | --individual-counts         FLAG: independent counts yield independent p-values

```bash
$ /home/user/programs/gene_burden_fisher_exact.py \
 --cases-file my_annotated_samples_5e-05_cadd_counts_table.tsv \
 --num-cases 100
 --controls-file my_annotated_controls_5e-05_cadd_counts_table.tsv \
 --num-controls 2700 \
 > /home/user/results/fisher_exact_5e-05_results.tsv
```

These columns are listed once for the total variant count of a gene, then listed again for the number of samples containing significant variants for that gene. Because a single sample may contribute several variants, these numbers can be very different.

At this point, you need to choose a header or set of headers that you want to work with. Generally, my previous work focused on "Total_Missense", "Total_Deleterious", "Total_LoF" for samples (the second set of columns). This decision focuses on the number of samples contributing variants to a gene.

Once this decision is made, you need to use this column to filter out all genes that do not have a p-value above your significance threshold. The threshold we used was 2.5E-6. This was achieved by dividing 0.05 (a normal p-value threshold) by 20,000 (the approximate number of human genes). Thus, this value is corrected for multiple testing.
