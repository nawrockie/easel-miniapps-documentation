EPN, Thu Nov 10 2016

README 
Easel mini-applications (miniapps) for sequence file handling and
calculating statistics.

Organization of this file:
INTRODUCTION
GENERAL INSTRUCTIONS (APPLY TO ALL MINIAPPS)
USAGE EXAMPLES FOR SELECTED MINIAPPS
 - esl-seqstat  :: show simple statistics on a sequence file
 - esl-sfetch   :: retrieve sequence(s) from a file
 - esl-selectn  :: select n lines randomly from a file
 - esl-alistat  :: show summary statistics for a multiple sequence alignment file
 - esl-reformat :: convert between sequence file formats
 - esl-alimanip :: manipulate a multiple sequence alignment
 - esl-alipid   :: calculate pairwise %id for each seq pair in an MSA

============
INTRODUCTION
============

Easel is a sequence analysis C library that is used by the HMMER and
Infernal software packages. Easel contains 21 executable programs
called 'miniapps' that perform various tasks related to sequence file
and alignment file handling and manipulation. All of these
applications are pretty fast and use very efficient C code written
(mainly) by Sean Eddy.

The 21 miniapps are in /usr/local/infernal/1.1.1/bin/:

(This directory is accessible from the cbbdev machines and I think
also from the iebdev machines, but I'm not sure because I can't login
to those.)

All 21 miniapps are: 
1.  esl-afetch     :: retrieve multiple sequence alignment(s) from a file
2.  esl-alimanip   :: manipulate a multiple sequence alignment
3.  esl-alimanip   :: manipulate a multiple sequence alignment
4.  esl-alimap     :: map two alignments to each other
5.  esl-alimask    :: remove columns from a multiple sequence alignment
6.  esl-alimerge   :: merge alignments based on their reference (RF) annotation
7.  esl-alipid     :: calculate pairwise %id for each seq pair in an MSA
8.  esl-alistat    :: show summary statistics for a multiple sequence alignment file
9.  esl-cluster    :: clusters tabular data file
10. esl-compalign  :: compare two multiple sequence alignments
11. esl-compstruct :: calculate accuracy of RNA secondary structure predictions
12. esl-histplot   :: collate a data histogram, output xmgrace datafile
13. esl-mask       :: mask sequences in a sequence file
14. esl-reformat   :: convert between sequence file formats
15. esl-reformat   :: convert between sequence file formats
16. esl-selectn    :: select n lines randomly from a file
17. esl-seqrange   :: determine range of sequences for one of many parallel processes
18. esl-sfetch     :: retrieve sequence(s) from a file
19. esl-shuffle    :: shuffling or generating random sequences
20. esl-ssdraw     :: draw postscript secondary structure diagrams
21. esl-weight     :: calculate sequence weights for an alignment

Only some of these (the ones I find the most generally useful) are
described in further detail in the USAGE EXAMPLES FOR SELECTED MINIAPPS
section.

============================================
GENERAL INSTRUCTIONS (APPLY TO ALL MINIAPPS)
============================================
General instruction 1. To determine the usage and available command
line options for any of these miniapps, use the -h option, like this:

$ esl-seqstat -h
# esl-seqstat :: show simple statistics on a sequence file
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-seqstat    [options] <seqfile>

 where general options are:
  -h             : help; show brief info on version and usage
  -a             : report per-sequence info line, not just a summary
  -c             : count and report residue composition
  --informat <s> : specify that input file is in format <s>
  --rna          : specify that <seqfile> contains RNA sequence
  --dna          : specify that <seqfile> contains DNA sequence
  --amino        : specify that <seqfile> contains protein sequence

-------------------
General instruction 2. Most of these miniapps take as input a sequence
or alignment file, but you can also pipe in standard input by using
the '-' character in place of the sequence name. This allows
sophisticated users to string commands together, using the output of
one miniapp command as input to another one.  There are examples of
this below in the sections esl-seqstat and esl-sfetch.

====================================
USAGE EXAMPLES FOR SELECTED MINIAPPS
====================================
esl-seqstat :: show simple statistics on a sequence file

This script is useful for getting summary statistics on a file if run
in default mode:

$ esl-seqstat /net/snowman/vol/export2/mcveigh/arb_db/arb_128/arb-silva.de_2016-10-26_id380334_tax_silva.fasta
Format:              FASTA
Alphabet type:       RNA
Number of sequences: 582392
Total # residues:    824637340
Smallest:            330
Largest:             3469
Average length:      1415.9

--------------------------------
Some useful options are: 
  -a             : report per-sequence info line, not just a summary
  -c             : count and report residue composition


The -a option prints the name, length and description of each sequence
in the file:

$ esl-seqstat -a /net/snowman/vol/export2/mcveigh/arb_db/arb_128/arb-silva.de_2016-10-26_id380334_tax_silva.fasta
= A16379.1.1485            1485 Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;Pasteurellaceae;Haemophilus;[Haemophilus] ducreyi
= A93610.1.1301            1301 Archaea;Euryarchaeota;Thermococci;Thermococcales;Thermococcaceae;Thermococcus;unidentified
= AAAA02020710.72.1593     1522 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Pantoea;Oryza sativa Indica Group (long-grained rice)
= AAAA02029391.4444.5952     1509 Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Comamonadaceae;Acidovorax;Oryza sativa Indica Group (long-grained rice)
...output truncated...

Each sequence information line starts with the = character.
You can use this to get a list of all sequences sorted by length,
like this: 
$ esl-seqstat -a /net/snowman/vol/export2/mcveigh/arb_db/arb_128/arb-silva.de_2016-10-26_id380334_tax_silva.fasta | grep ^\= | sort -nk 3
= LFKQ01000006.1.330        330 Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus hominis subsp. hominis
= LNMG01000063.1.331        331 Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium
= LNAL01000006.1292456.1292788      333 Bacteria;Bacteroidetes;Cytophagia;Cytophagales;Cytophagaceae;Hymenobacter;Cytophagales bacterium MC1A
= CP011799.6213695.6214029      335 Bacteria;Actinobacteria;Actinobacteria;Streptomycetales;Streptomycetaceae;Streptomyces;Streptomyces sp. PBH53
...output truncated...

$ esl-seqstat -c /net/snowman/vol/export2/mcveigh/arb_db/arb_128/arb-silva.de_2016-10-26_id380334_tax_silva.fasta 
Format:              FASTA
Alphabet type:       RNA
Number of sequences: 582392
Total # residues:    824637340
Smallest:            330
Largest:             3469
Average length:      1415.9

Residue composition:
residue: A    208091795  0.2523
residue: C    189603341  0.2299
residue: G    259421501  0.3146
residue: U    167280509  0.2029
residue: R        17412  0.0000
residue: Y        15062  0.0000
residue: M         6997  0.0000
residue: K         9231  0.0000
residue: S        10799  0.0000
residue: W         7156  0.0000
residue: H          285  0.0000
residue: B          354  0.0000
residue: V          318  0.0000
residue: D          358  0.0000
residue: N       172222  0.0002

======================================================
esl-sfetch   :: retrieve sequence(s) from a file

This miniapp is helpful for fetching sequences or subsequences from a
file. It can be used in 4 modes, before any of these can be used, the
sequence file you want to fetch from must be 'indexed' like this:

$ esl-sfetch arb-silva.de_2016-10-26_id380334_tax_silva.fasta 

This will create a file called
arb-silva.de_2016-10-26_id380334_tax_silva.fasta.ssi
which esl-sfetch will use to subsequently fetch sequences quickly.
 
--------------------
Four common modes of usage:

A. fetch a single sequence from the file:
esl-sfetch [options] <sqfile> <name>        (one seq named <name>)

example:
$ esl-sfetch arb-silva.de_2016-10-26_id380334_tax_silva.fasta LFKQ01000006.1.330
>LFKQ01000006.1.330 Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus hominis subsp. hominis
AGAGUUUGAUCCUGGCUCAGGAUGAACGCUGGCGGCGUGCCUAAUACAUGCAAGUCGAGCGAACAGACGAGGAGCUUGCU
CCUUCGACGUUAGCGGCGGACGGGUGAGUAACACGUAGGUAACCUACCUAUAAGACUGGGAUAACUUCGGGAAACCGGAG
CUAAUACCGGAUAAUAUUUCGAACCGCAUGGUUCGAUAGUGAAAGAUGGCUUUGCUAUCACUUAUAGAUGGACCUGCGCC
GUAUUAGCUAGUUGGUAAGGUAACGGCUUACCAAGGCAACGAUACGUAGCCGACCUGAGAGGGUGAUCGGCCACACUGGA
ACUGAGACAC

--------------------
B. fetch a single subsequence from the file, using the -c option
combination:
Usage: esl-sfetch [options] -c <s> <sqfile> <name>

  On command line, subseq coords are separated by any nonnumeric, nonspace character(s).
  for example, -c 23..100 or -c 23/100 or -c 23-100 all work.

  Additionally, to retrieve a suffix to the end, omit the end coord; -c 23: will work.
  By default, the subseq will be named <source name>/<from>-<to>. To assign a name of
  your choice, use -n <newname>.

example:
$ esl-sfetch -c 30..100 arb-silva.de_2016-10-26_id380334_tax_silva.fasta LFKQ01000006.1.330 
>LFKQ01000006.1.330/30-100 Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus hominis subsp. hominis
UGGCGGCGUGCCUAAUACAUGCAAGUCGAGCGAACAGACGAGGAGCUUGCUCCUUCGACG
UUAGCGGCGGA

--------------------
C. fetch a list of sequences from the file, using the -f option. 
Usage: esl-sfetch [options] -f <sqfile> <namefile> (all seqs in <namefile>)

example:
$ cat 4.list 
LFKQ01000006.1.330
LNMG01000063.1.331
LNAL01000006.1292456.1292788
CP011799.6213695.6214029

$ esl-sfetch -f arb-silva.de_2016-10-26_id380334_tax_silva.fasta 4.list
>LFKQ01000006.1.330 Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus hominis subsp. hominis
AGAGUUUGAUCCUGGCUCAGGAUGAACGCUGGCGGCGUGCCUAAUACAUGCAAGUCGAGCGAACAGACGAGGAGCUUGCU
CCUUCGACGUUAGCGGCGGACGGGUGAGUAACACGUAGGUAACCUACCUAUAAGACUGGGAUAACUUCGGGAAACCGGAG
CUAAUACCGGAUAAUAUUUCGAACCGCAUGGUUCGAUAGUGAAAGAUGGCUUUGCUAUCACUUAUAGAUGGACCUGCGCC
GUAUUAGCUAGUUGGUAAGGUAACGGCUUACCAAGGCAACGAUACGUAGCCGACCUGAGAGGGUGAUCGGCCACACUGGA
ACUGAGACAC
>LNMG01000063.1.331 Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium
UUUUUAUGAGAGUUUGAUCCUGGCUCAGGACGAACGCUGGCGGCGUGCCUAAUACAUGCAAGUCGAACGCUUCUUUUUCC
ACCGGAGCUUGCUCCACCGGAAAAAGAAGAGUGGCGAACGGGUGAGUAACACGUGGGUAACCUGCCCAUCAGAAGGGGAU
AACACUUGGAAACAGGUGCUAAUACCGUAUAACAAUCAAAACCGCAUGGUUUUGAUUUGAAAGGCGCUUUCGGGUGUCGC
...output truncated...

------------------
D. fetch a list of subsequences from the file, and specify the names
of the new sequences, using the combination of the -C and -f options:
Usage: esl-sfetch [options] -Cf <sqfile> <namefile with subseq coords> (all subseqs in <namefile>)

  In retrieving subsequences listed in a file (-C -f, or just -Cf), each line of the file
  is in GDF format: <newname> <from> <to> <source seqname>, space/tab delimited.

example:
$ cat 4subseq.list
LFKQ01000006.1.330-subseq           30 100 LFKQ01000006.1.330
LNMG01000063.1.331-subseq           40 110 LNMG01000063.1.331
LNAL01000006.1292456.1292788/50-120 50 120 LNAL01000006.1292456.1292788
CP011799.6213695.6214029/60-130     60 130 CP011799.6213695.6214029

$ esl-sfetch -Cf arb-silva.de_2016-10-26_id380334_tax_silva.fasta 4subseq.list
>LFKQ01000006.1.330-subseq Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus hominis subsp. hominis
UGGCGGCGUGCCUAAUACAUGCAAGUCGAGCGAACAGACGAGGAGCUUGCUCCUUCGACG
UUAGCGGCGGA
>LNMG01000063.1.331-subseq Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium
GCGGCGUGCCUAAUACAUGCAAGUCGAACGCUUCUUUUUCCACCGGAGCUUGCUCCACCG
GAAAAAGAAGA
>LNAL01000006.1292456.1292788/50-120 Bacteria;Bacteroidetes;Cytophagia;Cytophagales;Cytophagaceae;Hymenobacter;Cytophagales bacterium MC1A
AAUACAUGCAAGUCGAACGGGCGCAGCAAUGCGUCAGUGGCGCACGGGUGCGUAACGCGU
AGGCAAUCUGC
>CP011799.6213695.6214029/60-130 Bacteria;Actinobacteria;Actinobacteria;Streptomycetales;Streptomycetaceae;Streptomyces;Streptomyces sp. PBH53
CAAGUCGAACGAUGAACCUCCUUCGGGAGGGGAUUAGUGGCGAACGGGUGAGUAACACGU
GGGCAAUCUGC

Here is the complete usage:

$ esl-sfetch -h
# esl-sfetch :: retrieve sequence(s) from a file
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-sfetch [options] <sqfile> <name>        (one seq named <name>)
Usage: esl-sfetch [options] -f <sqfile> <namefile> (all seqs in <namefile>)
Usage: esl-sfetch [options] --index <sqfile>       (index <sqfile>)

 where general options are:
  -h     : help; show brief info on version and usage
  -o <f> : output sequences to file <f> instead of stdout
  -O     : output sequence to file named <key>
  -n <s> : rename the sequence <s>
  -r     : reverse complement the seq(s)

 Options for retrieving subsequences:
  -c <s> : retrieve subsequence coords <from>..<to>
  -C     : <namefile> in <f> contains subseq coords too

  On command line, subseq coords are separated by any nonnumeric, nonspace character(s).
  for example, -c 23..100 or -c 23/100 or -c 23-100 all work.

  Additionally, to retrieve a suffix to the end, omit the end coord; -c 23: will work.
  By default, the subseq will be named <source name>/<from>-<to>. To assign a name of
  your choice, use -n <newname>.

  In retrieving subsequences listed in a file (-C -f, or just -Cf), each line of the file
  is in GDF format: <newname> <from> <to> <source seqname>, space/tab delimited.

  When <start> coordinate is greater than <end>, for DNA or RNA, the reverse complement is
  retrieved; in protein sequence, this is an error. The -r option is another way to revcomp.

 other options:
  --informat <s> : specify that input file is in format <s>

=============================================================
esl-selectn  :: select n lines randomly from a file

This is useful if you want a random sample of anything, where the
items you want to choose from can be arranged as one-per-line in a
text file.

The random number generator algorithm used is the Mersenne Twister
algorithm.

The usage: 
$ esl-selectn -h 
# esl-selectn :: select n lines randomly from a file
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-selectn [-options] <n> <file>

 where general options are:
  -h         : show brief help on version and usage
  --seed <n> : set random number generator's seed to <n>  [0]


An example: to select 2 random lines from a the file 4.list:
$ esl-selectn 2 4.list 
LFKQ01000006.1.330
CP011799.6213695.6214029

Another example: say you want to choose 10 sequences randomly from a
file, you can combine the esl-seqstat -a option with esl-selectn like
this: 

$ esl-seqstat -a arb-silva.de_2016-10-26_id380334_tax_silva.fasta | grep ^\= | awk '{ print $2 }' | esl-selectn 10 -
KP204837.1.1423
HQ804831.1.1450
AB806097.1.1427
KF668479.1.1491
JN986325.1.1452
FM992729.1.1328
CT573572.2.1347
HQ783929.1.1421
JF150632.1.1352
JX240747.1.1486

Breakdown of this command:
'esl-seqstat -a arb-silva.de_2016-10-26_id380334_tax_silva.fasta' :
               calls esl-seqstat on the sequence file 

'| grep ^\=' : takes output of esl-seqstat and pipes it into grep to
               select only lines that begin with =

'| awk '{ print $2 }' : takes grep output and pipes into awk and
                        prints only 2nd token of each line (the
                        sequence name)

'| esl-selectn 10 -' : takes awk output and pipes into esl-selectn,
                       and selects 10 lines at random. The '-'
                       informs esl-seqstat that the input is not a
                       file but standard input

============================================================
esl-alistat  :: show summary statistics for a multiple sequence alignment file

Usage:
# esl-selectn :: select n lines randomly from a file
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-selectn [-options] <n> <file>

 where general options are:
  -h         : show brief help on version and usage
  --seed <n> : set random number generator's seed to <n>  [0]
<[(16_1108_esl_readme_rich_mcveigh)]> esl-alistat -h
# esl-alistat :: show summary statistics for a multiple sequence alignment file
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-alistat [options] <msafile>

 where options are:
  -h             : help; show brief info on version and usage
  -1             : use tabular output, one line per alignment
  --informat <s> : specify that input file is in format <s>
  --amino        : <msafile> contains protein alignments
  --dna          : <msafile> contains DNA alignments
  --rna          : <msafile> contains RNA alignments

 small memory mode, requires --amino,--dna, or --rna and --informat pfam:
  --small : use minimal RAM (RAM usage will be independent of aln size)

 optional output files:
  --list <f>   : output list of sequence names in alignment(s) to file <f>
  --icinfo <f> : print info on information content alignment column
  --rinfo <f>  : print info on # of non-gap residues in each column to <f>
  --pcinfo <f> : print per-column   posterior probability info to <f>
  --psinfo <f> : print per-sequence posterior probability info to <f>
  --iinfo <f>  : print info on # of insertions b/t all non-gap RF cols to <f>
  --cinfo <f>  : print per-column residue counts to <f>
  --noambig    : with --cinfo, do not count ambiguous residues
  --bpinfo <f> : print per-column base-pair counts to <f>
  --weight     : with --*info files, weight counts using WT annotation from msa


Example:
(The file RF00005.stk is the tRNA 'seed' alignment from the RNA
database in stockholm format (https://en.wikipedia.org/wiki/Stockholm_format))

$ esl-alistat RF00005.stk
Alignment number:    1
Alignment name:      tRNA
Format:              Stockholm
Number of sequences: 954
Alignment length:    118
Total # residues:    70041
Smallest:            62
Largest:             93
Average length:      73.4
Average identity:    44%
//

Some other potentially useful options are:
  --list <f>   : output list of sequence names in alignment(s) to file <f>
  --icinfo <f> : print info on information content alignment column
  --rinfo <f>  : print info on # of non-gap residues in each column to <f>
  --pcinfo <f> : print per-column   posterior probability info to <f>
  --psinfo <f> : print per-sequence posterior probability info to <f>
  --iinfo <f>  : print info on # of insertions b/t all non-gap RF cols to <f>
  --cinfo <f>  : print per-column residue counts to <f>
  --bpinfo <f> : print per-column base-pair counts to <f>

I don't include any examples of those examples here.

=============================================================
esl-reformat :: convert between sequence file formats

This script is only really useful for alignment formats. It allows you
to convert between the following formats: a2m, afa (aligned FASTA),
clustal, clustal-like, pfam, phylip, psiblast, selex and stockholm.

<[(16_1108_esl_readme_rich_mcveigh)]> esl-reformat -h
# esl-reformat :: convert between sequence file formats
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-reformat [-options] <format> <seqfile>
  Output format choices: Unaligned      Aligned    
                         -----------    -------    
                         fasta          a2m        
                                        afa        
                                        clustal    
                                        clustallike
                                        pfam       
                                        phylip     
                                        phylips    
                                        psiblast   
                                        selex      
                                        stockholm  


  where options are:

  -d             : convert to DNA alphabet (U->T)
  -h             : help; print brief info on version and usage
  -l             : convert to lower case
  -n             : remove DNA IUPAC codes; convert ambig chars to N
  -o <s>         : send output to file <f>, not stdout
  -r             : convert to RNA alphabet (T->U)
  -u             : convert to upper case
  -x             : convert non-IUPAC chars (e.g. X) in DNA to N
  --gapsym <s>   : convert all gaps to character <c>
  --informat <s> : input sequence file is in format <s>
  --mingap       : remove columns containing all gaps (seqfile=MSA)
  --keeprf       : with --mingap, keep all nongap #=GC RF columns
  --nogap        : remove columns containing any gaps (seqfile=MSA)
  --wussify      : convert old RNA structure markup lines to WUSS
  --dewuss       : convert WUSS RNA structure markup to old format
  --fullwuss     : convert simple WUSS notation to full (output) WUSS
  --ignore <s>   : ignore input seq characters listed in string <s>
  --acceptx <s>  : accept input seq chars in string <s> as X
  --rename <s>   : rename and number each sequence <s>.<n>
  --replace <s>  : <s> = <s1>:<s2> replace characters in <s1> with those in <s2>
  --small        : use minimal RAM, input must be pfam, output must be afa or pfam

An example:

$ esl-reformat afa ssu.sto > ssu.afa

==============================================
esl-alimanip :: manipulate a multiple sequence alignment

Usage:
$ esl-alimanip -h
# esl-alimanip :: manipulate a multiple sequence alignment
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-alimanip [options] <msafile>

where basic options are:
  -h              : help; show brief info on version and usage
  -o <f>          : output the alignment to file <f>, not stdout
  --informat <s>  : specify that input file is in format <s>
  --outformat <s> : specify that output format be <s>
  --devhelp       : show list of undocumented developer options

options for removing/reordering/trimming sequences:
  --lnfract <x> : remove sequences w/length < <x> fraction of median length
  --lxfract <x> : remove sequences w/length > <x> fraction of median length
  --lmin <n>    : remove sequences w/length < <n> residues
  --lmax <n>    : remove sequences w/length > <n> residues
  --detrunc <n> : remove seqs w/gaps in >= <n> 5' or 3'-most non-gap #=GC RF cols
  --xambig <n>  : remove sequences with >= <n> ambiguous residues
  --seq-r <f>   : remove sequences with names listed in file <f>
  --seq-k <f>   : remove all seqs *except* those listed in <f>
  --small       : w/--seq-r or --seq-k use minimal RAM (no seq reordering)
  --k-reorder   : with --seq-k <f>, reorder sequences to order in <f>
  --seq-ins <n> : keep only seqs w/an insert after non-gap RF col <n>
  --seq-ni <n>  : w/--seq-ins require at least <n> residue insertions  [1]
  --seq-xi <n>  : w/--seq-ins require at most  <n> residue insertions  [1000000]
  --trim <f>    : trim aligned seqs in <msafile> to subseqs in <f>
  --t-keeprf    : w/--trim keep GC RF annotation in msa, if it exists
  --minpp <x>   : replace residues with posterior probabilities < <x> with gaps
  --tree <f>    : reorder MSA to tree order following SLC, save Newick tree to <f>
  --reorder <f> : reorder seqs to the order listed in <f>, all seqs must be listed

options for adding/removing alignment annotation:
  --mask2rf <f> : set #=GC RF as x=1, gap=0 from 1/0s in 1-line <f>
  --m-keeprf    : with --mask2rf, do not overwrite nongap RF characters with 'x'
  --num-all     : add annotation numbering all columns
  --num-rf      : add annotation numbering the non-gap RF columns
  --rm-gc <s>   : remove GC <s> markup, <s> must be RF|SS_cons|SA_cons|PP_cons
  --sindi       : annotate individual secondary structures by imposing consensus
  --post2pp     : convert infernal 0.72-1.0.2 POST posterior prob annotation to PP

options for specifying bio alphabet:
  --amino : <msafile> contains protein alignments
  --dna   : <msafile> contains DNA alignments
  --rna   : <msafile> contains RNA alignments

One example:

- Remove a set of sequences from an alignment:
  First get a list of sequences in an alignment using esl-alistat
  --list:
  $ esl-alistat --list RF00005.list RF00005.stk

  Then manipulate the sequences in RF00005.list to only include the
  sequences that you want to remove, and rename it RF00005.remove.list:

  $ esl-alimanip --seq-r RF00005.remove.list RF00005.stk > RF00005.sub.stk

  Alternatively, the --seq-k option will keep only those sequences
  listed in the list file and remove all others.

============================================
esl-alipid   :: calculate pairwise %id for each seq pair in an MSA

Usage:
$ esl-alipid -h 
# esl-alipid :: calculate pairwise %id for each seq pair in an MSA
# Easel i1.1.1 (July 2014)
# Copyright (C) 2014 HHMI Janelia Farm Research Campus
# Freely distributed under the Janelia Software License.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: esl-alipid [options] <msafile>

Options:
  -h              : help; show brief info on version and usage
  --informat <s>  : specify the input MSA file is in format <s>
  --outformat <s> : write the output MSA in format <s>  [Clustal]
  --dna           : use DNA alphabet
  --rna           : use RNA alphabet
  --amino         : use protein alphabet

Example:

$ esl-alipid RF00005.stk
AB003409.1/96-167          AB009835.1/1-71             28.17     20     71
AB003409.1/96-167          AB013372.1/8-81             56.94     41     72
AB003409.1/96-167          AB013373.1/3754-3825        62.50     45     72
...output truncated...

The 4 output columns are: 
column 1: sequence 1
column 2: sequence 2
column 3: percent identity between sequence 1 and sequence 2 given the alignment
column 4: number of identical aligned residues between sequence 1 and 2
column 5: minimum length of the two sequences is 

The percent identity is calculated as the number of identities divided
by the minimum length of the sequences multiplied by 100 ((column
4/column5) * 100).

One line is printed for every unique pair of sequences in the
alignment.








