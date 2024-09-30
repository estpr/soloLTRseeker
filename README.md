# soloLTRseeker: a new tool to identify soloLTR sequences in eukaryotic genomes 

SoloLTRs are Transposable Element (TE) derived sequences generated during unequal homologous recombination events in which, the internal region of an LTR retrotranposon (LTR-RT) and one of its two Long Terminal Repeats (LTRs) are deleted, producing a single LTR structure. We have developed soloLTRseeker, a computational pipeline, to systematically quantify and compare the intensity and impact of these structures. Although still **UNDER DEVELOPMENT**, we carried out a number of tests in several angiosperm and gymnosperm species and, so far, the results show that it annotates consistent soloLTRs. Benchmarking, error estimation and further documentation on how to run it and juice every single bit of info coming soon.

## PIPELINE OVERVIEW

**DISCLAIMER:** Essentially, this is all about parsing other tool’s outputs (cd-hit, TEsorter, BLAST, water, you name it!) using awk and some additional basic commands. Be that as it may, order and care matter when it comes to get relevant results. So, there we go:

The pipeline comprises two main modules, viz.:

1.	LTR library construction

Input LTR-RTs are filtered to keep only those entries fulfilling a set of user defined parameters. This initial effort is undertaken to generate a curated set of LTR sequences, which we find is paramount to optimise the analysis. Partial or low-quality queries, for instance, are likely to cause spurious hits surely pervading the final tally.  

Two aspects linked to the integrity of the LTR sequences are considered at this stage. The user might custom a lower limit for the **length** of any LTR incorporated in the LTRlibrary and a **maximum size divergence** between counterparts, switching on the appropriate flags on the command line (see USAGE). In this sense, patterns of sequence conservation are also assessed evaluating the presence of **TG..CA** at both 5 and 3 prime ends of all LTRs. By default, 3 out of 4 nucleotides must match that canonical pattern. LTR library is finally condensed to a set of **non-redundant representatives** to reduce the runtime while preserving the genetic diversity. Different clustering parameters, like percentage identity and coverage cutoffs might be adjusted at this point (i.e., p, aL and aS arguments that are all specified for the cd-hit run).

Additionally, TEsorter might be run within the pipeline too. If that is the case, lineage-based classification of full-length LTR-RTs is used to group soloLTRs candidates. This option might come in handy to further study TE dynamics private to each category.

2.	soloLTR mining 

**BLAST** is used to identify candidate loci across the target genome. Initial hits are collapsed and classify as per the query which showed the highest percentage identity. By default, hits must cover no less than 0.99 of the sequence that produced them, though more permissive thresholds can be set. Next, **signs of truncated TEs** are searched out by comparing the flanking sequences of each BLAST hit to the internal domain edges of all intact elements. Our tests indicate that the resolution of this homology search is maximized surveying 200 bp length sequences. Yet, this value might be adapted, if the user finds it convenient. Finally, upstream and downstream sequences are locally aligned to locate both Target Sites Duplication (**TSDs**). In general, TSD length ranges narrowly from 5 to 6 bp. Consequently, this step is configured to align pairs of 6-mers allowing 1 mismatch between them.


## INPUT FILES

genome file (fasta)  

full length and LTRs coordinates of intact LTR retrotransposons (gff3)  

OPTIONAL: to run a partial analysis, one chromosome name, as it appears in the fasta file header, can be provided ("string")  


## OUTPUT FILES 

`sample_soloLTR.gff3`. This file includes all hits that fulfilled the conditions specified for the run. Its format follows gff3 conventions: *seq.id source type start end score strand phase attributes*. Tags includes in the last field appear separated by semi-colon:

*ID*: a unique id is formed combining sample name, chromosome and start and end coordinates.  
*Name*: all soloLTRs are numbered sequentially to create a fixed id for each one.  
*Classification*: each soloLTR inherits the feature’s type term (as appears in the 3rd column of the input gff3 file) from the LTR query which produced the BLAST hit with the greatest identity.  
*Sequence_ontology*: SO:0001003 (http://www.sequenceontology.org).  
*Lineage*: If TEsorter flag is activated, then each soloLTR is assigned the same taxonomic classification as the LTR query which produced the BLAST hit with highest identity. Otherwise, a NA label will appear under this tag.   
*motif*: This is the di-nucleotide at the start and at the end of the soloLTR sequence.  
*TSD*: TSDs sequences are reported together with the number of mismatches found in the alignment.  

`sample_soloLTR.fa`: fasta file including the sequence of all annotated soloLTRs. ID tags reported in the gff3 file is used as fasta headers: *>soloLTR/sample/chr/start/end*

`sample_LTRRT_lengths.txt`: this file reports a number of relevant features of all LTR-RTs considered during the LTR library module: *LTR_query_id  fl_LTRRT_length lLTR_length  rLTR_length  idom_length  motif  b_n_filter  tg_ca_filter  cd_hit_step*. Particularly, the last three fields of the table show the outcome of the three-step filtering process carried out in that part of the pipeline (i.e., min length and max difference between pairs of LTRs, dinucleotide pattern assessment, and cd-hit clustering). Besides, LTR_query_id contain coordinates of the intact LTR-RT, feature’s type term, and its lineage classification, if applicable. Finally, apart from an itemized length report, motif field covers the di-nucleotides sites at the start and end of the lLTR sequence.

`sample_LTR_BLAST_overlap.txt`: each BLAST-hit initially found in the mining module is reported in this file as follows: *chr  start  end  queryID  strand  piden  locus  locus_type  BLAST_single_hit  homology  TSD*. Here, queryID is formatted as in sample_LTRRT_lengths.txt file. It is worth noting that locus_type is used to further describe how all hits overlap in a specific locus. In particular, complete refers to those loci whose all hits are enclosed by the longest one. Conversely, partial tag accounts for those loci whose ends bear single overhanging nucleotides in reference to the entire overlapping region. Again, the last three columns allow tracing back the fate of every hit following a binary key (i.e., PASS/FAIL).  

`sample_hit_count.txt`: A summary file gathering the overall counts at each step of the pipeline.  

`out.log`: this file is intended to report run settings and track its progress. Warnings and unexpected behaviours that don’t make the run crash, should be appended in this file, so it’s worth keeping an eye on it.


## USAGE

soloLTRseeker -h

    USAGE
      soloLTRseeker [-b] [-n] [-t] [-c] [-l] [-s] [-o] [-q] [-p] [-m] [-r] [-u] [-h] ann_file fasta_file [chromosome_name]

    DESCRIPTION
      soloLTRseeker v0.0

    INPUT FILES
      ann_file: annotation in gff3 format
      fasta_file

    OPTIONAL ARGUMENTS
      -b max length difference between LTR pairs;default = 50
      -n min LTR length;default = 1
      -t LTR-RT lineage classification via TEsorter
      -c percentage identity (cd-hit);default = 0.95
      -l longest seq coverage (cd-hit);default 0.9
      -s shortest seq coverage (cd-hit);default 0.3
      -o query coverage (BLAST);default = 0.99
      -q query length for homology search;default = 200
      -p query length for TSD annotation;default = 6
      -m mismatch allowed in TSD seq;default = 1
      -r split BLAST LTR query library into # files;default = 10
      -u analyse BLAST query files in parallel;default = T
      -h print USAGE, DESCRIPTION and OPTIONAL ARGUMENTS

      [chromosome_name] for partial analysis


## QUICK START

Provided your system has installed all necessary dependencies (see them below), it just takes executing the main script to get the pipeline started (be sure all the sh files are placed together in the same directory). You can either call it directly or after adding it to your $PATH, whatever fits better your habitual practices.

`soloLTRseeker   /abs/path/to/sample.gff3  /abs/path/to/sample.fasta`

Two arguments deserve further clarification,though. Particularly, -r ARG will adopt different meanings depending on -u. If the latter is set to F, then the former, -r, refers to the number of threads allocated for a single BLAST process. Conversely, when -u ARG is T, the LTR library will be split in the number of files specified by -r and they will be analysed in parallel as independent BLAST processes. 
Even though, -r and -u were both conceived with runtime in mind, we must admit now that this *parallelisation* option falls entirely within the realm of ad hoc implementations, as there are multiple factors influencing this variable.  
However, based on our little experience we might suggest that, if server’s capacity is not a real concern, then the faster option for analysing a small/medium genome might be the second one. On the contrary, if memory resources are restricted, it might be better going for the first alternative.

In any case, we are aware that shell scripting is a rather limited base for an ambitious endeavour. That being so, the pipeline does modestly seek to channel the result of much greater tools, reporting along the way as much information as possible.

## DEPENDENCIES

BEDtools v2.27.1 - Quinlan AR and Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 26, 6, pp. 841–842 (2010)  

BLAST 2.14.1+ - Altschul, S. F. et al. Basic local alignment search tool. *Journal of Molecular Biology* 215, 403-410 (1990)  

CD-HIT version 4.8.1 - Fu, Limin et al. CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* (Oxford, England) vol. 28,23 (2012)  

EMBOSS 6.6.0.0 - Rice, P., Longden, I. and Bleasby, A. EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16, 276–277 (2000)  

TEsorter 1.4.6 - Zhang, R.-G. et al. TEsorter: An accurate and fast method to classify LTR-retrotransposons in plant genomes. *Horticulture Research* 9, uhac017 (2022)  


## A WORD OF WARNING

To ensure its performance, we plan to upload a toy run, with light files and a dry run library, to let you test the pipeline first on your device. Please, do contact us if you need to run soloLTRseeker imminently.


