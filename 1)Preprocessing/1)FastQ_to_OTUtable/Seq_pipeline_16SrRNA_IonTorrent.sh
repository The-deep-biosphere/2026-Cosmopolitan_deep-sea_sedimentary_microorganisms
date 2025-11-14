# Visual Studio Code
# Comment: Ctrl + K + C
# Uncomment: Ctrl + K + U

#===================================================================================================================================================================================================
# REQUIREMENTS AND INFO ON USING THE PIPELINE
#===================================================================================================================================================================================================

### For non-original data:
### (skip to pipeline if you are using your own sequencing data)
### download data from NCBI/ENA archives



# fasterq-dump SRRXXXXXXXXXX SRRXXXXXXXX .... or
# for (( i = 579; i <= 622; i++ ));   do   fasterq-dump SRR17289$i; done
# fasterq-dump SRRXXXXXXXXXX SRRXXXXXXXX .... or
# for (( i = 579; i <= 622; i++ ));   do   fasterq-dump SRR17289$i; done
### or
# cat SRR_acc_list.txt | xargs fasterq-dump

# unzip fastq.gz files
# gunzip *.gz	

### Merge forward and reverse reads
#  pandaseq -f  *forward_1.fastq* -r *reverse_2.fastq* -p *forward primer sequence* -q *reverse primer sequence*  -w  *output.fastq* -F -t 0.6

### ** if you have both forward and reverse reads and want to merge. There are other mergers you can use or you can use just the forward reads to run through the pipeline and compare with Ion Torrent reads. 
### - See VSEARCH paper for merging pair-ends option in VSEARCH instead

### FastQC
### - open fastq files in FastQC to check quality and at what length quality deteriorates to determine the length to trim. (if using downloaded data check for adapter content and add to cutadapt if additional adapters need removed)


### Pipeline
### VSEARCH pipeline reference
### https://github.com/torognes/vsearch/wiki/Alternative-VSEARCH-pipeline
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5075697/

### Setup 
## - create directory to work in, name this folder after your data

## - move fastq files, **silvamod138pr2.fasta**, **filter_seq_by_OTUtable.py**, and **LULU.R** to the working directory
## - Silvamod138 can be downloaded here: http://services.cbu.uib.no/supplementary/crest

## - install the following: fastqc, cutadapt, vsearch, crest4 (https://github.com/xapple/crest4), ncbi-blast+, and Rstudio. (Update your Python, by default installed)

### Remove primers 
### (If using merged sequences, remove "--discard-untrimmed")

### Change forward primer to the primer sequence used. If using pair-end data (forward and reverse reads) you can add funtion "-a " with the reverse primer sequence


#===================================================================================================================================================================================================
# THE PIPELINE
#===================================================================================================================================================================================================


gunzip *.gz	

!/bin/bash

VSEARCH=$(which vsearch)
THREADS=6

# Amount reads in fastq file.  (take notes of how much sequences are cut out over the pipeline)
for f in I*.fastq; do  
 awk '{s++}END{print s/4}' $f  
done

### Trim forward primers (~15 bp)
# Only forward primers, not all sequences have been measured until the reverse primer.
for f in *.fastq; do
    s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'.' '{print "2-" $1}' <<< $s)
     echo $s

     # Primers used.
     # Forward primer: 519f (5-́CAGCMGCCGCGGTAA-3)
     # Referse primer: 805r (5-́GACTACHVGGGTATCTAATCC-3)


 cutadapt -j 3 -g CAGCMGCCGCGGTAA --no-indels --discard-untrimmed -o $s.noprimers.fastq $f
 
 done

#Print number of reads:
# Amount reads in fastq file.  
for f in 2-*.fastq; do  
 awk '{s++}END{print s/4}' $f  
done

### Trim the sequences at 220bp.  
# This can be shortened or lengthened based on what you saw in FastQC. 220 is a good length usually. 

 for f in 2-*.fastq; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "3-" $2}' <<< $s)

     $VSEARCH --threads $THREADS \
         --fastq_filter $f \
         --fastq_maxns 0 \
         --fastq_trunclen 220 \
         --fastqout $s.trimmed.fastq
 done

# Amount reads in fastq file.  
for f in 3-*.fastq; do  
 awk '{s++}END{print s/4}' $f  
done


### Quality filtering at maxee = 2.
 #This sets the maximum errors allowed. It correlates the the Qulaity/Phred score given to each base - this score is the probability the base is incorrect. 
 #https://www.drive5.com/usearch/manual/exp_errs.html 
#default is 2 for IonTorrent (1 for Illlumina), but this number can be lowered for more aggresive filtering 

  for f in 3-*.fastq; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "4-" $2}' <<< $s)
     l=$(awk -F'-' '{print $2}' <<< $s)
     l=$(awk -F'.' '{print $1}' <<< $l)
     $VSEARCH --threads $THREADS \
         --fastq_filter $f \
         --relabel $l. \
         --fastq_maxee 2 \
         --fastaout $s.fasta \
         --fasta_width 0
 done

for f in 4-*.fasta; do  
 awk '{s++}END{print s/2}' $f  
done
### Dereplicate at sample level and relabel with sample_n
#This removes identical sequences within a sample

 for f in 4-*.fasta; do
     s=$(awk 'BEGIN{FS=OFS="."}{NF--; print}' <<< $f)
     echo
     echo ====================================
     echo Processing sample $s
     echo ====================================    
     s=$(awk -F'-' '{print "5-" $2}' <<< $s)
     l=$(awk -F'-' '{print $2}' <<< $s)
     l=$(awk -F'.' '{print $1}' <<< $l)
     $VSEARCH --threads $THREADS \
         --derep_fulllength $f \
         --strand plus \
         --sizeout \
         --relabel $l. \
         --fasta_width 0 \
         --output $s.derep.dfa
 done 

for f in 5-*.dfa; do  
 awk '{s++}END{print s/2}' $f  
done


### Cat all fasta files
# combine all files into one fasta file

cat 5-*.dfa > all.fasta

### Remove unneeded files. (not nessecary)

#rm 1-*.fastq 2-*.fastq 3-*.fastq 4-*.fasta 5-*.dfab 


### Dereplication of the fasta file.
# removes identical sequences 

#**derep.uc** is a way of mapping which sequences goes to which sample and how many times. In other words, it maps how much sequence replicates are present per sample.
 
 $VSEARCH --derep_fulllength all.fasta \
     --threads $THREADS \
     --minuniquesize 2 \
     --sizein \
     --sizeout \
     --fasta_width 0 \
    --uc derep.uc \
    --output derep.fasta

### Clustering at 97% similarity.
# De Novo clustering of  sequences with 97% similarity as default (similarity level can be adjusted).  See Rognes et. al (2016) for more in depth explanation
# Centroids representative sequence of OTU (OTU = a pool of multiple sequences that have in this case 97% similarity)
# How it works: for each sample the most abundant sequence is determined and the other sequences are compared to the most abundant one.

 $VSEARCH --cluster_size derep.fasta \
     --threads $THREADS \
     --id 0.97 \
     --strand plus \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --centroids centroids.fasta


### Sort centroids and remove singletons.
#This sorts the sequences by abundance and removes the sequences found only once

 $VSEARCH --sortbysize centroids.fasta \
     --threads $THREADS \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --minsize 2 \
     --output sorted.fasta


### denovo chimera detection.
# Chimera: problem during PCR; after a cycle the amplification of a sequence can be terminated in the middle causing wrong copies in the next cycle. Chimera detection checks if two sequences have been merged.
# Looks for and removes chimeras formed from clustering

#chimeras are produced in elongation of PCR when only part of sequenced in amplified and connects to other short sequences forming sequences that are part of one sequences and part of another.

$VSEARCH --uchime_denovo sorted.fasta \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --nonchimeras denovo.nonchimeras.fasta \


### Reference chimera detection against SILVA138.
# Blasts sequences against SILVA138 database. Removes chimeras
# Here,  silvamod138pr2.fasta  is used, normally  silvamod138.fasta is used. This file consists of some eukaryotes. 

 $VSEARCH --uchime_ref denovo.nonchimeras.fasta \
     --threads $THREADS \
     --db silvamod138pr2.fasta \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --dbmask none \
     --nonchimeras nonchimeras.fasta


### Relabel OTUs.

 $VSEARCH --fastx_filter nonchimeras.fasta \
     --threads $THREADS \
     --sizein \
     --fasta_width 0 \
     --relabel OTU_ \
     --fastaout otus.fasta

 
### map sequences to OTU
# Map all .fasta (before 97%) to otus.fasta (after chimera detection)
# Check the representativeness of the most abundant OTU of each sample (before excluding the smaller underrepresented sequences).
 $VSEARCH --usearch_global all.fasta \
     --threads $THREADS \
     --db otus.fasta \
     --id 0.97 \
     --strand plus \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --qmask none \
     --dbmask none \
     --otutabout GS19_GC05_GC25-GS21GC09_all.otutab.txt 


#### Sort OTU table numberically 
# Sort the dataset on highest abundance

 sort -k1.5n GS19_GC05_GC25-GS21GC09_all.otutab.txt > GS19_GC05_GC25-GS21GC09_all.otutab.sorted.txt


## LULU cleanup
# make blast db in terminal using last .fasta file from "relabel OTUs" step

makeblastdb -in otus.fasta -dbtype nucl

blastn -db otus.fasta -outfmt '6 qseqid sseqid pident' -out LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 -query otus.fasta -num_threads 3

# Run LULU.R script and change "XX" to unique name used in output from mapping sequences to OTUs 
# Looking for OTUs that are considered as 'NEW', but actually consist of a sequencing error. 
# Searches for similarity where one is always a higher representative within a sample than the other, but are present at the same time.
# Set at 95% similarity.

# RUN LULU.r in R
R LULU.r

### Remove sequences from fasta file that got removed using LULU
#(in Python)

python filter_seq_by_OTUtable.py otus.fasta GS19_GC05_GC25-GS21GC09_table_curated.tsv > OTUs_curated.fasta


## CREST Classifying
# Based on the last common ancestor. Blast your sequences, and compares it to the data in silvamod138.
# Using CREST4 which is installed locally
# activateconda environment to run CREST

conda activate condacrest4

#Run CREST4
crest4 -f OTUs_curated.fasta

# Rename the OTU files to have unique name for datasets you are working on
#" **XX_OTU_table.csv "

