#!/bin/sh                                                                       

# ########################################################         
#     This BASH script was used to cluster OTUs using 
#          VSEARCH. Please makesure the sequencing quality
#     data is checked using programs like Qiime2. 
#    Please see the VSEARCH webpage for more information.
#    https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline
#    
# Kabin Xie
# 2018.8.20 
# 2019.9.6 updated
# Modified to process 1171-1196 soil samples data

THREADS=7
REF=/home/kabinxie/mystation/microbiome/gold.fasta
PERL=$(which perl)
VSEARCH=$(which vsearch)
SEQ_DIR=/home/kabinxie/mystation/microbiome/Cas-16S-seq
date

## sample ID and file name are listed in root_soil_sample-seq.info.vsearch.csv
## Make sure no - in sample id
##


mkdir vsearch.out
cd vsearch.out

# Process samples  
# The sample ID and seq file infomation retrived from csv file.                                                             

sample_file=~/mystation/microbiome/Cas-16S-seq/DADA2/root_soil_sample-seq.info.vsearch.csv
#sample_file=~/mystation/microbiome/Cas-16S-seq/DADA2/test.csv

sample_id=$(cat $sample_file | awk -F , '{print $1}')
  
for Sample in $sample_id; do
    
    f=$SEQ_DIR/$(grep ^$Sample  $sample_file | awk -F , '{print $2}')
    r=$SEQ_DIR/$(grep ^$Sample  $sample_file | awk -F , '{print $3}')

    echo
    echo ====================================
    echo Processing sample $Sample
    echo ====================================

    s=$Sample
    $VSEARCH --threads $THREADS --fastq_mergepairs $f  --reverse $r --fastq_minovlen 10  --fastq_maxdiffs 15   --fastqout $s.merged.fastq   --fastq_eeout

    # Remove f and r primers 2019.9.4           
    # using cutadapt -a AACMGGATTAGATACCCKG(For) -A CGTCATCCMCACCTTCCTCADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq   
    # The last 20nt removed due to low-quality
    # For unknown reason, the 
                                
    echo Removing adaptor
    cutadapt --cores=$THREADS -g AACMGGATTAGATACCCKG $s.merged.fastq -o $s.merged.trimmed1.fastq  
    cutadapt --cores=$THREADS -a GAGGAAGGTGKGGATGACG $s.merged.trimmed1.fastq -o $s.merged.trimmed2.fastq
    cutadapt -u -20 $s.merged.trimmed2.fastq  -o $s.merged.trimmed3.fastq

    echo
    echo Calculate quality statistics

    $VSEARCH --threads $THREADS \
        --fastq_eestats $s.merged.trimmed3.fastq \
        --output $s.stats
    echo
    echo Quality filtering
## set maxeln 550, make sure no mito sequences (470-490bp) are not filtered)

    $VSEARCH --threads $THREADS \
        --fastq_filter $s.merged.trimmed3.fastq \
        --fastq_maxee 1 \
        --fastq_minlen 250 \
        --fastq_maxlen 550 \
        --fastq_maxns 0 \
        --fastaout $s.filtered.fasta \
        --fasta_width 0

    echo
    echo Dereplicate at sample level and relabel with sample_n

    $VSEARCH --threads $THREADS \
        --derep_fulllength $s.filtered.fasta \
        --strand plus \
        --output $s.derep.fasta \
        --sizeout \
        --uc $s.derep.uc \
        --relabel $s. \
        --fasta_width 0   
done


echo Sum of unique sequences in each sample: $(grep -c "^>" *.derep.fasta)

# At this point there should be one fasta file for each sample                  
# It should be quality filtered and dereplicated.                               

echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

rm -f all.derep.fasta all.nonchimeras.derep.fasta
cat *.derep.fasta > all.fasta

echo
echo Dereplicate across samples and remove singletons

$VSEARCH --threads $THREADS \
    --derep_fulllength all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.derep.uc \
    --output all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" all.derep.fasta)

## unoise the reads 2019.09.04
echo
echo Unoise the sequences


echo
echo Precluster at 98% before chimera detection

$VSEARCH --threads $THREADS \
    --cluster_size all.derep.fasta \
    --id 0.98 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.preclustered.uc \
    --centroids all.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" all.preclustered.fasta)



echo
echo De novo chimera detection

$VSEARCH --threads $THREADS \
    --uchime_denovo all.preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras all.denovo.nonchimeras.fasta \

echo Unique sequences after de novo chimera detection: $(grep -c "^>" all.denovo.nonchimeras.fasta)

echo
echo Reference chimera detection

$VSEARCH --threads $THREADS \
    --uchime_ref all.denovo.nonchimeras.fasta \
    --db $REF \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras all.ref.nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" all.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

$PERL ../map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

$PERL ../map.pl all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

echo
echo Cluster at 97%, 99% and 100% and relabel with OTU_n, generate OTU table

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel OTU_ \
    --centroids all.otus.fasta \
    --otutabout all.otutab.txt

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 0.99 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel OTU_ \
    --centroids all.otus99.fasta \
    --otutabout all.otutab99.txt

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 1 \
    --strand plus \
    --minsize 8 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel OTU_ \
    --centroids all.otus100.fasta \
    --otutabout all.otutab100.txt
echo
echo Number of OTUs: $(grep -c "^>" all.otus*.fasta)

# Get the OTU of rice mitochondrion rDNA
$VSEARCH  --usearch_global  ~/mystation/microbiome/Os.Mito.16SrDNA.fasta -db all.otus.fasta  --id 0.99 -alnout mito-out.aln

echo
echo Done

date
















