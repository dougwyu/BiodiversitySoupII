#!/bin/bash
set -e
set -u
set -o pipefail
##########################################################################################################
##########################################################################################################
# script for metabarcoding: fastp to Begum pipeline
# tested on macOS Catalina (zsh)
##########################################################################################################
###########################################################################################################

# I make extensive use of GNU parallel to run on commands on multiple files at once. Sending parallel commands to the terminal is a bit difficult. Usually, it works best to select the whole parallel command before sending to terminal. However, if the command is particularly long, it might work better to send each line to terminal, and the terminal concatenates themto build up the whole command. Which works better seems to depend on one's shell and also the particulars of the terminal being used.
# Begum (and DAMe before it) requires a set of support files and folders in particular locations:
    # data/Primers_COILeray.txt # primer sequences
    # data/Tags_biosoupII_COI.txt # tag sequences
    # data/PSinfo_biosoupII_COI{A,B,C,D,E,F,G,H}.txt # sample name (e.g. sample name Hhmlbody, the associated tags, and PCR replicate ('pool'))
    # data/seqs/folder{A,B,C,D,E,F,G,H}/pools_info.txt # fastq pathname for each PCR replicate ('pool')
    # data/seqs/folder{A,B,C,D,E,F,G,H}/pool{1,2,3} # folders to hold the individual PCR replicate fastq files


# set clustering and taxonomic assignment parameters
OTUSIM=97 # OTU similarity in percentage (0-100, e.g. 97)
echo "OTU similarity percentage is" ${OTUSIM}"%." #
ARTHMINPROB=0.8 # minimum boostrap support of assignment to Arthropoda for vsearch --sintax

# set pathnames

# edit $HOMEFOLDER to your folder holding your scripts, data, and output folders
HOMEFOLDER="/Users/Negorashi2011/Xiaoyangmiseqdata/BiodiversitySoupII_repo/"  # do not have a ~ in the path
echo "Home folder is" ${HOMEFOLDER}

# fastq files in $HOMEFOLDER/data/seqs/
# analysis outputs in $HOMEFOLDER/analysis/
SEQS="data/seqs/"
ANALYSIS="analysis/"
# Begum parameter input files in $HOMEFOLDER/data/

# set DAMe and Begum pathnames and check that they work
DAME="/usr/local/bin/DAMe/bin/" # pathname for the binaries on DAMe
    echo "The DAMe binaries are in" ${DAME}
    python3 --version # tested with 3.8.3
    python3 ${DAME}convertToUSearch.py -h # should see help file
BEGUM="/Users/Negorashi2011/src/Begum/src/" # pathname for the binaries on Begum
    echo "The Begum binaries are in" ${BEGUM}
    python2 --version # tested with 2.7.17
    python2 ${BEGUM}Begum.py sort -h # should see help file
    python2 ${BEGUM}Begum.py filter -h # should see help file

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
# confirm that you are in seqs/ folder

#############################################################################################
#### Pre-process fastq files:  adapter removal, trimming, merging, subsampling, moving to Begum folders

# Read in sample list and make a bash array of the sample names (e.g. A1_S1)
find * -maxdepth 0 -name "*_L001_R1_001.fastq.gz" > samplelist.txt  # find all files ending with _L001_R1_001_fastq.gz
sample_info=samplelist.txt # put samplelist.txt into variable
sample_names=($(cut -f 1 "${sample_info}" | uniq)) # convert variable to array
    echo "${sample_names[@]}" # echo all array elements
    echo ${sample_names[1]} # echo first array element
echo "There are" "${#sample_names[@]}" "samples that will be processed:  " "${sample_names[@]}" # echo number of elements in the array

# In the original April run, we observed that Tag27_Tag27 in PCRA_pool2 failed. We resequenced this sample in the November MiSeq run, and we concatenated those sequences to the original April run fastq files. However, to make things fair, we first checked the number of seqs in the other samples in PCRA_pool2, and we subsampled the November run to a similar number before concatenating.
# Note: The concatenation has already been done. There is no need to run this code. It is here as documentation.
    # mean number of totseqs in the April PCRApool2 was 31497 (not counting Tag27_Tag27, which failed)
    # Tag27_Tag27 in the November PCRApool2 was 345853 totseqs
    # 31497/345853 = 9.1% # November run datasize needs to be reduced to 9.1% of its original size
# use seqkit sample to reduce the November fastq to 9.1% of original size and concatenate to the April fastq
    # cd /Users/Negorashi2011/Xiaoyangmiseqdata/MiSeq_20171121/data/seqs
    # seqkit sample -s 100 -p 0.091 A2_S2_L001_R1_001.fastq.gz | gzip > A2_S2_9.1_L001_R1_001.fastq.gz
    # seqkit sample -s 100 -p 0.091 A2_S2_L001_R2_001.fastq.gz | gzip > A2_S2_9.1_L001_R2_001.fastq.gz
    # zcat A2_S2_9.1_L001_R1_001.fastq.gz | less
    # zcat A2_S2_9.1_L001_R2_001.fastq.gz | less
    # echo ${HOMEFOLDER}${SEQS}
    # cp A2_S2_9.1_L001_R1_001.fastq.gz ${HOMEFOLDER}${SEQS}
    # cp A2_S2_9.1_L001_R2_001.fastq.gz ${HOMEFOLDER}${SEQS}
    # cd ${HOMEFOLDER}${SEQS}
    # cat A2_S2_L001_R1_001.fastq.gz A2_S2_9.1_L001_R1_001.fastq.gz > A2_S2_L001_R1_001_comb.fastq.gz
    # cat A2_S2_L001_R2_001.fastq.gz A2_S2_9.1_L001_R2_001.fastq.gz > A2_S2_L001_R2_001_comb.fastq.gz
    # rm A2_S2_L001_R1_001.fastq.gz
    # rm A2_S2_L001_R1_001.fastq.gz
    # rm A2_S2_9.1_L001_R1_001.fastq.gz
    # rm A2_S2_9.1_L001_R2_001.fastq.gz
    # mv A2_S2_L001_R1_001_comb.fastq.gz A2_S2_L001_R1_001.fastq.gz # Date Modified is 6 March 2020
    # mv A2_S2_L001_R2_001_comb.fastq.gz A2_S2_L001_R2_001.fastq.gz # Date Modified is 6 March 2020

# To run a loop interactively, select the entire loop and send to terminal.  Do not ctrl-Enter each line because this can send a command N times, as many lines as the command covers. So if the command wraps over 2 lines, ctrl-Entering each line sends the whole command twice, causing the command to be run twice per loop.

# 1. Use fastp to trim adapters. trim bad nucleotides, and merge reads
# fastp -h # 0.20.1
# takes ~6 minutes on a 2.8GHz Intel i7
for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array.  So loop over all samples
do
      sample_prefix="$( basename $sample "_L001_R1_001.fastq.gz")"
      echo ${sample_prefix}
      # fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
      fastp --merge -i ${sample_prefix}_L001_R1_001.fastq.gz -I ${sample_prefix}_L001_R2_001.fastq.gz --merged_out ${sample_prefix}_mergedtrimmed.fq.gz --html ${sample_prefix}fastp.html # -o ${sample_prefix}.R1.fq.gz -O ${sample_prefix}.R2.fq.gz
done

# clean up
mkdir fastp_outputs/
mv *fastp.html fastp_outputs/
rm fastp.json

# 2. count number of reads in merged fastq files and inspect output to find a common num_seqs value that is roughly equal across all samples.
ls {A,B,C,D,E,F,G,H}{1,2,3}*_mergedtrimmed.fq.gz # view file sizes
ls {A,B,C,D,E,F,G,H}{1,2,3}*_mergedtrimmed.fq.gz | wc -l # 24 files
seqkit --threads 4 stats *_mergedtrimmed.fq.gz > mergedtrimmed_read_counts.txt
cat mergedtrimmed_read_counts.txt # view num_seqs. I choose a value of 350 000, which all PCRs except PCR_F exceed (PCR_F files have around 300K seqs), but PCR_F is matched by its technical replicate PCR_G. Note that this value throws out a lot of data, which one normally would not do, but for comparison purposes, we need to equalise file sizes.

# subsample sequence files to the same size (350 000 reads).
parallel -j 3 -k "seqkit sample -s 100 -n 350000 {1}{2}*_mergedtrimmed.fq.gz | gzip > {1}{2}_350K_mergedtrimmed.fq.gz" ::: A B C D E F G H ::: 1 2 3
    # setting the random seed to -s 100 in order to make results reproducible
ls {A,B,C,D,E,F,G,H}{1,2,3}*_350K_mergedtrimmed.fq.gz | wc -l # 24 files
ls {A,B,C,D,E,F,G,H}{1,2,3}*_350K_mergedtrimmed.fq.gz # Should all have 40-60M filesizes. However, the parallel command occasionally randomly fails for some files, so i look for any output files with 20 byte file sizes, remove them, and run just those again by adjusting the options in the parallel command. For instance, in one run, E1 and C2 failed, so i removed the 20-byte 350K output files and ran again for just E1 and C2 (easier to just adjust the parallel command, even though i'm only running one file, e.g. parallel "cmd" ::: E ::: 1)
# if failure happens regularly, set -j 1, to run one file at a time. This is slower, but surer.
    # rm E1_350K_mergedtrimmed.fq.gz
    # rm C2_350K_mergedtrimmed.fq.gz
    # parallel "seqkit sample -s 100 -n 350000 {1}{2}*_mergedtrimmed.fq.gz | gzip > {1}{2}_350K_mergedtrimmed.fq.gz" ::: E ::: 1
    # parallel "seqkit sample -s 100 -n 350000 {1}{2}*_mergedtrimmed.fq.gz | gzip > {1}{2}_350K_mergedtrimmed.fq.gz" ::: C ::: 2

ls {A,B,C,D,E,F,G,H}{1,2,3}*_350K_mergedtrimmed.fq.gz | wc -l # 24 files
ls {A,B,C,D,E,F,G,H}{1,2,3}*_350K_mergedtrimmed.fq.gz # should all have 50-60M filesizes except for F3

# 3. place libraries in folder/pool structure (see folder hierarchy in seqs/)
# there should also be a pools_info.txt file in each folder, which i made by hand
parallel mv {1}{2}_350K_mergedtrimmed.fq.gz folder{1}/pool{2} ::: A B C D E F G H ::: 1 2 3
ls folder{A,B,C,D,E,F,G,H}/pool{1,2,3}

# 4. remove original mergedtrimmed files
ls {A,B,C,D,E,F,G,H}{1,2,3}*_mergedtrimmed.fq.gz # should not display any _350K_mergedtrimmed files
ls {A,B,C,D,E,F,G,H}{1,2,3}*_mergedtrimmed.fq.gz | wc -l
# remove "--dryrun" and run
    parallel --dryrun "rm {1}{2}*_mergedtrimmed.fq.gz" ::: A B C D E F G H ::: 1 2 3
ls {A,B,C,D,E,F,G,H}{1,2,3}*_mergedtrimmed.fq.gz # should not display any _mergedtrimmed files

#############################################################################################
#### Begum - using PCR replicates and read numbers to filter out bad reads

# 1. Begum.py sort, takes ~2 hours on a 2.8GHz Intel i7 if run with -j 2
python2 ${BEGUM}Begum.py sort -h
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder

parallel -j 4 --progress "cd ${HOMEFOLDER}${SEQS}folder{1}; python2 ${BEGUM}Begum.py sort -p ${HOMEFOLDER}data/Primers_COILeray.txt -t ${HOMEFOLDER}data/Tags_biosoupII_COI.txt -s ${HOMEFOLDER}data/PSinfo_biosoupII_COI{1}.txt -l pools_info.txt -pm 2 -tm 1 -o PCR_{1}" ::: A B C D E F G H


# 2. Begum.py filter, takes a < 1 min per filter run
python2 ${BEGUM}Begum.py filter -h
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder

# filter all outputs by min 1-3 PCRS and min 1-4 copies per PCR
parallel -j 6 -k "cd ${HOMEFOLDER}${SEQS}; python2 ${BEGUM}Begum.py filter -d folder{1} -i folder{1}/PCR_{1} -s ${HOMEFOLDER}data/PSinfo_biosoupII_COI{1}.txt -p {2} -m {3} -l 300 -o Filter_min{2}PCRs_min{3}copies" ::: A B C D E F G H ::: 0.3 0.6 0.9 ::: 1 2 3 4

# change output filenames
cd ${HOMEFOLDER}${SEQS}
# rename for macOS is installed via homebrew
rename 's/min0.3PCR/min1PCR/' folder{A,B,C,D,E,F,G,H}/*.fna
rename 's/min0.6PCR/min2PCR/' folder{A,B,C,D,E,F,G,H}/*.fna
rename 's/min0.9PCR/min3PCR/' folder{A,B,C,D,E,F,G,H}/*.fna

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder

# 3. Organise Begum filter outputs

# # create text file BegumFilters, one line for each filter run above. Adjust as required
rm BegumFilters # remove any previous versions of the file
parallel -k "echo Filter_min{1}PCRs_min{2}copies >> BegumFilters" ::: 1 2 3 ::: 1 2 3 4
parallel -k echo :::: BegumFilters # should have only these 12 lines:
    # Filter_min1PCRs_min1copies
    # Filter_min1PCRs_min2copies
    # Filter_min1PCRs_min3copies
    # Filter_min1PCRs_min4copies
    # Filter_min2PCRs_min1copies
    # Filter_min2PCRs_min2copies
    # Filter_min2PCRs_min3copies
    # Filter_min2PCRs_min4copies
    # Filter_min3PCRs_min1copies
    # Filter_min3PCRs_min2copies
    # Filter_min3PCRs_min3copies
    # Filter_min3PCRs_min4copies

# make directories for Begum filter outputs
parallel mkdir folder{1}/{2}_{1} ::: A B C D E F G H :::: BegumFilters

# move Begum filter outputs into corresponding directories
parallel mv folder{1}/{2}.fna folder{1}/{2}_{1}/ ::: A B C D E F G H :::: BegumFilters

# 4. Convert Begum filter output files to usearch header format. The output fasta files (*.fna) have been filtered for erroneous sequences (at the different Begum stringency levels), and we now need one function from the original DAMe pipeline to put the number of sequences back into the header line. (e.g. '>mmmmbody Tag9.Tag9_Tag27.Tag27_Tag34.Tag34 61_27_24' gets converted to '>mmmmbody;size=112'). We plan to move this function to Begum.

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
python3 ${DAME}convertToUSearch.py -h
parallel -k "cd folder{1}/{2}_{1}; \
    python3 ${DAME}convertToUSearch.py -i {2}.fna -lmin 300 -lmax 330 -u; \
    seqkit replace -is -p \"n\" -r \"\" FilteredReads.forusearch.fna > FilteredReads.forusearch_noN.fna; \
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters
    # remove the Ns from the ends of the sequences put there by ${DAME}convertToUSearch.py

# At this point, the fasta files can be processed by a variety of metabarcoding pipelines. We continue below with a vsearch pipeline.

#############################################################################################
#### Generate OTUs and OTU tables

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder

# 1. use the vsearch --cluster_size method to generate OTUs, and use --usearch_global + --otutabout to generate OTU tables by mapping the FilteredReads.forusearch_noN.fna reads to the OTU representative sequences. size= information automatically read
# This step takes ~ 2 hrs, entirely because of the Filter_min1PCRs_min1copies data, which we normally do not create.
parallel -j 4 -k --progress "cd folder{1}/{2}_{1};\
    vsearch --derep_fulllength FilteredReads.forusearch_noN.fna --sizein --sizeout --fasta_width 0 --threads 0 --output FilteredReads_derep.fas --relabel "OTU";\
    vsearch --sortbysize FilteredReads_derep.fas --output FilteredReads_derep_sorted.fas;\
    vsearch --uchime_denovo FilteredReads_derep_sorted.fas --nonchimeras FilteredReads_derep_sorted_nonchimeras.fas;\
    vsearch --cluster_size FilteredReads_derep_sorted_nonchimeras.fas --sizein --sizeout --id 0.97 --sizeorder --centroids FilteredReads_derep_sorted_nonchimeras_vsearch97.fas;\
    vsearch --sortbysize FilteredReads_derep_sorted_nonchimeras_vsearch97.fas --output FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted.fas;\
    vsearch --threads 6 --usearch_global FilteredReads.forusearch_noN.fna --db FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted.fas --id .97 --otutabout "table_BioSoupII_{1}_${OTUSIM}.txt";\
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters


# At this stage, the user might choose to run the R package {lulu} to merge parent-child OTUs that have not successfully been clustered by vsearch. We did not run lulu because all our OTUs were added to all our samples, which means that we cannot use OTU co-occurrence in subsets of samples as a signal that two similar sequences are parent-child OTUs.
# # 2. generate matchlist.txt for lulu (compare OTUs pairwise for % similarities)
# cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
# parallel -j 4 -k "cd folder{1}/{2}_{1}/; \
#     vsearch --usearch_global FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted.fas --db FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted.fas --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10; \
#     cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters
#
# # 3. lulu generate lulu-curated OTU table
# parallel -j 4 -k "cd folder{1}/{2}_{1}/; \
#     Rscript --vanilla --verbose ${HOMEFOLDER}scripts/LULU_20201109.R "table_BioSoupII_{1}_${OTUSIM}.txt"; \
#     rm lulu.log_*; \
#     rm matchlist.txt; \
#     rm "table_BioSoupII_{1}_${OTUSIM}.txt"; \
#     cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters
#
# # 4. remove any child sequences from the OTU list
# parallel -k "cd folder{1}/{2}_{1}/; \
#     seqkit replace -p \";\" -r \" ;\" FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted.fas > FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_.fas; \
#     seqtk subseq FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_.fas <(cut -f 1 table_BioSoupII_97_lulu.txt) > FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_lulu_{1}.fas
#     cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters
#         # seqkit replace to add space before the semicolon
#         # <(cut -f 1 table_300test_${sample}_${sim}_lulu.txt) produces the list of sequences to keep


# If you have not run lulu, this command changes filenames for downstream processing, by adding the folder to the end of FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted.fas
# FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_B.fas
parallel -k "cd folder{1}/{2}_{1}/; \
    mv FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted.fas FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.fas; \
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters

# 5. remove vsearch working files
parallel -k "cd folder{1}/{2}_{1}/; \
    rm FilteredReads.forusearch.fna; \
    rm FilteredReads.forusearch_noN.fna; \
    rm FilteredReads_derep.fas; \
    rm FilteredReads_derep_sorted.fas; \
    rm FilteredReads_derep_sorted_nonchimeras.fas; \
    rm FilteredReads_derep_sorted_nonchimeras_vsearch97.fas; \
    # rm FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_.fas; \ # if you ran lulu, uncomment this
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters

# 6. Taxonomic assignment of OTU sequences
# At this stage, we are only checking taxonomies down to Arthropoda, so we use something simple like vsearch --sintax to the MIDORI COI database: MIDORI_UNIQUE_20180221_COI_SINTAX.fasta.zip, available at http://www.reference-midori.info/download.php#.  The file is large (744 MB), so we unzip before use, and discard the unzipped file after use.

cd ${HOMEFOLDER}scripts/MIDORI # cd into the sequence folder and unzip the MIDORI file before using. The zipped file remains, so there is no need to rezip afterwards
# gunzip < MIDORI_UNIQUE_20180221_COI_SINTAX.fasta.gz > MIDORI_UNIQUE_20180221_COI_SINTAX.fasta
MIDORIDB=${HOMEFOLDER}scripts/MIDORI/MIDORI_UNIQUE_20180221_COI_SINTAX.fasta
seqkit head ${MIDORIDB} # check that the unzipped file is available

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
# this takes ~5 hrs, mostly b/c of the Filter_min1PCRs_min1copies data. The other datasets take a few mins each
# tried with --strand both, which is slower, and all hits were + strand
parallel -j 1 -k --progress "cd folder{1}/{2}_{1}/; \
    vsearch --threads 7 --sintax FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.fas \
    --db ${MIDORIDB} --strand plus --sintax_cutoff ${ARTHMINPROB} \
    --tabbedout FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.out;
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters


# remove unzipped MIDORI database
ls ${MIDORIDB}
rm ${MIDORIDB}
ls ${MIDORIDB} # should say "No such file or directory"

# 7. filter OTU representative sequences to keep only Arthropoda OTUs with prob >= ARTHMINPROB (set to 0.80).
## if you are running on macOS with the zsh shell, using awk with GNU parallel is awkward. first put the awk function inside a function and export it. then run the function with env_parallel, not parallel
DOAWK() {
   FOLDER=$1
   awk -F'[;,\t]' '$4 ~ /Arthropoda/  { print }' FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_$FOLDER.out > FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_$FOLDER.out
}; export -f DOAWK

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
# use env_parallel to run user-defined functions with GNU parallel
# to install env_parallel run these two commands; orig instructions at https://stackoverflow.com/a/27559537
    # echo '. `which env_parallel.zsh`' >> $HOME/.zshenv
    # source HOME/.zshenv
env_parallel -k --dryrun "cd folder{1}/{2}_{1}; DOAWK {1}" ::: A B C D E F G H :::: BegumFilters
env_parallel -k "cd folder{1}/{2}_{1}; DOAWK {1}" ::: A B C D E F G H :::: BegumFilters

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
parallel -k "cd folder{1}/{2}_{1}/; \
    seqtk subseq FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.fas <(cut -f 1 FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.out) > FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.fas; \
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters

parallel "ls folder{1}/{2}_{1}" ::: A B C D E F G H :::: BegumFilters


## if you are running under a bash shell, everything can be done in a single parallel command.
# 7 orig. filter OTU representative sequences to keep only Arthropoda OTUs with prob >= ARTHMINPROB (set to 0.80).
# cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
# parallel -k "cd folder{1}/{2}_{1}/; \
#     awk '$(printf '$4') ~ /Arthropoda/  { print }' FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.out > FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.out; \
#     seqtk subseq FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.fas <(cut -f 1 FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.out) > FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.fas; \
#     cd ${HOMEFOLDER}${SEQS}" ::: B :::: BegumFilters
# The awk cmd filters in sintax output for rows where 4th column contains "Arthropoda" (which are the ones with likelihood > 0.80)
# if awk command throws a syntax error, try using $(printf '$4') in place of $($4)
# seqtk subseq filters for OTUs where OTU number is in the Arthropoda-filtered sintax output

# sanity check. number of rows and number of seqs should be the same
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
parallel -k "seqkit stats folder{1}/{2}_{1}/FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.fas; \
cat folder{1}/{2}_{1}/FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.out | wc -l" ::: A B C D E F G H :::: BegumFilters

# 8. filter OTU table table_BioSoupII_97_Arthropoda_{1}.txt to keep only only Arthropoda OTUs with prob >= ARTHMINPROB (set to 0.80)
    # gsed '1i <string>' file > outfile.  1i adds 'OTU_ID' to first line;  sort, add, then redirect
    # https://shapeshed.com/unix-join/
    # -t "$(printf '\t')" # to generate a tab character on the fly within join
    # only joins lines with matching first column fields
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
parallel -k "cd folder{1}/{2}_{1}/; \
    gsed 's/#OTU ID/OTU_ID/' table_BioSoupII_{1}_97.txt > tmp; \
    mv tmp table_BioSoupII_{1}_97.txt; \
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters

parallel -k "cd folder{1}/{2}_{1}/; \
    gsed '1i OTU_ID' <(cut -d ';' -f 1 FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.out | sort) > Arthropoda_OTUs.txt
    join Arthropoda_OTUs.txt table_BioSoupII_{1}_97.txt > table_BioSoupII_97_Arthropoda_{1}.txt; \
    gsed 's/ /\t/g' table_BioSoupII_97_Arthropoda_{1}.txt > table_BioSoupII_97_Arthropoda_{1}2.txt; \
    mv table_BioSoupII_97_Arthropoda_{1}2.txt table_BioSoupII_97_Arthropoda_{1}.txt; \
    rm Arthropoda_OTUs.txt; \
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters


# 9. vsearch OTUs to mock-soup fasta database

# vsearch variant
OTUSIMvsearch=94 # to allow in some low-pct matches for inspection (blastn does this anyway)
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
parallel -k "cd folder{1}/{2}_{1}/; \
    vsearch --usearch_global FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.fas --db ${HOMEFOLDER}data/MTB/S1_MTBFAS.fasta --id .${OTUSIMvsearch} --blast6out table_BioSoupII_{1}_Arthropoda.vsearchMTB.txt; \
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters

# blastn variant
# run once to make a blast database
    # cd ${HOMEFOLDER}data/MTB # cd into the sequence folder
    # makeblastdb -in S1_MTBFAS.fasta -dbtype nucl -parse_seqids
echo .${OTUSIM} # .97
# OTUSIM=97 # if variable not present
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
parallel -k "cd folder{1}/{2}_{1}/; \
    blastn -db ${HOMEFOLDER}data/MTB/S1_MTBFAS.fasta -num_threads 7 -qcov_hsp_perc .90 -perc_identity .${OTUSIM} \
    -max_target_seqs 1 -outfmt 6 -out table_BioSoupII_{1}_Arthropoda.blastnMTB.txt \
    -query FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_{1}.fas; \
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H :::: BegumFilters

############################################################################################################
#### Single Pool analysis

# The output header lines of ${DAME}convertToUSearch.py (e.g. in Filter_min1PCRs_min1copies.fna) have the format:
    # >mmmmbody Tag9.Tag9_Tag27.Tag27_Tag34.Tag34 61_27_24
    # seqkit replace to replace \tTag[0-9]+.Tag[0-9]+_Tag[0-9]+.Tag[0-9]+_Tag[0-9]+.Tag[0-9]+\t with ;size=
    # seqkit replace to parse 61_27_24 into three files, i.e.
        # >mmmmbody;size=61  # for pool1
        # >mmmmbody;size=27  # for pool2
        # >mmmmbody;size=24  # for pool3
    # seqkit grep to remove all size=0 seqs
# this creates 3 fna files for each Filter_min1PCRs_min1copies.fna, one for each pool
# then generate singlepool OTU tables using vsearch --otutabout to compare these reads against the
    # representative OTUs for that folder (FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.fas)

cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
# seqkit grep -n search by full name, -v invert-match, -r --use-regexp
# this parallel command needs to be entered one line at a time
parallel -j 1 -k "cd folder{1}/Filter_min1PCRs_min1copies_{1}; \
    echo {1};\
    seqkit -j 4 replace -i -p \"\tTag[0-9]+.Tag[0-9]+_Tag[0-9]+.Tag[0-9]+_Tag[0-9]+.Tag[0-9]+\t\" -r \";size=\" Filter_min1PCRs_min1copies.fna > Filter_min1PCRs_min1copies_notag.fna; \
    seqkit -j 4 replace -i -p \"=([0-9]+)_([0-9]+)_([0-9]+)\" -r '=$(printf '$1')' Filter_min1PCRs_min1copies_notag.fna > Filter_min1PCRs_min1copies_pool1.fna;\
    seqkit -j 4 replace -i -p \"=([0-9]+)_([0-9]+)_([0-9]+)\" -r '=$(printf '$2')' Filter_min1PCRs_min1copies_notag.fna > Filter_min1PCRs_min1copies_pool2.fna;\
    seqkit -j 4 replace -i -p \"=([0-9]+)_([0-9]+)_([0-9]+)\" -r '=$(printf '$3')' Filter_min1PCRs_min1copies_notag.fna > Filter_min1PCRs_min1copies_pool3.fna;\
    seqkit -j 4 grep -n -v -r -p \"size=0\" Filter_min1PCRs_min1copies_pool1.fna > Filter_min1PCRs_min1copies_pool1_min1.fna;\
    seqkit -j 4 grep -n -v -r -p \"size=0\" Filter_min1PCRs_min1copies_pool2.fna > Filter_min1PCRs_min1copies_pool2_min1.fna;\
    seqkit -j 4 grep -n -v -r -p \"size=0\" Filter_min1PCRs_min1copies_pool3.fna > Filter_min1PCRs_min1copies_pool3_min1.fna;\
    cd ${HOMEFOLDER}${SEQS}" ::: A B C D E F G H ::: 1 2 3

parallel -j 1 -k "rm folder{1}/Filter_min1PCRs_min1copies_{1}/Filter_min1PCRs_min1copies_notag.fna" ::: A B C D E F G H
parallel -j 1 -k "rm folder{1}/Filter_min1PCRs_min1copies_{1}/Filter_min1PCRs_min1copies_pool{2}.fna" ::: A B C D E F G H ::: 1 2 3

# make otu table by matching the single pool reads against the OTUs generated for that folder
# using the representative OTUs that include the non-Arthropoda because we want to see as many OTUs as possible
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
parallel -j 1 -k "cd folder{1}/Filter_min1PCRs_min1copies_{1}; \
    echo folder{1} pool{2}; \
    vsearch --threads 6 --usearch_global Filter_min1PCRs_min1copies_pool{2}_min1.fna --db FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_{1}.fas --id .97 --otutabout table_BioSoupII_97_SnglPl_{1}{2}.txt" ::: A B C D E F G H ::: 1 2 3


############################################################################################################
#### Move output tables to ${HOMEFOLDER}${ANALYSIS}

## single pool outputs
cd ${HOMEFOLDER}${SEQS} # cd into the sequence folder
mkdir ${HOMEFOLDER}${ANALYSIS}singlepools/
parallel -k "cp folder{1}/Filter_min1PCRs_min1copies_{1}/table_BioSoupII_97_SnglPl_{1}{2}.txt ${HOMEFOLDER}${ANALYSIS}singlepools/" ::: A B C D E F G H ::: 1 2 3

# slightly change filename and change "#OTU ID" to OTU_ID (to fit downstream R code a bit better)
cd ${HOMEFOLDER}${ANALYSIS}singlepools/
parallel "mv table_BioSoupII_97_SnglPl_{1}{2}.txt table_BioSoupII_97_{1}{2}_SnglPl.txt" ::: A B C D E F G H ::: 1 2 3
parallel "gsed -i 's/#OTU ID/OTU_ID/g' table_BioSoupII_97_{1}{2}_SnglPl.txt" ::: A B C D E F G H ::: 1 2 3

# Begum outputs
# copy OTU, blast, and vsearch tables to analysis/ folder
cat ${HOMEFOLDER}${SEQS}BegumFilters # should show the twelve filter settings
# make analysis folders
    cd ${HOMEFOLDER}${ANALYSIS}
    parallel mkdir {} :::: ${HOMEFOLDER}${SEQS}BegumFilters
# copy tables to analysis folders
parallel "cp ${HOMEFOLDER}${SEQS}folder*/{1}_*/table_BioSoupII_97_Arthropoda_*.txt ${HOMEFOLDER}${ANALYSIS}{1}/" :::: ${HOMEFOLDER}${SEQS}BegumFilters

parallel "cp ${HOMEFOLDER}${SEQS}folder*/{1}_*/table_BioSoupII_*_Arthropoda.blastnMTB.txt ${HOMEFOLDER}${ANALYSIS}{1}/" :::: ${HOMEFOLDER}${SEQS}BegumFilters

parallel "cp ${HOMEFOLDER}${SEQS}folder*/{1}_*/table_BioSoupII_*_Arthropoda.vsearchMTB.txt ${HOMEFOLDER}${ANALYSIS}{1}/" :::: ${HOMEFOLDER}${SEQS}BegumFilters

parallel "cp ${HOMEFOLDER}${SEQS}folder*/{1}_*/FilteredReads_derep_sorted_nonchimeras_vsearch97_sorted_Arthropoda_*.out ${HOMEFOLDER}${ANALYSIS}{1}/" :::: ${HOMEFOLDER}${SEQS}BegumFilters

# END

# Go to 3_SoupII_ecological_analysis.Rmd and then 4_SoupII_singlepools_OTU_table.Rmd
