# Biodiversity Soup II paper

See Wiki page for required folder structure.  

This folder contains all the files created after the BioSoupII pipeline has finished creating OTU tables and assigning taxonomies. If you want to run the pipeline from scratch, remove the following folders. The remaining files are the information files that DAMe/Begum uses to run the sorting and filtering and the reference files for blasting. 
```
cd BiodiversitySoupII_repo/  
# there should be three folders:  analysis/, data/, and scripts/
rm -rf analysis/Filter_min1PCRs_min1copies/  
rm -rf analysis/Filter_min2PCRs_min4copies/  
rm -rf analysis/Filter_min3PCRs_min3copies/  
rm -rf analysis/singlepools/  
rm data/seqs/BegumFilters  
rm data/seqs/mergedtrimmed_read_counts.txt  
rm data/seqs/samplelist.txt  
rm -rf data/seqs/fastp_outputs/  
rm -rf data/seqs/folder{A,B,C,D,E,F,G,H}/Filter_min1PCRs_min1copies_{A,B,C,D,E,F,G,H}  
rm -rf data/seqs/folder{A,B,C,D,E,F,G,H}/Filter_min2PCRs_min4copies_{A,B,C,D,E,F,G,H}  
rm -rf data/seqs/folder{A,B,C,D,E,F,G,H}/Filter_min3PCRs_min3copies_{A,B,C,D,E,F,G,H}  
rm -rf data/seqs/folder{A,B,C,D,E,F,G,H}/PCR_{A,B,C,D,E,F,G,H}_pool{1,2,3}.summaryCounts  
rm -rf data/seqs/folder{A,B,C,D,E,F,G,H}/PCR_{A,B,C,D,E,F,G,H}_pool{1,2,3}.tagInfo  
rm -rf data/seqs/folder{A,B,C,D,E,F,G,H}/pool{1,2,3}/*_350K_mergedtrimmed.fq.gz  
```



**These are the scripts used in our Biodiversity Soup II manuscript.**  

**1_SoupII_software_install_on_macos_20191010.sh**  

This pipeline is tested for macOS and processes Illumina HiSeq/MiSeq files for metabarcoding. Software installation information is in this script.  

**2_SoupII_fastp_to_Begum_Metabarcoding_pipeline.sh**  

1. Create merged reads using fastp (Chen et al. 2018)  

2. Demultiplex, Filter, and Cluster the merged reads using the Begum pipeline, which is a rewrite of the DAMe pipeline (Zepeda-Mendoza et al. 2016, BMC Research Notes) (http://github.com/shyamsg/Begum)  
a)  sort.py to remove tags from merged reads and sum by tag pair to visualise tag jumping and PCR/sequencing success  
b)  filter.py to keep only reads that exceed Minimum_PCR and Minimum_copy_number_per_PCR thresholds, defined by the user  
c)  vsearch -uchime_denovo chimera removal  
d)  vsearch clustering and OTU table creation  
e)  vsearch --sintax taxonomic assignment on the MIDORI COI database  
f)  Filter OTU tables and OTU representative sequences fasta by taxonomy (keeping Arthropoda ≥ 0.80 probability)  

**3_SoupII_ecological_analysis.Rmd**  

1. Community analysis  
a) filter out artefactually small OTUs via phyloseq method (https://github.com/joey711/phyloseq)  
b) NMDS ordination to view results  
c) estimate proportion of OTUs that are missing (drop-outs) and are erroneous (false-positives, drop-ins) as a function of Begum filtering, PCR condition, and taxonomy  
d) analyse whether OTU size contains abundance information    

**4_SoupII_singlepools_OTU_table.Rmd**   

1. Tag-bias analysis  
a) Exploit the replicate PCRs using same and different tags and compare for differences in community composition  

Note that the Begum method works by independently PCRing each sample 3 times.  Begum filters out erroneous reads by keeping only those that appear in ≥M PCRs (e.g. ≥2 of 3 total PCRs) with more than ≥N copies per PCR (e.g. ≥4 copies per PCR). Begum also works best if both the forward and reverse reads have been tagged using the same tag (what we call 'twin tags').  
