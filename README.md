# Biodiversity Soup II paper

See Wiki page for assumed folder structure.  

**This pipeline is designed for my laboratory to use, so clone and employ in your environment with care!**

1_SoupII_software_install_on_macos_20191010.sh  

This pipeline is tested for macOS and processes Illumina HiSeq/MiSeq files for metabarcoding. Software installation information is in this script.  

2_SoupII_fastp_to_Begum_Metabarcoding_pipeline.sh  
1. Create merged reads using fastp (Chen et al. 2018)  

2. Demultiplex, Filter, and Cluster the merged reads using the Begum pipeline, which is a rewrite of the DAMe pipeline (Zepeda-Mendoza et al. 2016, BMC Research Notes) (http://github.com/shyamsg/Begum)  
a)  sort.py to remove tags from merged reads and sum by tag pair to visualise tag jumping and PCR/sequencing success  
b)  filter.py to keep only reads that exceed Minimum_PCR and Minimum_copy_number_per_PCR thresholds, defined by the user  
c)  vsearch -uchime_denovo chimera removal  
d)  vsearch clustering and OTU table creation  
e)  vsearch --sintax taxonomic assignment on the MIDORI COI database  
f)  Filter OTU tables and OTU representative sequences fasta by taxonomy (keeping Arthropoda ≥ 0.80 probability)  

3_SoupII_ecological_analysis.Rmd  
1. Community analysis  
a) filter out artefactually small OTUs via phyloseq method (https://github.com/joey711/phyloseq)  
b) NMDS ordination to view results  
c) estimate proportion of OTUs that are missing (drop-outs) and are erroneous (false-positives, drop-ins) as a function of Begum filtering, PCR condition, and taxonomy  
d) analyse whether OTU size contains abundance information    

4_SoupII_singlepools_OTU_table.Rmd   
1. Tag-bias analysis  
a) Exploit the replicate PCRs using same and different tags and compare for differences in community composition  

Note that the Begum method works by independently PCRing each sample 3 times, us assumes that each set of samples (a 96-well plate) has been separately PCRd three times.  DAME keeps reads that appear in ≥M PCRs (e.g. ≥2 of 3 total PCRs) and in each PCR appears in ≥N copies (e.g. ≥4 copies per PCR). DAMe also works best if both the forward and reverse reads have been tagged using the same tag (what we call 'twin tags').
