# The following was performed in Unix
# Retrieving SRA FASTQ files and basic filtering using NCBI SRA Toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software). This will download the paired-end sequences in 2 FASTQ files for forward and reverse reads.
for i in `cat SRA_HMP.txt`; do fasterq-dump $i --skip-technical --split-3 --min-read-len 50 --outdir $i -e 36; done

# Joining reads using Fastq-join (https://github.com/brwnj/fastq-join).
for i in `cat SRA_HMP.txt`; do fastq-join ${i}/${i}_1.fastq ${i}/${i}_2.fastq -o ${i}/${i}_%.fastq; done

# Quality filtering using PRINSEQ++ and converting FASTQ to FASTA (https://github.com/Adrian-Cantu/PRINSEQ-plus-plus). This cleans out the low quality reads.
for i in */*join.fastq; do prinseq++ -fastq $i -min_qual_mean 20 -ns_max_n 0 -derep -trim_qual_right=20 -lc_entropy -min_len 50 -threads 36 -out_format 1 -out_name $i; done

# Taxonomic profiling using MetaPhyler (http://metaphyler.cbcb.umd.edu/). This step identifies the phylogenetic marker genes present, classify them and compute the taxonomy profile.
for i in `cat SRA_HMP.txt`; do runMetaphyler.pl ${i}_join.fastq_good_out.fasta blastn $i 1; done

# Use DIAMOND (https://github.com/bbuchfink/diamond) to align the FASTA files to the reference database.
for i in `cat SRA.txt`; do diamond blastx -d diamond_ref/reference -q $i -o ${i}_matches.m8 -p 36; done

# The follwoing was performed in R
library(plyr)

tab=read.table("diamond_results",header=F)

colnames(tab)[1]="query.acc."

colnames(tab)[2]="reference.acc."

tab$reference.acc.="BG"

tab$seq=sub("..$","",tab$query.acc.)

tab2=tab[!duplicated(tab[,"seq"]),]

tab2$SRA=sub("\\..*","",tab2$seq)

df <- count(tab2, c('SRA','reference.acc.'))

taxa=read.table("result_taxa.txt",header=F)

colnames(taxa)=c("SRA","phyla")

final=merge(df, taxa, all.x = TRUE)

metadata=read.csv("metadata.csv",header=T)

final2=merge(final, metadata, all.x = TRUE)

final2$norm_count=((final2$freq)*1000/final2$phyla)

final2 <- na.omit(final2)

final2$trans_norm_count=log10(final2$norm_count)

write.csv(final2,"betaglucuronidase.csv")


