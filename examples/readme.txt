##GFF3 file: curl ftp://ftp.ensembl.org/pub/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh38.84.gff3.gz --output Homo_sapiens.GRCh38.84.gff3.gz
##Region to check: chr6:146483132-149217381 ##Chosen arbitrarily because it is around the YAP gene and doesn't have gene features that are split

##Pull out only chromosome 6 genes
zless Homo_sapiens.GRCh38.84.gff3.gz | grep -P "^6\t" | sed 's/^/chr/g' | gt gff3 -tidy -sort -retainids > Homo_sapiens.GRCh38_chr6.84.gff3

##Pull out features only in the subregion
../../gffkit/gffkit.py subgff -g Homo_sapiens.GRCh38_chr6.84.gff3 -r 146483132:149217381 | gt gff3 -tidy -sort -retainids > Homo_sapiens.GRCh38_chr6-subset.84.gff3

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr6.fa.gz .
seqkit subseq -r 146483132:149217381 chr6.fa.gz > hg38-chr6-subseq1.fa
seqkit seq --reverse --complement > hg38-chr6-subseq2.fa
