# download viral protein from ncbi
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
gunzip viral* -c |cat - > virus.fa

# create usearch index of viruses
usearch -makeudb_usearch virus.fa -output idx.udb

# split hhmer peptides into chunks
awk
L=$(( $(wc -l <LANDALE.pep) / 10 )
split -l $L virus.sl.fa.
