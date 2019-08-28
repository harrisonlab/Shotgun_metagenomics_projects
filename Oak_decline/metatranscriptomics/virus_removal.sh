# download viral protein from ncbi
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz
gunzip viral* -c |cat - > virus.fa

# create usearch index of viruses
usearch -makeudb_usearch virus.fa -output idx.udb

# single line protein file
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < LANGDALE.pep| tail -n +2 > LANGDALE.sl.pep

# split hhmer peptides into chunks
L=$(( $(wc -l <LANDGALE.sl.pep) / 10 )
split -l $L out.pep.

for f in out.pep*; do
  usearch -usearch_global $f -db idx.udb -id 0.25 -blast6out ${f}.out; 
done

