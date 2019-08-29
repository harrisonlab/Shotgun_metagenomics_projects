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

# local alignment
for f in out.pep*; do
  for f in pep*; do
  usearch -ublast $f -db ublast.udb -evalue 1e-7 -maxaccepts 1 -userout $f.hits\
  -userfields query+target+ql+qs+ts+alnlen+id+evalue+bits 
  # usearch -usearch_global $f -db idx.udb -id 0.25 -blast6out ${f}.out; 
done


# deduplicate results
awk -F"\t" '{split($1,N,"_");print N[1]"_"N[2]}'} *.hits|sort|uniq > viral_peptides.txt


#### NEW IDEA ####
grep -E "Virus|virus|Viral|viral|Virion|virion|phage|Phage" ~/pipelines/common/resources/pfam/names.txt > all_virus.txt

```(R)
library(data.table)
dom <- fread("LANGDALE.domains")
pf <- fread("~/pipelines/common/resources/pfam/all_virus.txt",header=F)
dd <- pf[dom,on=c("V2"="V1")]
fwrite(as.data.table(unique(sub("_[0-9]*$","",dd[complete.cases(dd),i.V2],perl=T))),"all_viral_domains.txt")
```