# use ORFfinder to find - bacterial ORFS
ORFfinder -in test.fa -g 11 -s 1 -ml 180 -n true -outfmt 3 -out temp.out
grep ">" temp.out |
sed 's/.*_c/c/'|
awk -F":" '{count++;if ($3>$2){d="f"}else{d="r";t=$2;$2=$3;$3=t};p=1;$3=$3+1;if(($3-$2)%3==0){p=0};print count,$1,$2,$3,d,p,"program","v1.0"}' OFS="\t" > t2.out

# generate anvio database - 2 methods depending whether there are gene calls already available

## add header for first method
sed -i '1 i\gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\t\source\tversion' test.out
anvi-gen-contigs-database -f test.fa -o test.db -n "test database"  --external-gene-calls test.out

## no header for second method
anvi-gen-contigs-database -f test.fa -o test.db -n "test database" --skip-gene-calling 
## import ORFs into sqlite database (use  Control-v <TAB> for separator)
sqlite3 -separator ' ' test.db  ".import test.out genes_in_contigs"
## test sql import has worked
sqlite3 test.db "select * from genes_in_contigs limit 10"
## extract and translate 
anvi-get-sequences-for-gene-calls -c test.db -o test.prot.fa
