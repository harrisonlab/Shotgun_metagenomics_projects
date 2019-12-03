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
sqlite3 -separator 'Control-v <TAB>' test.db  ".import test.out genes_in_contigs"
## test sql import has worked
sqlite3 test.db "select * from genes_in_contigs limit 10"
## extract/translate /import
anvi-get-sequences-for-gene-calls -c test.db -o test.prot.fa
mkfifo mytemp
java -jar ~/programs/bin/macse_v2.03.jar -prog translateNT2AA -gc_def 11 -seq test.prot.fa -out_AA mytemp &
tr -d '>'<mytemp|paste - - > x.import.fa
sqlite3 -separator 'Control-v <TAB>' test.db  ".import x.import.fa gene_amino_acid_sequences"
sqlite3 test.db "select * from gene_amino_acid_sequences limit 10"
## update the self table to tell it we have gene calls
sqlite3 test.db "update self set value = 1 where key = 'genes_are_called'"
