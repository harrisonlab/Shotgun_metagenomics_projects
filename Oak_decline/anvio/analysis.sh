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

# add functional annotations (I have these already - will look at how to import)
anvi-run-hmms -c test.db --num-threads 20
anvi-run-ncbi-cogs -c test.db --num-threads 20

# add alignment information
## bam files - need to be sorted and indexed (eugh)
/home/greg/emr-dstore1/Cluster/data/data2/scratch2/deakig/Oak/aligned/A723_NDME02225_HFFF7DMXX_L1.bam

samtools view -h /data/data2/scratch2/deakig/Oak/aligned/A723_NDME02225_HFFF7DMXX_L1.bam|
head -n 19345292|
samtools sort -O bam -o sub.bam
samtools index sub.bam

head names.txt -n 5000|awk -F" " '{$1="";print $0}' > fasta
anvi-profile -i sub.bam -c test.db -S "A723" --contigs-of-interest fasta

# argh the headers need to be free from extra characters
samtools view -h  sub.bam|sed 's/ flag.*len=[0-9]*//'|samtools view -S -b  >sub2.bam
samtools index sub2.bam

# sort bam files
for f in *.bam; do
PREFIX=$(echo $f|sed -e 's/\..*//')
sbatch --mem-per-cpu 2000M -c 10 \
~/pipelines/metagenomics/scripts/slurm/sub_bam_sort.sh \
10 /data/data2/scratch2/deakig/Oak/sorted $PREFIX /data/data2/scratch2/deakig/Oak/aligned/$f
done

# will need to remove extra characters from headers before indexing
samtools view -h A402_NDME02229_HFFF7DMXX_L1_1|sed -E 's/ flag=[0-9]* multi=[0-9]*\.[0-9]* len=[0-9]*//g'|samtools view -S -b  >A402_NDME02229_HFFF7DMXX_L1_1.bam

for f in *_1; do
samtools view -h $f|sed -E 's/ flag=[0-9]* multi=[0-9]*\.[0-9]* len=[0-9]*//g'|samtools view -S -b > $f.bam
samtools index $f.bam
done
