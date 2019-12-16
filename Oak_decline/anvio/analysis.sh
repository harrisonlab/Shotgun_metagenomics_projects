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
## extract/translate /import - the below will work but is appalingly memory hungary
#anvi-get-sequences-for-gene-calls -c test.db -o test.prot.fa
## may be better to write my own code - the t-sql below will extract all the required info
# then some perl to extract seq and revComp if required.

# rm mytemp1 mytemp2
# mkfifo mytemp1 mytemp2
# sqlite3 langdale.db 'SELECT start,stop,direction,sequence 
# FROM genes_in_contigs 
# LEFT JOIN contig_sequences 
# WHERE genes_in_contigs.contig=contig_sequences.contig' |
#perl -e '
#my $count=1;
#while(<>)
#{
#  chomp;
#  my @x = split /\|/,$_;
#  my $seq =substr $x[3],$x[0],($x[1]-$x[0]+1);
#  if ($x[2] eq r) {
#    $seq =~tr/ATCG/TAGC/;  
#    $seq = reverse $seq;
#  }
#  print ">$count\n";
#  print "$seq\n";
#  $count++
#}' > mytemp1 & 
# translate 
#java -jar ~/programs/bin/macse_v2.03.jar -prog translateNT2AA -gc_def 11 -seq mytemp1 -out_AA mytemp2 &
#tr -d '>'<mytemp2|paste - - > x.import_4.fa

# the above is slow, slow,slow - I think it's the java bit mostly

sqlite3 langdale.db 'SELECT start,stop,direction,sequence 
FROM genes_in_contigs 
LEFT JOIN contig_sequences 
WHERE genes_in_contigs.contig=contig_sequences.contig' | perl -e '
my $count=1;
while(<>)
{
  chomp;
  my @x = split /\|/,$_;
  my $seq =substr $x[3],$x[0],($x[1]-$x[0]+1);
  if ($x[2] eq r) {
    $seq =~tr/ATCG/TAGC/;  
    $seq = reverse $seq;
  }
  print ">$count\n";
  print "$seq\n";
  $count++
}' | ~/pipelines/common/scripts/translate.pl | tr -d '>'|paste - - > x.import_4.fa



sqlite3 -separator 'Control-v <TAB>' langdale.db  ".import x.import_4.fa gene_amino_acid_sequences"
sqlite3 langdale.db "select * from gene_amino_acid_sequences limit 10"
## update the self table to tell it we have gene calls
sqlite3 langdale.db "update self set value = 1 where key = 'genes_are_called'"

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

# profile bam files 
anvi-profile -i bam_file --min-contig-length 2000 --output-dir ./profiles --sample-name sample -c langdale.db

# Taxonomy
#sqlite3 langdale.db 'select contig,sequence,length(sequence) from contig_sequences where length(sequence)>1500 order by length(sequence) desc;'|
sqlite3 langdale.db 'SELECT start,stop,direction,sequence,contig_sequences.contig 
FROM genes_in_contigs 
LEFT JOIN contig_sequences 
WHERE genes_in_contigs.contig=contig_sequences.contig' | perl -e '
my $count=1;
while(<>)
{
  chomp;
  my @x = split /\|/,$_;
  my $seq =substr $x[3],$x[0],($x[1]-$x[0]+1);
  if ($x[2] eq r) {
    $seq =~tr/ATCG/TAGC/;  
    $seq = reverse $seq;
  }
  print ">$x[4]\n";
  print "$seq\n";
  $count++
}' > fortaxa.fa



awk -F"|" '{print ">"$1"\n"$2}' > fortaxa.fa
### kaiju
kaiju-makedb -s nr_euk
