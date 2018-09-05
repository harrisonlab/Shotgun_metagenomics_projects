# set variables
$PROJECT_FOLDER=~/projects/OAK_decline/metagenomics
PREFIX=BIGWOOD # and etc.
P1=${PREFIX:0:1}

# identify pfam domains in assembly
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
 1 $PROJECT_FOLDER/data/assembled/megahit/$PREFIX \
 $PREFIX.contigs.fa \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03

# concatenate annotation output
find -type f -name X.gff|head -n1|xargs -I% head -n1 % >$PREFIX.gff
find -type f -name X.gff|xargs -I% grep -v "##" % >>$PREFIX.gff
find -type f -name X.pep|xargs -I% cat % >$PREFIX.pep
find -type f -name X.hmmout|head -n1|xargs -I% head -n3 % >$PREFIX.hmmout   
find -type f -name X.hmmout|xargs -I% grep -v "#" % >>$PREFIX.hmmout
find -type f -name X.hmmout|head -n1|xargs -I% tail -n10 % >>$PREFIX.hmmout

grep -v "#" $PREFIX.hmmout|awk -F" " '($21~/^[0-9]+$/) && ($20~/^[0-9]+$/) {print $4,$1,$20,$21,$3,$7}' OFS="\t" > $PREFIX.hmm.cut
awk -F"\t" '{print $1}' $PREFIX.hmm.cut|sort|uniq > $PREFIX.domains # this is better method as some domains may not be present in gff due to filtering
# cut -f9 $PREFIX.gff|sort|uniq|sed 's/ID=//'|tail -n +2 > $PREFIX.2.domains # the tail bit gets rid of the first line of output an d the grep removes errors in the output

# mapping

bbmap.sh ref=$PREFIX.contigs.fa.gz usemodulo=t #k=11 

for FR in $PROJECT_FOLDER/data/fastq/$P1*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  16 blacklace[01][0-9].blacklace \
  $PROJECT_FOLDER/data/assembled/aligned/megahit \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.contigs.fa.gz \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  path=$PROJECT_FOLDER/data/assembled/megahit/$PREFIX/ 
  usemodulo=T 
done

# count overlapping features
for BAM in $PROJECT_FOLDER/data/assembled/aligned/megahit/$P1*.bam; do
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c coverage -p bam_count \
  blacklace[01][0-9].blacklace \
  $BAM \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/binning/$PREFIX \
  cov
done

Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/cov_count.R "." "$P1.*\\.cov" "$PREFIX.countData"

# Sub binning - convert cov to tab
for F in $P1*.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{sub("ID=","|",$(NF-1));OUT=$1$(NF-1)":"$4":"$5":"$7;print OUT,$NF}' OFS="\t" $F > $O
done 

for F in $P1*.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{sub("ID=","",$(NF-1));OUT=$1"_"$(NF-1)"_"$4"_"$5;print OUT,$(NF-1),$4,$5,$NF}' OFS="\t" $F > $O
done 


# then create the required metadata file
echo -e \
"Name\tGroup\tReference\tAnnotation\tCounts\tDomain\n"\
"$PREFIX\tSTATUS\t$PREFIX.pep\t$PREFIX.hmm.cut\tEMPTY\t$PREFIX.domains" > metadata.txt

# The extractSequences module of clusterBinsToSubbins is inpractical for large pep files - the R script subbin_fasta_extractor.R should be used in it's palace

Rscript subbin_fasta_extractor.R $PREFIX.hmm.cut $PREFIX.pep $PREFIX_hirbin_output
# for some reason this R script didn't work correctly on the latest round. The below awk scripts will quickly fix it, but very odd...
# it was due to a problem in reordeing the rows - this is now fixed, but below was used for Chestnuts

for F in *.fasta; do
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $F|sed -e '1d' > ${F}.2
done

for F in *.2; do
  G=$(sed 's/\..*//' <<<$F) 
  grep ">.*$G" -A 1 --no-group-separator $F >${F}.3; 
done

for F in *.2; do
  G=$(sed 's/\..*//' <<<$F) 
  awk -F" " -v G=$G '($1~/^>/)&&($2!~G){line=$0;OUTF=$2".fasta.2.3";getline;print line >> OUTF;print >> OUTF}' $F
done

rename 's/\..*/.fasta/' *.3

# Clustering
clusterBinsToSubbins.py -m $PREFIX.metadata.txt -id 0.7 --reClustering --onlyClustering -f -o  $PREFIX_hirbin_output# clustering without sub bin extraction (no parsing)
awk -F"\t" '($1~/[HS]/){print $2, $9, $10}' *.uc| \
awk -F" " '{sub(/_[0-9]+$/,"",$2);sub(/_[0-9]+$/,"",$6);A=$2"_"$3"_"$4"_"$5;if($6~/\*/){B=A}else{B=$6"_"$7"_"$8"_"$9};print A,B}' OFS="\t" > reduced.txt

# parse 
Rscript subbin_parser.R reduced.txt *.tab $PREFIX.countData.sub_bins

