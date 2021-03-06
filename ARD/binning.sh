# set variables
PROJECT_FOLDER=~/projects/ARD/metagenomics
PREFIX=COMB
P1=${PREFIX:0:1}

# identify pfam domains in assemblies
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
 1 25 $PROJECT_FOLDER/data/assembled \
 ARD_COMB.contigs.fa.gz \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03

# concatenate annotation output
find -type f -name X.gff|head -n1|xargs -I% head -n1 % >$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.gff
find -type f -name X.gff|xargs -I% grep -v "##" % >>$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.gff
find -type f -name X.pep|xargs -I% cat % >$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.pep
find -type f -name X.hmmout|xargs -I% grep -v "#" % | \
awk -F" " '($21~/^[0-9]+$/) && ($20~/^[0-9]+$/) {print $4,$1,$20,$21,$3,$7}' OFS="\t"| \
$PROJECT_FOLDER/metagenomics_pipeline/scripts/subbin_domain_extractor.pl \
> $PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.domains

# extract sub bins from proteins file with subbin_fasta_extractor.R - last two args control memory usage (and exection time)
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/subbin_fasta_extractor.R $PREFIX.domains $PREFIX.pep "${PREFIX}_clustering/forClustering" 100 T

# Clustering
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c cluster_super_fast \
  blacklace[01][0-9].blacklace 100 \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/forClustering \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/clust0.7 \
  0.7
  
# concatenate clustering output
cat $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/clust0.5/*.uc > $PROJECT_FOLDER/data/binning/$PREFIX/reduced.txt

# create bbmap assembly reference
bbmap.sh ref=ARD_COMB.contigs.fa.gz usemodulo=t 

# map to assembly (large assembly ace11 only - or 1 with less processors)
# REMEMBER: SET PREFIX (though it's not used)
for FR in $PROJECT_FOLDER/data/cleaned/*_1.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  24 blacklace11.blacklace \
  $PROJECT_FOLDER/data/assembled/aligned \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/FINAL_COMBINED.assembly.fa \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  usemodulo=t \
  path=$PROJECT_FOLDER/data/assembled
done

# count overlapping features
for BAM in $PROJECT_FOLDER/data/assembled/aligned/*.bam; do
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c coverage -p bam_count \
  blacklace[01][0-9].blacklace \
  $BAM \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/binning/$PREFIX \
  cov
done

# count bins ~ sub bin parsing will also output bin count table (this is still here as counting sub-bins is v. slow)
# Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/cov_count.R "." ".*\\.cov" "$PREFIX.countData"

# Sub binning - convert cov to tab
for F in $P1*.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{sub("ID=","",$(NF-1));OUT=$1"_"$(NF-1)"_"$4"_"$5;print OUT,$(NF-1),$4,$5,$NF}' OFS="\t" $F > $O
done 

# get bin and sub_bin counts (final argument is number of cores to use, can be memory hungry)
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/subbin_parser_v2.R reduced.txt . $PREFIX.countData 8
