# set variables
PROJECT_FOLDER=~/projects/Oak_decline/metatranscriptomics
PREFIX=ATTINGHAM
# PREFIX=LANGDALE
P1=${PREFIX:0:1}

# identify pfam domains in assemblies
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
 1 25 $PROJECT_FOLDER/data/assembly/$PREFIX \
 ${PREFIX}_COMB.contigs.fa.gz \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03

# concatenate annotation output
find -type f -name X.gff|head -n1|xargs -I% head -n1 % >$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.gff
find -type f -name X.gff|xargs -I% grep -v "##" % >>$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.gff &
find -type f -name X.pep|xargs -I% cat % >$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.pep &
find -type f -name X.hmmout|xargs -I% grep -v "#" % |tee $PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.hmmout |
awk -F" " '($21~/^[0-9]+$/) && ($20~/^[0-9]+$/) {print $4,$1,$20,$21,$3,$7}' OFS="\t"| \
$PROJECT_FOLDER/metagenomics_pipeline/scripts/subbin_domain_extractor.pl \
> $PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.domains &

# extract domains from proteins file with subbin_fasta_extractor.R - last three args control memory usage (and exection time)
# 100 = no. of chunks to cut the protein file into (larger will reduce memory footprint)
# T = use parallel processing 
# 8 = no. cores for parallel processing
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/subbin_fasta_extractor.R $PREFIX.domains $PREFIX.pep "${PREFIX}_clustering/forClustering" 100 T 8 &

# Clustering
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c cluster_super_fast \
  blacklace[01][0-9].blacklace 100 \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/forClustering \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/clust0.7 \
  0.7
  
# concatenate clustering output
# NOTE FOR SELF - check output is in new format otherwise run the awk script instead
cat $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/clust0.7/*.uc > $PROJECT_FOLDER/data/binning/$PREFIX/reduced.txt

# mapping
bbmap.sh ref=$PREFIX.contigs.fa.gz usemodulo=t #k=11 


# count overlapping features
for BAM in $PROJECT_FOLDER/data/aligned/$P1*.bam; do
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c coverage -p bam_count \
  blacklace[01][0-9].blacklace \
  $BAM \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/binning/$PREFIX/map
done

# count bins (depricated - parsing will do this now)
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/cov_count.R "." "$P1.*\\.cov" "$PREFIX.countData"

# Sub binning - convert cov to tab
for F in *.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{sub("ID=","",$(NF-1));OUT=$1"_"$(NF-1)"_"$4"_"$5;print OUT,$(NF-1),$4,$5,$NF}' OFS="\t" $F > $O &
done 

# parsing
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/subbin_parser_v2.R \
 $PROJECT_FOLDER/data/binning/$PREFIX/reduced.txt \
 $PROJECT_FOLDER/data/binning/$PREFIX/map \
 $PROJECT_FOLDER/analysis/$PREFIX \
 10 &
