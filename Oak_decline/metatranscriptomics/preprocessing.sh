# add project folders
PROJECT_FOLDER=~/projects/myproject
ln -s ~/pipelines/metatranscriptomics $PROJECT_FOLDER/metatranscriptomics_pipeline

mkdir $PROJECT_FOLDER/data
mkdir $PROJECT_FOLDER/data/cluster
mkdir $PROJECT_FOLDER/data/fastq
mkdir $PROJECT_FOLDER/data/trimmed
mkdir $PROJECT_FOLDER/data/filtered
mkdir $PROJECT_FOLDER/data/normalised
mkdir $PROJECT_FOLDER/data/cleaned
mkdir $PROJECT_FOLDER/data/corrected
mkdir $PROJECT_FOLDER/data/merged

# adapter/phix/rrna filtering
for FR in $PROJECT_FOLDER/data/fastq/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c MEGAFILT \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/adapters/truseq.fa \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/contaminants/ribokmers.fa.gz \
  $PROJECT_FOLDER/data/filtered \
  $FR \
  $RR
done  


# Human contaminant removal (BBMap)
for FR in $PROJECT_FOLDER/data/filtered/*_1.fq.gz.filtered.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/contaminants/bbmap_human \
  $PROJECT_FOLDER/data/cleaned \
  $FR \
  $RR \
  minid=0.95 \
  maxindel=3 \
  bwr=0.16 \
  bw=12 \
  quickmatch \
  fast \
  minhits=2 \
  t=8
done

# error correct and normalise
for FR in $PROJECT_FOLDER/data/cleaned/*_1.fq.gz.filtered.fq.gz.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c normalise -p bbnorm \
  $PROJECT_FOLDER/data/corrected \
  $FR \
  $RR  \
  target=100 \
  min=2 \
  ecc=t \
  passes=1 \
  bits=16 prefilter
done

# rename files (o.k this could have been implemeneted in each of the above scripts - maybe at some time)
find $PROJECT_FOLDER/data -type f -n *.fq.gz|rename 's/(.*_[12]).*(\.[a-zA-Z]+\.fq\.gz$)/$1$2/'
