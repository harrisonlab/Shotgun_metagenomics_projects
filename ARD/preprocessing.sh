# add project folders
mkdir ~/projects/ARD/data/metagenomics
$PROJECT_FOLDER=~/projects/ARD/metagenomics
ln -s ~/pipelines/metagenomics $PROJECT_FOLDER/metagenomics_pipeline

mkdir $PROJECT_FOLDER/data
mkdir $PROJECT_FOLDER/data/fastq
mkdir $PROJECT_FOLDER/data/trimmed
mkdir $PROJECT_FOLDER/data/filtered
mkdir $PROJECT_FOLDER/data/normalised
mkdir $PROJECT_FOLDER/data/merged

#Adapter trimming (trimmomatic)
for FR in $PROJECT_FOLDER/data/fastq/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c trim \
  $FR \
  $RR \
  $PROJECT_FOLDER/data/trimmed \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/adapters/truseq.fa \
  4
done

# phix filter
for FR in $PROJECT_FOLDER/data/trimmed/*_1.fq.gz.trimmed.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbduk \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/data/filtered \
  $FR \
  $RR  \
  k=31 \
  hdist=1
done

# human contaminant filter
for FR in $PROJECT_FOLDER/data/filtered/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/bbmap_human \
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


# normalition and error correction
for FR in $PROJECT_FOLDER/data/cleaned/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c normalise -p bbnorm \
  $PROJECT_FOLDER/data/corrected \
  $FR \
  $RR  \
  target=100 \
  min=5 \
  ecc=t \
  passes=1 \
  bits=16 prefilter
done
