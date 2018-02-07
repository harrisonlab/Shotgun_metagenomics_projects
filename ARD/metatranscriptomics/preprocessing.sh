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

#### Adapter trimming
for FR in $PROJECT_FOLDER/data/fastq/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c trim \
  $FR \
  $RR \
  $PROJECT_FOLDER/data/trimmed \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/adapters/truseq.fa \
  4
done

#### Synthetic construct/contaminant filtering 
for FR in $PROJECT_FOLDER/data/trimmed/*_1.fq.gz.trimmed.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c filter -p bbduk \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/data/filtered \
  $FR \
  $RR \
  k=31 \
  hdist=1 \
  t=4
done


#### Human contaminant removal (BBMap)
for FR in $PROJECT_FOLDER/data/filtered/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz; do
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

# error correct and normalise (this was performed while there were write errors to /home hence output has gone to /Data drive
for FR in $PROJECT_FOLDER/data/cleaned/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c normalise -p bbnorm \
  /data/scratch/deakig/metatranscriptomics/corrected \
  $FR \
  $RR  \
  target=100 \
  min=2 \
  ecc=t \
  passes=1 \
  bits=16 prefilter
done

# merge pairs 
for FR in /data/scratch/deakig/metatranscriptomics/corrected/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz.cleaned.fq.gz.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c merge -p bbmerge-auto \
  /data/scratch/deakig/metatranscriptomics/merged \
  $FR \
  $RR  \
  rem k=62 \
  extend2=50 \
  t=12 \
  vstrict
done  

# rename files
rename 's/_1.*[^n]merged/\.merged/' *.gz
rename 's/_1.*unmerged/_1\.unmerged/' *.gz
rename 's/_2.*unmerged/_2\.unmerged/' *.gz
