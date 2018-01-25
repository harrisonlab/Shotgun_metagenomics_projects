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
for R1 in $PROJECT_FOLDER/data/fastq/*_1.fq.gz; do
  R2=$(echo $R1|sed 's/_1/_2/')
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c trim \
  $R1 \
  $R2 \
  $PROJECT_FOLDER/data/$RUN/trimmed \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/adapters/truseq.fa \
  4
done
