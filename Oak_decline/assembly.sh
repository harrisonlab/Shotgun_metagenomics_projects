# To speed up assembly cutting data into site specific

# megahit 
for FR in $PROJECT_FOLDER/data/corrected/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
  $PROJECT_FOLDER/data/assembled/megahit \
  $PREFIX \
 -1 $FR -2 $RR -r $MR \
 --k-min=27 --k-step 10 --k-max 127
done
