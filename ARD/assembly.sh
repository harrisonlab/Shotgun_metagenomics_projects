# metaspades
for FR in $PROJECT_FOLDER/data/corrected/*_1.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p metaspades \
  $PROJECT_FOLDER/data/assembled \
  $FR \
  $RR  \
  $PREFIX \
  -k 21,33,55,77
done

# megahit (using pre-merged reads) - this is dross, dropping 
for FR in $PROJECT_FOLDER/data/merged/*_1.unmerged.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit \
  $PROJECT_FOLDER/data/assembled/megahit/megahit_merged \
  $PREFIX \
  -r $MR,$FR,$RR \
  --k-min=27 --k-step 10 --k-max 127
done

# megahit (using unmerged reads)
for FR in $PROJECT_FOLDER/data/corrected/*_1.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit \
  $PROJECT_FOLDER/data/assembled/megahit \
  $PREFIX \
 -1 $FR -2 $RR -r $MR \
 --k-min=27 --k-step 10 --k-max 127
done

# assembly with spades
for FR in $PROJECT_FOLDER/data/corrected/*_1.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR) #t his line is specific to the file naming convention
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p metaspades \
  $PROJECT_FOLDER/data/assembled/spades \
  $FR \
  $RR  \
  $PREFIX \
  -k 21,33,55,77
done

# assemblies are compressed, bbmap requires them to be uncompressed for alignment
cd $PROJECT_FOLDER/data/assembly
Find . -maxdepth=3 type=f -name *.gz|xargs -I% pigz -dp 4 %

# align reads to assembly - spades
for FR in $PROJECT_FOLDER/data/cleaned/*_1.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/data/assembled/spades/$PREFIX/scaffolds.fasta \
  $PROJECT_FOLDER/data/assembly_checks/spades \
  $FR \
  $RR \
  nodisk=t \
  kfilter=22 \
  subfilter=15 \
  maxindel=80 \
  unpigz=t \
  touppercase=t \
  t=8
done

# align reads to assembly - megahit
for FR in $PROJECT_FOLDER/data/cleaned/*_1.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.contigs.fa \
  $PROJECT_FOLDER/data/assembly_checks/megahit \
  $FR \
  $RR \
  nodisk=t \
  kfilter=22 \
  subfilter=15 \
  maxindel=80 \
  unpigz=t \
  touppercase=t \
  t=8
done

# align unaligned megahit with spades
for FR in $PROJECT_FOLDER/data/assembly_checks/megahit/*_1.cleaned.fq.gz.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR) #this line is specific to the file naming convention
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p metaspades \
  $PROJECT_FOLDER/data/assembled/spades_of_megahit \
  $FR \
  $RR  \
  $PREFIX \
  -k 21,33,55,77
done

# align unaligned spades with megahit
for FR in $PROJECT_FOLDER/data/assembly_checks/spades/*_1.cleaned.fq.gz.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit \
  $PROJECT_FOLDER/data/assembled/megahit_of_spades \
  $PREFIX \
 -1 $FR -2 $RR -r $MR \
 --k-min=27 --k-step 10 --k-max 127
done
