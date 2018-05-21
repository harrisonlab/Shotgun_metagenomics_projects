# metaspades
for FR in $PROJECT_FOLDER/data/corrected/*_1.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'M[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c assemble -p metaspades \
  $PROJECT_FOLDER/data/assembled \
  $FR \
  $RR  \
  $PREFIX \
  -k 21,33,55,77
done

# megahit (using pre-merged reads)
for FR in $PROJECT_FOLDER/data/merged/*_1.unmerged.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'M[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit \
  $PROJECT_FOLDER/data/assembled/megahit_merged \
  $PREFIX \
  -r $MR,$FR,$RR \
  --k-min=27 --k-step 10 --k-max 127
done

# megahit (using unmerged reads)
for FR in $PROJECT_FOLDER/data/corrected/*_1.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'M[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit \
  $PROJECT_FOLDER/data/assembled/megahit_unmerged \
  $PREFIX \
 -1 $FR -2 $RR -r $MR \
 --k-min=27 --k-step 10 --k-max 127 --bubble-level 0
done

# align with bbmap 
# assemblies are compressed, bbmap requires them to be uncompressed for alignment
cd $PROJECT_FOLDER/data/assembly
Find . -maxdepth=3 type=f -name *.gz|xargs -I% pigz -dp 4 %

# align reads to assembly - spades
for FR in $PROJECT_FOLDER/data/cleaned/*_1.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'M[0-9]+.' <<<$FR)
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

# assemble assemblies
megahit -o output -t 20 --kmin-1pass --out-prefix all_spades --k-min=27 --k-step 20 --k-max 127 -r \
$PROJECT_FOLDER/data/assembled/spades/M1S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M2S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M4S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M4H/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M11S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M14S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M14H/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M1S/M1S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M2S/M2S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M4S/M4S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M4H/M4H.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M11S/M11S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M14S/M14S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M14H/M14H.contigs.fa.gz &

# assemble assemblies - useful samples
megahit -o output -t 20 --kmin-1pass --out-prefix all_spades --k-min=27 --k-step 20 --k-max 127 -r \
$PROJECT_FOLDER/data/assembled/spades/M1H/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M1S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M2H/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/M2S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M1H/M1H.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M1S/M1S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M2H/M2H.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/M2S/M2S.contigs.fa.gz 
