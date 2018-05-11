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

# assemble assemblies
megahit -o output \
-t 12 \
--kmin-1pass \
--out-prefix megahit_of_spades \
--k-min=27 --k-step 20 --k-max 127 \
-r $PROJECT_FOLDER/data/assembled/megahit_of_spades/N11S/N11S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/N14H/N14H.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/N14S/N14S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/N1H/N1H.contigs.fa.gz,\
$PROJECT_FOLDERs/data/assembled/megahit_of_spades/N1S/N1S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/N2H/N2H.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/N2S/N2S.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/N3H/N3H.contigs.fa.gz,\
$PROJECT_FOLDER/data/assembled/megahit_of_spades/N4H/N4H.contigs.fa

megahit  -o output \
-t 12 \
--kmin-1pass \
--out-prefix megahit \
--k-min=27 --k-step 20 --k-max 127 \
-r $PROJECT_FOLDER/data/assembled/megahit/N11S/N11S.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N14H/N14H.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N14S/N14S.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N1H/N1H.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N1S/N1S.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N2H/N2H.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N2S/N2S.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N3H/N3H.contigs.fa,\
$PROJECT_FOLDER/data/assembled/megahit/N4H/N4H.contigs.fa

megahit  -o output \
-t 12 \
--kmin-1pass \
--out-prefix spades \
--k-min=27 --k-step 20 --k-max 127 \
-r $PROJECT_FOLDER/data/assembled/spades/N11S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/N14H/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/N14S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/N1H/scaffolds.fasta,\
/$PROJECT_FOLDER/data/assembled/spades/N1S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/N2H/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/N2S/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/N3H/scaffolds.fasta,\
$PROJECT_FOLDER/data/assembled/spades/N4H/scaffolds.fasta

# CAP 3 assembly - note cap3 can only assembleupto about 1G file size

Cap3 input.fa 
mv input.contigs.fa $PROJECT_FOLDER/data/assembled/Final.fa


