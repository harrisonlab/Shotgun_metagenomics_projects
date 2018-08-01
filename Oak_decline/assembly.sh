# To speed up assembly cutting data into site specific and assembling with megahit (accepts multiple fq input)
# but needs minimum of 40Gb memory (upped min kmer to 31, still required 34Gb with the smallest dataset - annoying)

# Bigwood
f=$(ls -m $PROJECT_FOLDER/data/corrected/B*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
# blacklace[01][06789].blacklace \
10 blacklace01.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
BIGWOOD \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77

# Attingham - produces tmp files too large for tmp directory on nodes
f=$(ls -m $PROJECT_FOLDER/data/corrected/A*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
24 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
ATTINGHAM \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77 --tmp-dir /data/scratch/deakig/tmp/ATTINGHAM

# Chestnuts
f=$(ls -m $PROJECT_FOLDER/data/corrected/C*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
24 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
CHESTNUTS \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77

# Gt_Monk
f=$(ls -m $PROJECT_FOLDER/data/corrected/G*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
24 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
GT_MONK \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77 

# Langdale - produces tmp files too large for tmp directory on nodes
f=$(ls -m $PROJECT_FOLDER/data/corrected/L*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
24 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
LANGDALE \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77 --tmp-dir /data/scratch/deakig/tmp/LANGDALE

# Speculation
f=$(ls -m $PROJECT_FOLDER/data/corrected/S*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
24 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
SPECULATION \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77

# Winding
f=$(ls -m $PROJECT_FOLDER/data/corrected/W*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
10 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
WINDING \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77


# align reads to assembly - BIGWOOD
for FR in $PROJECT_FOLDER/data/corrected/B*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=BIGWOOD
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.contigs.fa \
  $PROJECT_FOLDER/data/assembled/unassembled/megahit \
  $FR \
  $RR \
  nodisk=t \
  kfilter=22 \
  subfilter=15 \
  maxindel=80 \
  unpigz=t \
  touppercase=t \
  t=10
done

# align reads to assembly - WINDING
for FR in $PROJECT_FOLDER/data/corrected/W*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=WINDING
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.contigs.fa \
  $PROJECT_FOLDER/data/assembled/unassembled/megahit \
  $FR \
  $RR \
  nodisk=t \
  kfilter=22 \
  subfilter=15 \
  maxindel=80 \
  unpigz=t \
  touppercase=t \
  t=10
done


# assemble unaligned megahit with megahit
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

# assemble assemblies - megahit

