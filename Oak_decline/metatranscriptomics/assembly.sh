# combined assembly of metgenomics assembly + metatranscriptomic read

# Attingham
cd $TMP
cp $PROJECT_FOLDER/../metagenomics/data/assembled/megahit/ATTINGHAM/ATTINGHAM.contigs.fa.gz .
pigz -d ATTINGHAM.contigs.fa.gz
f=$(ls -m $PROJECT_FOLDER/data/corrected/A*_1.corrected.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.corrected/2.corrected/g' <<<$f)
megahit -o OUTPUT -t 21 --kmin-1pass --out-prefix ATTINGHAM_COMB -1 $f -2 $r -r ATTINGHAM.contigs.fa --k-min=27 --k-step 10 --k-max 77 --tmp-dir /data/scratch/deakig/tmp/ATTINGHAM

# Langdale
cd $TMP
cp $PROJECT_FOLDER/../metagenomics/data/assembled/megahit/LANGDALE/LANGDALE.contigs.fa.gz .
pigz -d LANGDALE.contigs.fa.gz
f=$(ls -m $PROJECT_FOLDER/data/corrected/L*_1.corrected.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.corrected/2.corrected/g' <<<$f)
megahit -o OUTPUT -t 24 --kmin-1pass --out-prefix LANGDALE_COMB -1 $f -2 $r -r LANGDALE.contigs.fa --k-min=27 --k-step 10 --k-max 77 --tmp-dir /data/scratch/deakig/tmp/LANGDALE

# index assemblies
bbmap.sh ref=ATTINGHAM_COMB.contigs.fa.gz usemodulo=t
bbmap.sh ref=LANGDALE_COMB.contigs.fa.gz usemodulo=t

# map reads
PREFIX=ATTINGHAM
for FR in $PROJECT_FOLDER/data/fastq/A*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  16 blacklace[01][0-9].blacklace \
  $PROJECT_FOLDER/data/assembly/aligned \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/$PREFIX/${PREFIX}_COMB.contigs.fa.gz \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  usemodulo=t \
  path=$PROJECT_FOLDER/data/assembled/$PREFIX
done

PREFIX=LANGDALE
for FR in $PROJECT_FOLDER/data/fastq/L*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  16 blacklace[01][0-9].blacklace \
  $PROJECT_FOLDER/data/assembly/aligned \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/$PREFIX/${PREFIX}_COMB.contigs.fa.gz \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  usemodulo=t \
  path=$PROJECT_FOLDER/data/assembled/$PREFIX
done


# single assemblies
# Attingham - produces tmp files too large for tmp directory on nodes
f=$(ls -m $PROJECT_FOLDER/data/corrected/A*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
24 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
ATTINGHAM \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77 --tmp-dir /data/scratch/deakig/tmp/ATTINGHAM

# Langdale - produces tmp files too large for tmp directory on nodes
f=$(ls -m $PROJECT_FOLDER/data/corrected/L*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
24 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
LANGDALE \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77 --tmp-dir /data/scratch/deakig/tmp/LANGDALE
