# functional binning with HirBin
# I've had to hack some of the HirBin scripts (specifically clusterbinstosubbins.py) as it doesn't work in current format

# HirBin does an hmm alignment of an assembly to a protein domain database 
# this will take a looong time unless the assembly and preferbly the hmm database is divided into chunks

# it has three/four steps

# annotate uses functionalAnnotaion.py, but splits input file into 20,000 droid chunks for running on cluster (25 concurrent jobs)
#functionalAnnotation.py -m METADATA_FILE -db DATABASE_FILE -e EVALUE_CUTOFF -n N -p MAX_ACCEPTABLE_OVERLAP
PREFIX=BIGWOOD
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
  1 $PROJECT_FOLDER/data/assembled/megahit/$PREFIX \
 $PREFIX.contigs.fa \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03

# concatenate annotate output
find -type f -name X.gff|head -n1|xargs -I% head -n1 % >$PREFIX.gff
find -type f -name X.gff|xargs -I% grep -v "##" % >>$PREFIX.gff
find -type f -name X.pep|xargs -I% cat % >$PREFIX.pep
find -type f -name X.hmmout|head -n1|xargs -I% head -n3 % >$PREFIX.hmmout   
find -type f -name X.hmmout|xargs -I% grep -v "#" % >>$PREFIX.hmmout
find -type f -name X.hmmout|head -n1|xargs -I% tail -n10 % >>$PREFIX.hmmout
 
# mapping
# mapping is not implemented very well in HirBin, will do this seperately with bbmap + HirBin tools to get BAM in correct format
# align reads to assembly - will need to index first
bbmap ref=$PREFIX.contigs.fa.gz usemodulo=T k=11

P1=${PREFIX:0:1}
for FR in $PROJECT_FOLDER/data/fastq/$P1*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  16 blacklace[01][0-9].blacklace \
  $PROJECT_FOLDER/data/assembled/aligned/megahit \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.contigs.fa \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  path=$PROJECT_FOLDER/data/assembled/megahit/$PREFIX/ usemodulo=T k=11
done
 
