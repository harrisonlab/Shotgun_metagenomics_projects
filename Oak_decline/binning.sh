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
 
 
