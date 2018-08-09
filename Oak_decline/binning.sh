# functional binning with HirBin
# I've had to hack some of the HirBin scripts (specifically clusterbinstosubbins.py) as it doesn't work in current format

# HirBin does an hmm alignment of an assembly to a protein domain database 
# this will take a looong time unless the assembly and preferbly the hmm database is divided into chunks

# it has three/four steps

# annotate 
#functionalAnnotation.py -m METADATA_FILE -db DATABASE_FILE -e EVALUE_CUTOFF -n N -p MAX_ACCEPTABLE_OVERLAP
PREFIX=BIGWOOD
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
  4 $PROJECT_FOLDER/data/assembled/megahit/$PREFIX \
 $PREFIX.contigs.fa \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03
 
 # and etc.
 
 # mapping is not implemented very well in HirBin, will do this seperately with bbmap + HirBin tools to get BAM in correct format
 
 
