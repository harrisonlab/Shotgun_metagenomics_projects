# set variables
PROJECT_FOLDER=~/projects/Oak_decline/metatranscriptomics
PREFIX=ATTINGHAM
# PREFIX=LANGDALE
P1=${PREFIX:0:1}

# identify pfam domains in assemblies
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
 1 25 $PROJECT_FOLDER/data/assembly/$PREFIX \
 ${PREFIX}_COMB.contigs.fa.gz \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03
