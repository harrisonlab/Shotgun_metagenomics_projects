f=$(ls -m $PROJECT_FOLDER/data/corrected/B*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.cleaned/2.cleaned/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
# blacklace[01][06789].blacklace \
10 blacklace01.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
BIGWOOD \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 77


# assemble spades MT assemblies + MT and MG corrected reads - everything assembly
for F in $PROJECT_FOLDER/assembly/mt_mg_data/*.gz; do
  R=$(echo $R,$F)
done  
R=$(sed '1s/^,//' <<<$R)
qsub -l h=blacklace11 qs.sh -pe smp 24 $PROJECT_FOLDER/assembly/everything spades \
--k-min=27 --k-step 20 --k-max 127 \
-r $R

# align MG + MT cleaned reads to everything assembly

for F in $PROJECT_FOLDER/assembly/everything/checks/*.gz; do
  U=$(echo $U,$F)
done  
U=$(sed '1s/^,//' <<<$U)
qsub -l h=blacklace11  -pe smp 24  qs.sh $PROJECT_FOLDER/assembly/everything/unaligned mega_un \
--k-min=27 --k-step 20 --k-max 127 \
-r $U
