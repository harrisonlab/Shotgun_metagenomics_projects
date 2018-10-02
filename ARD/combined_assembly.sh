# assemble all metatranscriptome samples
f=$(ls -m $PROJECT_FOLDER/data/cleaned/M*_1.cleaned.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.cleaned/2.cleaned/g' <<<$f)
megahit -o OUTPUT -t 24 --kmin-1pass --out-prefix ARD_TRANS -1 $f -2 $r --k-min=27 --k-step 10 --k-max 77

# assemble all metagenomics samples with transcriptome 
f=$(ls -m $PROJECT_FOLDER/data/cleaned/N*_1.cleaned.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.cleaned/2.cleaned/g' <<<$f)
megahit -o OUTPUT -t 24 --kmin-1pass --out-prefix ARD_COMB -1 $f -2 $r -r ARD_TRANS.contigd.fa --k-min=27 --k-step 10 --k-max 77

# map to combined (ARD_COMB) assembly with bbmap
bbmap.sh ref=ARD_COMB.contigs.fa.gz usemodulo=t 

for FR in $PROJECT_FOLDER/data/cleaned/N*_1.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  24 blacklace11.blacklace \
  $PROJECT_FOLDER/data/aligned/COMB \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/ARD_COMB.contigs.fa.gz \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  path=$PROJECT_FOLDER/data/assembled/ \
  usemodulo=t 
done

# Pseudo mapping of RNA samples
salmon index ...

for FR in $PROJECT_FOLDER/data/cleaned/M*_1.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  OUTDIR=$(echo $FR|awk -F"/" '{print $NF}'|sed 's/_.*//')
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p salmon \
  24 blacklace11.blacklace \
  /data2/scratch2/deakig/SALMON_COMB \
  $PROJECT_FOLDER/data/counts/$OUTDIR \
  $FR \
  $RR \
  --numBootstraps 1000 \
  --dumpEq \
  --seqBias \
  --gcBias
done  


##### OLD STUFF BELOW #####

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
