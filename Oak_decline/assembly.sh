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
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Attingham
f=$(ls -m $PROJECT_FOLDER/data/corrected/A*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
12 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
ATTINGHAM \
-m 0.5 -1 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Chestnuts
f=$(ls -m $PROJECT_FOLDER/data/corrected/C*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
12 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
CHESTNUTS \
-m 0.5 -1 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Gt_Monk
f=$(ls -m $PROJECT_FOLDER/data/corrected/G*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
12 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
GT_MONK \
-m 0.5 -1 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Langdale
f=$(ls -m $PROJECT_FOLDER/data/corrected/L*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
12 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
LANGDALE \
-m 0.5 -1 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Speculation
f=$(ls -m $PROJECT_FOLDER/data/corrected/S*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
12 blacklace11.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
SPECULATION \
-m 12 -1 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Winding
f=$(ls -m $PROJECT_FOLDER/data/corrected/W*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
10 blacklace01.blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
WINDING \
-1 $f -2 $r --k-min=27 --k-step 10 --k-max 127
