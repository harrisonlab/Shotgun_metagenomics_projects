# To speed up assembly cutting data into site specific and assembling with megahit (accepts multiple fq input)

# Bigwood
f=$(ls -m $PROJECT_FOLDER/data/corrected/B*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
1 blacklace[01][06789].blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
BIGWOOD \
-24 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Attingham
f=$(ls -m $PROJECT_FOLDER/data/corrected/A*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
1 blacklace[01][06789].blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
ATTINGHAM \
-24 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Chestnuts
f=$(ls -m $PROJECT_FOLDER/data/corrected/C*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
1 blacklace[01][06789].blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
CHESTNUTS \
-24 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Gt_Monk
f=$(ls -m $PROJECT_FOLDER/data/corrected/G*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
1 blacklace[01][06789].blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
GT_MONK \
-24 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Langdale
f=$(ls -m $PROJECT_FOLDER/data/corrected/L*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
1 blacklace[01][06789].blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
LANGDALE \
-24 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Speculation
f=$(ls -m $PROJECT_FOLDER/data/corrected/S*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
1 blacklace[01][06789].blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
SPECULATION \
-24 $f -2 $r --k-min=27 --k-step 10 --k-max 127

# Winding
f=$(ls -m $PROJECT_FOLDER/data/corrected/W*_1.fq.gz|tr -d ' '|tr -d '\n')
r=$(sed 's/1\.fq\.gz/2.fq.gz/g' <<<$f)
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit2 \
1 blacklace[01][06789].blacklace \
$PROJECT_FOLDER/data/assembled/megahit \
WINDING \
-24 $f -2 $r --k-min=27 --k-step 10 --k-max 127
