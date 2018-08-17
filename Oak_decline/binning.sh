PREFIX=BIGWOOD # and etc.
P1=${PREFIX:0:1}

# functional binning with HirBin
# I've had to hack some of the HirBin scripts (specifically clusterbinstosubbins.py) as it doesn't work in current format
# also it uses usearch for clustering, while this is good the free 32bit version will almost certainly run out of memory for any sort of 
# soil metagenome assembly. I will probably change this to use vsearch which will give output the same as usearch

# HirBin does an hmm alignment of an assembly to a protein domain database 
# this will take a looong time unless the assembly and preferbly the hmm database is divided into chunks

# it has three/four steps

# annotate uses functionalAnnotaion.py, but splits input file into 20,000 droid chunks for running on cluster (25 concurrent jobs)
#functionalAnnotation.py -m METADATA_FILE -db DATABASE_FILE -e EVALUE_CUTOFF -n N -p MAX_ACCEPTABLE_OVERLAP

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

grep -v "#" BIGWOOD.hmmout|awk -F" " '{print $4,$1,$20,$21,"+",$7}' OFS="\t" > BIGWOOD.hmm.cut
awk -F"\t" '{print $1}' BIGWOOD.hmm.cut|sort|uniq > BIGWOOD.domains

# mapping
# mapping is not implemented very well in HirBin, will do this seperately with bbmap
# align reads to assembly - will need to index first
bbmap.sh ref=$PREFIX.contigs.fa.gz usemodulo=T #k=11


for FR in $PROJECT_FOLDER/data/fastq/$P1*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  16 blacklace[01][0-9].blacklace \
  $PROJECT_FOLDER/data/assembled/aligned/megahit \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.contigs.fa.gz \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  path=$PROJECT_FOLDER/data/assembled/megahit/$PREFIX/ 
  usemodulo=T 
done

# bedtools code is inefficient at getting over-lapping counts (if min overlap is set to 1)
# I've written something in perl which is way less memory hungry and takes about a millionth of the time to run
# output is not a cov file but just counts per domain - not certain the sub-binning is worth while (could modify bam_count to return a cov/tab file to implement this step)
# takes about ten minutes on a single core to run, could easily get it to produce a cov file
# bam_scaffold_count.pl will output a cov file rather than counts per domain
samtools view bam_file|~/pipelines/metagenomics/scripts/bam_count.pl $PREFIX.gff > bam_file.txt

for BAM in $PROJECT_FOLDER/data/assembled/aligned/megahit/$P1*.bam; do
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c coverage -p bam_count \
  blacklace[01][0-9].blacklace \
  $BAM \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/assembled/counts/megahit
done


# Sub binning - if required
# I've hacked around with a few of the HirBin settings
# ParsePFamTGRFAM.py accepts a chopped up version of the hhm output and the domains in a seperate file
# This was done for speed reasons - could probably use something similar to bam_count, if I can get it to work under python
# also the hmm file is better if cut to include only the necessary fields (prevents having to do various checking) 
# will require a cov file from bam_scaffold_count.pl
awk -F"\t" '{sub("ID=","|",$(NF-1));OUT=$1$(NF-1)":"$4":"$5":"$7;print OUT,$NF}' OFS="\t" $PREFIX.cov > $PREFIX.tab
grep "#" -v $PREFIX.hmmout|awk -F" " '{print $4,$1,$20,$21,"+",$7}' OFS="\t" > $PREFIX.cut.hmm
cut -f9 $PREFIX.gff|sort|uniq|sed 's/ID=//'|tail -n +2 > $PREFIX.domains # the tail bit gets rid of the first line of output
# then create the required metadata file
echo -e \
"Name\tGroup\tReference\tAnnotation\tCounts\Domain\n"\
"$PREFIX\tSTATUS\t$PREFIX.pep\t$PREFIX.hmm.cut\t$PREFIX.tab\t$PREFIX.domains" > metadata.txt

clusterBinsToSubbins.py -m metadata.txt -id 0.7 --onlyClustering  # this will create the sub bins
clusterBinsToSubbins.py -m metadata.txt -id 0.7 --onlyParsing # this will make count files for $PREFIX.tab to the bins and sub bins 
clusterBinsToSubbins.py -m metadata.txt -id 0.95 --reClustering # recluster at a different identity plus parsing
clusterBinsToSubbins.py -m metadata.txt -id 0.95 --reClustering --onlyClustering # as above but without the parsing

