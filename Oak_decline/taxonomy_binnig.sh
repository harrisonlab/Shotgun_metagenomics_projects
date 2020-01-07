# binning using metabat
# pipeline needs sorted bam files

# sort bam files
for f in *.bam; do
 PREFIX=$(echo $f|sed -e 's/\..*//')
 sbatch --mem-per-cpu 2000M -c 10 \
 ~/pipelines/metagenomics/scripts/slurm/sub_bam_sort.sh \
 10 /data/data2/scratch2/deakig/Oak/sorted $PREFIX /data/data2/scratch2/deakig/Oak/to_sort/$f
done

# run metabat
# runMetaBat.sh -i assembly.fa.gz --unbinned -o assembly_name -m 1500 -x 0 --minCVSum 0.5 bam_files 

runMetaBat.sh ATTINGHAM.fa --unbinned -o assembly_name -m 1500 -x 0 --minCVSum 0.5 bam_files 
