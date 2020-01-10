# binning using metabat
# pipeline needs sorted bam files

# sort bam files
for f in *.bam; do
 PREFIX=$(echo $f|sed -e 's/\..*//')
 sbatch --mem-per-cpu 2000M -c 10 \
 ~/pipelines/metagenomics/scripts/slurm/sub_bam_sort.sh \
 10 /data/data2/scratch2/deakig/Oak/sorted $PREFIX /data/data2/scratch2/deakig/Oak/to_sort/$f
done

# get list of bam files for each assembly
A=$(for f in ../sorted/A*; do echo $f; done|tr  '\n' ' ')
G=$(for f in ../sorted/G*; do echo $f; done|tr  '\n' ' ')
L=$(for f in ../sorted/L*; do echo $f; done|tr  '\n' ' ')
W=$(for f in ../sorted/W*; do echo $f; done|tr  '\n' ' ')


# run metabat
# runMetaBat.sh -i assembly.fa.gz --unbinned -o assembly_name -m 1500 -x 0 --minCVSum 0.5 bam_files 

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/ATTINGHAM/ATTINGHAM.contigs.fa $A &

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/GTMONK/GTMONK.contigs.fa $G &

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/LANGDALE/LANGDALE.contigs.fa $L &

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/WINDING/WINDING.contigs.fa $W &

jgi_summarize_bam_contig_depths --outputDepth GTMONKdepth.txt ../sorted/G*
metabat2 -i ~/projects/Oak_decline/metagenomics/data/assembled/megahit/GTMONK/GTMONK.contigs.fa -a GTMONKdepth.txt -o bins_dir/bin 
