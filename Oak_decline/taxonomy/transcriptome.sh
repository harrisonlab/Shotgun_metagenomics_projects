# Make salmon index from each of the combined bin files
salmon index -t ~/projects/Oak_decline/metagenomics/data/taxonomy/ATTINGHAM_BINS/ATTINGHAM.bins.fa \
-i ~/projects/Oak_decline/metagenomics/data/taxonomy/ATTINGHAM_BINS/salmon

salmon index -t ~/projects/Oak_decline/metagenomics/data/taxonomy/GTMONK_BINS/GTMONK.bins.fa \
-i ~/projects/Oak_decline/metagenomics/data/taxonomy/GTMONK_BINS/salmon

salmon index -t ~/projects/Oak_decline/metagenomics/data/taxonomy/LANGDALE_BINS/LANGDALE.bins.fa \
-i ~/projects/Oak_decline/metagenomics/data/taxonomy/LANGDALE_BINS/salmon

salmon index -t ~/projects/Oak_decline/metagenomics/data/taxonomy/WINDING_BINS/WINDING.bins.fa \
-i ~/projects/Oak_decline/metagenomics/data/taxonomy/WINDING_BINS/salmon
