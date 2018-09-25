
awk -F"\t" '{gsub(/ID=/,"");print $1"_"++count[$1],$4,$5,$7,$9}' OFS="\t" COMB.gff > ../../assembled/COMB_trasncriptome.map

