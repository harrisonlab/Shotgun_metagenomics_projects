
# for some reason this R script didn't work correctly very odd...
# it was due to a problem in reordeing the rows - this is now fixed, but below was used for Chestnuts
for F in *.fasta; do
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $F|sed -e '1d' > ${F}.2
done

for F in *.2; do
  G=$(sed 's/\..*//' <<<$F) 
  grep ">.*$G" -A 1 --no-group-separator $F >${F}.3; 
done

for F in *.2; do
  G=$(sed 's/\..*//' <<<$F) 
  awk -F" " -v G=$G '($1~/^>/)&&($2!~G){line=$0;OUTF=$2".fasta.2.3";getline;print line >> OUTF;print >> OUTF}' $F
done

rename 's/\..*/.fasta/' *.3
