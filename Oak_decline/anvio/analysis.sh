

# use ORFfinder to find - bacterial ORFS

ORFfinder -in LANGDALE.fa -g 11 -s 1 -ml 180 -n true -outfmt 3 -out temp.out
grep ">" temp.out |
sed 's/.*_c/c/'|
awk -F":" '{count++;if ($3>$2){d="f"}else{d="r";t=$2;$2=$3;$3=t};p=1;if(($3-$2+1)%3==0){p=0};print count,$1,$2,$3,d,p,"program","v1.0"}' OFS="\t" > langdale.out
