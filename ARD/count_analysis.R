```shell
awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" N11S/quant.sf
```

