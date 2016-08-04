#!/bin/bash

INdir=$1
OUTdir=$2

ENDmain1=1; ENDindex=2; ENDmain2=3 

for INmain1 in "${INdir}"/s_?_"${ENDmain1}"_????_qseq.txt.gz; do
  INbn=`basename "${INmain1}"` 
    INindex="${INdir}/${INbn:0:4}${ENDindex}${INbn:5}" 
          INmain2="${INdir}/${INbn:0:4}${ENDmain2}${INbn:5}" 
  OUTbn=`basename "${INmain1}" _qseq.txt.gz` 
    OUT="${OUTdir}/${OUTbn:0:4}${OUTbn:6}".txt
  echo @@@ "${INindex} ${INmain1} ${INmain2} -> ${OUT}...  `date`"
  paste <(zcat "${INindex}" | cut -f 9     | tr . N)  \
        <(zcat "${INmain1}" | cut -f 9     | tr . N)  \
        <(zcat "${INmain2}" | cut -f 9     | tr . N)  \
        <(zcat "${INindex}" | cut -f 10            )  \
        <(zcat "${INmain1}" | cut -f 10            )  \
        <(zcat "${INmain2}" | cut -f 10,11         )  \
        <(zcat "${INindex}" | cut -f 1,2,3,4,5,6 | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}') > "${OUT}"
done; echo @@@ "Done.  `date`"