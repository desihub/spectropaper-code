#!/bin/bash

rm tmp
for r in `seq 20` ; do
    echo $r
    ./auto_calibration.py | grep RESULT | awk '{print $4,$5,$8,$9}' >> tmp

done

cat tmp | awk '{n +=1; sum1 += $1; sum2 += $2 ; sum3 += $3 ; sum4 += $4 }END{print sum1/n,sum2/n,sum3/n,sum4/n}'

