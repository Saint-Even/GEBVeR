#!/bin/bash

#creates sym links in target dir

sou="/home/holdens/holdensQNAP/LIBS/gbs_snp/cassetteGBS/synergySRG/2021/001"
tar="/home/holdens/holdensQNAP/PROC/GEBVeR/mergeVCF/starIntersect/input"
ddd="NS.*"
f="variantsImputed.vcf"

for d in  ${sou}/${ddd}
do
  echo "Enter: ${d}"
  (cd "${d}" || exit
  
  echo "Build name for: ${f}"
  n=$(basename ${d} | sed s:N.*2021:2021:g | sed s:_::g)
  n+=".vcf"
  echo ${n}

  echo "Create renamed links in ${tar}"
  ln -ifs ${d}/${f} ${tar}/${n}
  echo ""
  )
done

exit 0
