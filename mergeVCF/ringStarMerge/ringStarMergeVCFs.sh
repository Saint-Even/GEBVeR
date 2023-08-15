#!/bin/bash
##           (\.---./)
##            /.-.-.\
##           /| 0_0 |\
##          |_`-(v)-'_|
##          \`-._._.-'/      .-.
##    -~-(((.`-\_|_/-'.)))-~' <_
##           `.     .'
##             `._.'
##
##    -----~--~---~~~----~-`.-;~
##            GEBVeR

#requires an input dir with uncompressed vcf files

#activae an environment with:
#parallel
#rtgtools
#bcftools

echo "Begin: common marker intersect"
#clear data and output
for d in data output
do
	rm -r ./${d}
	mkdir ${d}
	(cd ${d} || exit
	touch .gitkeep)
done

echo "copy inputs to data"
(cd input || exit
files=$(ls -1 *.vcf)
parallel \
	echo {}';' \
	cp  {} ../data \
::: ${files}
)

echo "compress and index files"
(cd data || exit
files=$(ls -1 *.vcf)
parallel \
	echo {}';' \
	rtg bgzip {}';' \
	bcftools index -cf --threads 4 {}.gz \
::: ${files}
)

echo "ring intersect"
(cd data || exit
IFS=$'\n'
files=($(find . -name "*.gz"))
unset IFS

len=${#files[@]}
last=$((len - 1))
for i in $(seq 0 ${last})
do
	#only once assign arg from array
	if [ ${i} = 0 ]
	then
		a=${files[0]}
	#later assign arg from previous result
	else
		a="./intersect.vcf.gz"
	fi

	#do nothing on last pass, no next value to intersect
	if [ ${i} -lt ${last} ]
	then
		next=$((i+1))
		b=${files[${next}]}
		echo a:${a}
		echo b:${b}
		bcftools isec \
			-c snps \
			-O z \
			-p isec \
			--threads 4 \
			${a} \
			${b}
		rm ./intersect.vcf.gz
		mv isec/0003.vcf.gz ./intersect.vcf.gz
		rm -r isec
		bcftools index -cf --threads 4 intersect.vcf.gz
	fi
done
)


echo "star intersect"
(cd data || exit
mv intersect.vcf.gz intersect.vcf.hidden
files=$(ls -1 *.gz)
mv intersect.vcf.hidden intersect.vcf.gz

parallel \
	echo {}';' \
	bcftools isec \
		-c snps \
		-O z \
		-p {/.} \
		--threads 4 \
		intersect.vcf.gz \
		{}';' \
	cd {/.}';' \
	mv 0003.vcf.gz {}.common';' \
	bcftools index -cf --threads 4 {}.common';' \
	mv {}.common* ..';' \
	cd .. ';' \
	rm -r {/.} \
::: ${files}
)

echo "merge"
(cd data || exit
files=$(ls *.vcf.gz.common)
echo ${files}

name=""
clean=$(echo ${files} | sed s:.vcf.gz.common::g)
for n in ${clean}
do
	name+="${n}_"
done
name+="common.vcf"

bcftools merge \
	--force-samples \
	--merge none \
	-o ${name} \
	-O v \
	--threads 4 \
	${files}

#copy result to output
mv ${name} ../output
)
echo "Complete"
