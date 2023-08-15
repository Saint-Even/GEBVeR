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
#one vcf file must be named ending with .vcf.pivot
#at least one other file ending with .vcf

#activate an environment with:
#parallel
#rtgtools
#bcftools

echo "Begin: Star"
echo "The output will be one file for each file ending .vcf compared to the file ending .vcf.pivot"
echo "    COMMON markers"
echo "    ONLY .vcf lines"


tag="STARINT"
:<<'COMMENT'
COMMENT

#clear data and output
for d in data output
do
	rm -r ./${d}
	mkdir ${d}
	(cd ${d} || exit
	touch .gitkeep)
done

echo "<<< copy inputs to data >>>"
(cd input || exit
files=$(ls -1 *.vcf*)
parallel \
	echo {}';' \
	cp  {} ../data \
::: ${files}
)

echo "<<< compress and index files >>>"
(cd data || exit
files=$(ls -1 *.vcf*)
parallel \
	echo {}';' \
	rtg bgzip {}';' \
	bcftools index -cf --threads 4 {}.gz \
::: ${files}
)

echo "<<< intersect markers >>>"
(cd data || exit
piv=$(ls -1 *.vcf.pivot.gz)
mv ${piv} hidden
files=$(ls -1 *.gz)
mv hidden ${piv}

name2=$(echo ${piv} | sed s:.vcf.pivot.gz::g)
parallel \
	echo {}';' \
	bcftools isec \
		-c snps \
		-O v \
		-p {/.} \
		--threads 4 \
		${piv} \
		{} ';' \
	cd {/.}';' \
	mv 0003.vcf {/.}_${tag}_${name2}.vcf';' \
	mv {/.}_${tag}_${name2}.vcf ../../output/';' \
	cd ..';' \
	rm -r {/.} \
::: ${files}
)

(cd output || exit
for f in *.vcf_${tag}_*
do
	n=$(echo ${f} | sed "s:.vcf_${tag}_:_${tag}_:g")
	mv ${f} ${n}
done
)

echo "Complete"
