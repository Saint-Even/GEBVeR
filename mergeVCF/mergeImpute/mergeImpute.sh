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
# designed to take the merged.vcf files from cassetteGBS output
# merge plates to a large training data set
# this super plate is then imputed as a whole
#
# !!!The output of this pipeline gave a null result when run through star intersect, investigate before using!!!

#activate an environment with:
#parallel
#rtgtools
#bcftools
#vcftools
#beagle

echo "Begin: all marker merge"

#first construct final name
cd input/
files=$(ls *.vcf)
echo ${files}
clean=$(echo ${files} | sed s:.vcf::g)

name=""
for n in ${clean}
do
	  name+="${n}_"
done
name+="combined.vcf"
cd ..

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
files=$(ls -1 *.vcf)
parallel \
	echo {}';' \
	cp  {} ../data \
::: ${files}
)

echo "<<< compress and index files >>>"
(cd data || exit
files=$(ls -1 *.vcf)
parallel \
	echo {}';' \
	rtg bgzip {}';' \
	bcftools index -cf --threads 4 {}.gz \
::: ${files}
)

echo "<<< convert >>>"
(cd data || exit
 files=$(ls -1 *.vcf.gz | sed 's/.vcf.gz//g')
 parallel \
	   echo {}';' \
     bcftools convert --threads 4 -O b -o {}.bcf {}.vcf.gz \
     ::: ${files}
)

echo "<<< sort >>>"
(cd data || exit
 files=$(ls -1 *.bcf)
 parallel \
	   echo {}';' \
     bcftools sort -O b -o {}.sort {}';' \
	   bcftools index -cf --threads 4 {}.sort \
     ::: ${files}
)

echo "<<< merge >>>"
#bcftools concat --threads 12 -o variantsMerged.vcf -O v ${files}
(cd data || exit
files=$(ls -1 *.bcf.sort)
bcftools merge --force-samples --threads 12 -o variantsMerged.vcf -O v ${files}
)

echo "<<< cleanup >>>"
(cd data || exit
 rm *.bcf*
 rm *.vcf.*
)

echo "<<< filter >>>"
(cd data || exit
output="variantsFiltered.vcf"
vcftools \
    --vcf variantsMerged.vcf \
    --remove-filtered-all \
    --max-missing 0.2 \
    --maf 0.01 \
    --remove-indels \
    --mac 1 \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --out ${output}

mv ${output}.recode.vcf ${output}
)

echo "<<< renameVariants >>>"
(cd data || exit
fileName="variantsFiltered.vcf"
outName="variantsFilteredRenamed.vcf"
python ../scripts/renameSNP.py ${fileName} ${outName}
)

echo "<<< impute >>>"
(cd data || exit
#input="variantsFilteredRenamed.vcf"
input="variantsFiltered.vcf"
output="variantsImputed.vcf"
beagle -Xmx16G \
    gt=${input} \
    out=${output}

rm *.log
gzip -d ${output}.vcf.gz
mv ${output}.vcf ${name}
echo ${name}
)

#copy result to output
(cd data || exit
name=$(ls *combined.vcf )
cp ${name} ../output
)

echo "Complete"
exit
