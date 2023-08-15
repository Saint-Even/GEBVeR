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
# designed to take the *.vcf files from starIntersect or starCombine
# filter for repeat lines
# drop all empty lines
# drop all customDrop lines

#activate an environment with:
#rtgtools
#bcftools
#vcftools
#parallel

#USER: enter the list of lines that will be dropped, space separated with no preceeding or trailing space
# eg.
# customDrop="A_Line AC_dc Special420-42" or
customDrop=""
#customDrop="AAC_Connect AC_Bountiful AC_Oxbow CDC_Austenson CDC_Copper CDC_Kendall CDC_Reserve Cerveza Major Newdale Norman"

echo "Begin: Line Filter"

:<<'MASK'
MASK

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

echo "<<< filter for repeat names >>>"
(cd data || exit
 for f in *.vcf
 do
     echo ${f}
     allNames=$(bcftools query -l ${f})
     cleanNames=""
     colonNames=""
     dropNames=""
     keepNames=""

     #subset potential duplicates
     for i in ${allNames}
     do
         if [[ "${i}" == *:* ]]; then
             colonNames+="${i} "
         else
             cleanNames+="${i} "
         fi
     done

     #test for duplication
     for i in ${cleanNames}
     do
         for j in ${colonNames}
         do
             # i as rootname, j as [digit]:rootname
             # ie. i gets treated and tested against j as is
             if [[ "${j}" == *:${i} ]]; then
                 #then j goes to dropName
                 dropNames+="${j} "
             fi
         done
     done

     # add to all plates
     dropNames+="empty Empty EMPTY"
     dropNames+=" ${customDrop}"

     #remove dropNames from allNames to get keepNames
     for i in ${allNames}
     do
         drop="false"
         for j in ${dropNames}
         do
             #test if kicking out
             if [[ "${i}" == "${j}" ]]; then
                 drop="true"
                 echo "dropping: ${i}"
             fi
         done
         if [ "${drop}" == "false" ]; then
             keepNames+="${i},"
         fi
     done
     #drop last comma
     keepNames=${keepNames::-1}

     #process original .vcf to retain keepNames
	   rtg bgzip ${f}
	   bcftools index -cf --threads 4 ${f}.gz
     bcftools view \
              -O v \
              -o ${f}.noDups \
              --threads 4 \
              -s ${keepNames} \
              ${f}.gz

     #replace original with clean file...
     mv ${f}.noDups ${f}

 done # end files loop
)

#copy result to output
(cd data || exit
files=$(ls -1 *.vcf)
parallel \
    bcftools stats -v --threads 4 {} ">" {/.}_bcf_stats.txt';' \
    rtg vcfstats {} ">" {/.}_rtg_stats.txt';' \
    cp {/.}_rtg_stats.txt {/.}_bcf_stats.txt {} ../output \
    ::: ${files}
)

echo "Complete"
exit
