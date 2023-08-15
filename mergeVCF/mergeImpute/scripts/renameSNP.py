#!/usr/bin/env python3

# %%
import pandas
import sys
import re
import os
import glob

# %%
#standalone usage: renameSNP.py <fileName.vcf> <outfileName.vcf>
fileName = sys.argv[1]
outName = sys.argv[2]

#get count of chromosomes
numChromosomes = 7

# %%
#take header off the top
maxHeader = 70
header = list()
i=0
with open(fileName, 'r') as reader:
    while i < maxHeader:
        line = reader.readline()
        if line.startswith('#'):
            header.append(line)
        else:
            i += 1

# %%
#keep body without header
#...change method soon
f = open(fileName, 'r')
vcf = pandas.read_table(f, skiprows= len(header),
    delim_whitespace=True, skip_blank_lines=True)

# %%
#identify the chromosome indicator character position

#known: in this case col 0 is 1 based and a numpy.int64 type

# %%
#convert SNP name col to only chr num
vcf.head(1)
print("early x")
exit()
vcf.iloc[:,0] = vcf.iloc[:,0].astype('int')

# %%
#make full SNP name
vcf.iloc[:,2] = "Chr" + vcf.iloc[:,0].astype(str) + "H_" + vcf.iloc[:,1].astype(str)

# %%
#write header and dataframe to new file
with open(outName, 'w') as writer:
    for line in header:
        writer.write(line)
    writer.write(vcf.to_csv(sep="\t", index=False, header=False))
