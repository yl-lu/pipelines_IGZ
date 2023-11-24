# Usage: bash pipeline.macs2.filter.igz.bash
bamFilesDir=$HOME/projects/KETCHUP/bwa_output
macs2OutputDir=$HOME/projects/KETCHUP/macs2Output/SlELF3

sraIndexFile=$HOME/projects/KETCHUP/SlELF3.235C.index
genomeDir=$HOME/genomeFiles/sly/ITAG4.0_release

IPsampleNamesArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk '$0!~/input/ && $0!~/Input/ && $0!~/INPUT/{printf $2" "}'))
allsampleNamesArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk '{printf $2" "}'))

cd $macs2OutputDir

# qvalue < 0.05 as the cutoff for filtering peaks
for sampleName in ${IPsampleNamesArray[@]}
do
	grep -v "#" ${sampleName}_peaks.xls | awk -F '\t' '$9>1.30103{OFS="\t"; print $1,$2-1,$3,$4,$5,$6,$7,$8,$9,$10}' > ${sampleName}_peaks.filter.bed
	cat ${sampleName}_peaks.filter.bed | awk -F '\t' '{OFS="\t"; print $1,$2,$3}' | sort -k1,1 -k2,2n | uniq > ${sampleName}_peaks.filter.uniq.bed
	bedtools intersect -a $genomeDir/ITAG4.0_PGD_2kb.bed -b ${sampleName}_peaks.filter.bed -c | awk -F '\t' '$7>0{OFS="\t"; print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n | uniq > ${sampleName}.targetGenes.bed
done

