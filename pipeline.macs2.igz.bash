# Usage: bash pipeline.macs2.igz.bash $sraIndexFile
bamFilesDir=$HOME/projects/KETCHUP/bwa_output
macs2OutputDir=$HOME/projects/KETCHUP/macs2Output/SlELF3
effectiveGenomeSize=7.8e8 # SL4.0
sraIndexFile=$1
IPsampleNamesArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk '$0!~/INPUT/ && $0!~/Input/ && $0!~/input/{printf $2" "}'))

sraIndexArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk -F '\t' '{printf $1" "}'))
sampleNameArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk -F '\t' '{printf $2" "}'))
layoutArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk -F '\t' '{printf $3" "}'))
sampleNumber=$(grep -E "SINGLE" $sraIndexFile | grep "ChIP-seq" | wc -l)

cd $macs2OutputDir

inputBam=$(echo "${bamFilesDir}/Input_SlMM_ZT16_30C_rep1.sort.filter.bam ${bamFilesDir}/Input_SlMM_ZT16_30C_rep2.sort.filter.bam ${bamFilesDir}/Input_SlMM_ZT16_37C_rep2.sort.filter.bam")
for sampleName in ${IPsampleNamesArray[@]}
do
	macs2 callpeak -t ${bamFilesDir}/${sampleName}.sort.filter.bam -c $inputBam -f BAM -g ${effectiveGenomeSize} --outdir ${macs2OutputDir} -n ${sampleName} --call-summits -q 0.05
done


