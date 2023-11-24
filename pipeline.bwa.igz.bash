# Usage: bash pipeline.bwa.igz.bash
sraIndexFile=$HOME/projects/KETCHUP/SlELF3.235C.index

genomeFilesDir=$HOME/genomeFiles/sly/ITAG4.0_release
genomeFastaFile=$genomeFilesDir/S_lycopersicum_chromosomes.4.00.fa
geneGffPath=$genomeFilesDir/ITAG4.0_geneLocus.gff
bwaIndexDir=$HOME/projects/KETCHUP/bwa_index
cleanDataDir=/mnt/petroselinum-2.int.igzev.de/petro-2_data-1_BiSc_2/IGZ_data/raw_data
bwaIndexDir=$HOME/projects/KETCHUP/bwa_index
bwaOutputDir=$HOME/projects/KETCHUP/bwa_output
bamCoverageDir=$HOME/projects/KETCHUP/deeptools/bamCoverage

threads=1

genomeVersion=SL4.0
effectiveGenomeSize=782475302 # SL4.0

if [ -f $bwaIndexDir/$genomeVersion.bwt ];
	then
		echo 'bwa index files already exist!'
	else
		echo 'index files can not be found, now building index in command line as follows...'
		cd $bwaIndexDir
		bwa index -p $genomeVersion $genomeFastaFile
fi

cd $bwaOutputDir


sraIndexArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk -F '\t' '{printf $1" "}'))
sampleNameArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk -F '\t' '{printf $2" "}'))
layoutArray=($(cat $sraIndexFile | grep "ChIP-seq" | awk -F '\t' '{printf $3" "}'))
sampleNumber=$(grep "SINGLE" $sraIndexFile | grep "ChIP-seq" | wc -l)

for ((i=0;i<${sampleNumber};i++))
do
	if [ -f $bwaOutputDir/${sampleNameArray[$i]}.sam ]; then
		echo "$bwaOutputDir/${sampleNameArray[$i]}.sam has already exist in the bwa output dir. Programme is existing to prevent overwriting files......"
		exit
	fi
	if [ ${layoutArray[$i]} == SINGLE ]; then
		bwa mem -t $threads -M $bwaIndexDir/$genomeVersion $cleanDataDir/${sraIndexArray[$i]}_R1.fq.gz > $bwaOutputDir/${sampleNameArray[$i]}.sam
	fi
done

# -M mark shorter split hits as secondary
	# A read may map ambiguously to multiple locations, e.g. due to repeats. Only one of the multiple read alignments is considered primary, and this decision may be arbitrary. All other alignments have the secondary alignment flag.

#-----------------------------------
for ((i=0;i<${sampleNumber};i++))
do
	samtools view -h -bS -@ $threads ${bwaOutputDir}/${sampleNameArray[$i]}.sam > ${bwaOutputDir}/${sampleNameArray[$i]}.unsort.bam
	samtools sort -@ $threads -o ${bwaOutputDir}/${sampleNameArray[$i]}.sort.bam ${bwaOutputDir}/${sampleNameArray[$i]}.unsort.bam
	samtools flagstat -@ $threads ${bwaOutputDir}/${sampleNameArray[$i]}.sort.bam > ${bwaOutputDir}/${sampleNameArray[$i]}.sort.bam.flagstat
	if [ ${layoutArray[$i]} == SINGLE ]; then
		samtools view -h -b -q 30 -F 3588 -@ ${threads} ${bwaOutputDir}/${sampleNameArray[$i]}.sort.bam > ${bwaOutputDir}/${sampleNameArray[$i]}.sort.filter.bam	
		samtools index -@ $threads ${bwaOutputDir}/${sampleNameArray[$i]}.sort.filter.bam
	fi
	rm ${bwaOutputDir}/${sampleNameArray[$i]}.sam
	rm ${bwaOutputDir}/${sampleNameArray[$i]}.unsort.bam
done

cd $bamCoverageDir

for sampleName in ${sampleNameArray[@]}
do
bamCoverage \
-b ${bwaOutputDir}/${sampleName}.sort.filter.bam \
-o ${sampleName}.bigwig \
-of bigwig \
--scaleFactor 1.0 \
--binSize 20 \
--skipNonCoveredRegions \
--smoothLength 100 \
--numberOfProcessors ${threads} \
--effectiveGenomeSize ${effectiveGenomeSize} \
--normalizeUsing RPKM \
--extendReads 100 \
--ignoreDuplicates \
--centerReads
done

