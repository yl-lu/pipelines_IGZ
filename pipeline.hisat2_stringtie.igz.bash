# Usage: bash pipeline.hisat2_stringtie.igz.bash $relpathSraIndexFile

genomeFilesDir=/net/home/igz/yllu/genomeFiles/sly/ITAG4.0_release
genomeFastaFile=$genomeFilesDir/S_lycopersicum_chromosomes.4.00.fa
geneGffPath=$genomeFilesDir/ITAG4.0_gene_models.gff
locusGffPath=$genomeFilesDir/ITAG4.0_geneLocus.gff
ssFile=$genomeFilesDir/ITAG4.0.ss
exonFile=$genomeFilesDir/ITAG4.0.exon
cleanDataDir=/mnt/petroselinum-2.int.igzev.de/petro-2_data-1_BiSc_2/IGZ_data/raw_data/216R/raw
hisatIndexDir=/net/home/igz/yllu/projects/KETCHUP/hisat2_index
hisatOutputDir=/net/home/igz/yllu/projects/KETCHUP/hisat2_output
stringtieOutputDir=/net/home/igz/yllu/projects/KETCHUP/stringtie_output
genomeVersion=SL4.0

threads=1

if [ -f $hisatIndexDir/$genomeVersion.1.ht2 -o -f $hisatIndexDir/$genomeVersion.1.ht2l ];
	then
		echo 'hisat2 index files were found.'
	else
		echo 'index files can not be found, please build index first...'
		echo 'Command line may look like follows:'
		echo "hisat2-build -p $threads $genomeFastaPath $hisatIndexDir/$genomeVersion"
fi

relpathSraIndexFile=$1
sraIndexFile=$(readlink -e ${relpathSraIndexFile})
echo "relpathSraIndexFile: $relpathSraIndexFile"
echo "sraIndexFile: $sraIndexFile"
if [ ! -f "$sraIndexFile" ]; then
	echo "Can not find index file ${sraIndexFile}!"
	exit
fi
sraIndexArray=($(cat $sraIndexFile | grep "RNA-seq" | awk -F '\t' '{printf $1" "}'))
sampleNameArray=($(cat $sraIndexFile | grep "RNA-seq" | awk -F '\t' '{printf $2" "}'))
layoutArray=($(cat $sraIndexFile | grep "RNA-seq" | awk -F '\t' '{printf $3" "}'))
sampleNumber=$(grep "RNA-seq" $sraIndexFile | wc -l)

cd $hisatOutputDir

for ((i=0;i<${sampleNumber};i++));
do
	if [ -f $hisatOutputDir/${sampleNameArray[$i]}.sort.filter.bam ]; then
		echo "$hisatOutputDir/${sampleNameArray[$i]}.sort.filter.bam already exist! hisat2 alignment will NOT run..."
		break
	else
		if [ ${layoutArray[$i]} == PAIRED ]; then
			echo "hisat2 alignment for paired reads..."
			hisat2 -p $threads --summary-file $hisatOutputDir/${sampleNameArray[$i]}.hisat2align_summary --dta --no-mixed --no-discordant -t -x $hisatIndexDir/$genomeVersion -1 $cleanDataDir/${sraIndexArray[$i]}_R1.fq.gz -2 $cleanDataDir/${sraIndexArray[$i]}_R2.fq.gz -S $hisatOutputDir/${sampleNameArray[$i]}.sam
		fi
	fi
done

# --no-mixed
#	 By default, when hisat2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.
# --no-discordant
#	 By default, hisat2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior.

# bam sort
cd $hisatOutputDir

for ((i=0;i<${sampleNumber};i++));
do
	if [ -f ${hisatOutputDir}/${sampleNameArray[$i]}.sort.filter.bam ]; then
		echo "${hisatOutputDir}/${sampleNameArray[$i]}.sort.filter.bam already exist! bamSort and bamFilter will NOT run..."
		break
	else
		samtools view -h -bS -@ $threads ${hisatOutputDir}/${sampleNameArray[$i]}.sam > ${hisatOutputDir}/${sampleNameArray[$i]}.unsort.bam
		samtools sort -@ $threads -o ${hisatOutputDir}/${sampleNameArray[$i]}.sort.bam ${hisatOutputDir}/${sampleNameArray[$i]}.unsort.bam
		samtools flagstat -@ $threads ${hisatOutputDir}/${sampleNameArray[$i]}.sort.bam > ${hisatOutputDir}/${sampleNameArray[$i]}.sort.bam.flagstat
		if [ ${layoutArray[$i]} == PAIRED ]; then
			samtools view -h -b -q 30 -f 2 -F 1292 -@ $threads ${hisatOutputDir}/${sampleNameArray[$i]}.sort.bam > ${hisatOutputDir}/${sampleNameArray[$i]}.sort.filter.bam
			samtools index -@ $threads ${hisatOutputDir}/${sampleNameArray[$i]}.sort.filter.bam
		fi
	    if [ ${layoutArray[$i]} == SINGLE ]; then
			samtools view -h -b -q 30 -F 3588 -@ ${threads} ${hisatOutputDir}/${sampleNameArray[$i]}.sort.bam > ${hisatOutputDir}/${sampleNameArray[$i]}.sort.filter.bam
			samtools index -@ $threads ${hisatOutputDir}/${sampleNameArray[$i]}.sort.filter.bam
		fi
	fi
	rm ${hisatOutputDir}/${sampleNameArray[$i]}.sam
	rm ${hisatOutputDir}/${sampleNameArray[$i]}.unsort.bam
done

# read counts
bamFiles=$(for sampleName in ${sampleNameArray[@]}; do echo ${sampleName}.sort.filter.bam; done)

cd $hisatOutputDir


# reshape readsCount matrix for edgeR

# stringtie
cd $stringtieOutputDir

for ((i=0;i<${sampleNumber};i++));
do
	if [ -f $stringtieOutputDir/${sampleNameArray[$i]}.cov_refs.gtf ]; then
		echo "$stringtieOutputDir/${sampleNameArray[$i]}.cov_refs.gtf already generated, stringtie will NOT run..."
		break
	else
		if [ ${layoutArray[$i]} == PAIRED ]; then
			stringtie -p ${threads} $hisatOutputDir/${sampleNameArray[$i]}.sort.filter.bam -o $stringtieOutputDir/${sampleNameArray[$i]}.assembled_transcripts.gtf -G $geneGffPath -l ${sampleNameArray[$i]} -A $stringtieOutputDir/${sampleNameArray[$i]}.gene_abund.tab -C $stringtieOutputDir/${sampleNameArray[$i]}.cov_refs.gtf -e
		fi
	fi
done
