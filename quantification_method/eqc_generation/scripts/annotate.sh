#!/bin/sh

if [ -z "$1" ]; then
    echo "Run only from preprocess.sh" >&2
    exit 1
fi

if [ ! -z "$2" ]; then
    GFILE="$2"
fi

BAM="$1"
BASE=$(basename ${BAM} \.bam)

#GENEANNOT="/grid/mhammell/data/tam/ALS/TE_snRNA_simulation/talitha/TEsolo_processing/db/hg38_refSeqCurated_20191216_transcripts.gtf"
#TEANNOT="/grid/mhammell/data/tam/databases/annotations/GTF/hg38_rmsk_TE.gtf"

#GENEANNOT="/grid/mhammell/data/tam/databases/cole_paper/src/gencode.v40.primary_assembly.annotation_CRfiltered_genes.gtf"
#TEANNOT="/grid/mhammell/data/tam/databases/cole_paper/src/GRCh38_primaryAssembly_rmsk_TEsubfam.gtf"

##this file is just genes, no transcript or exon annotations (not sure why...)
GENEANNOT="/grid/mhammell/data/tam/databases/cole_paper/src/gencode.v40.primary_assembly.annotation_CRfiltered_genes.gtf"
TEANNOT="../../backup/GRCh38_GENCODE_primaryAssembly_rmsk_TEsubfam.gtf"

if [ -z "$GFILE" ]; then
    intersectBed -bed -wa -wb -s -split -a ${BAM} -b ${GENEANNOT} | sed 's/gene_id \"//;s/\";.*$//' | awk -v OFS="	" '{print $4,$21}' > ${BASE}_geneAnnot.txt &
    intersectBed -bed -wa -wb -s -split -a ${BAM} -b ${TEANNOT} | sed 's/gene_id.*transcript_id \"//;s/\".*family_id \"/:/;s/\".*class_id \"/:/;s/\";.*$//' | awk -v OFS="	" '{print $4,$21}' > ${BASE}_TEAnnot.txt &
else
    intersectBed -bed -wa -wb -s -split -a ${BAM} -b ${GENEANNOT} -g ${GFILE} | sed 's/gene_id \"//;s/\";.*$//' | awk -v OFS="	" '{print $4,$21}' > ${BASE}_geneAnnot.txt &
    intersectBed -bed -wa -wb -s -split -a ${BAM} -b ${TEANNOT} -g ${GFILE} | sed 's/gene_id.*transcript_id \"//;s/\".*family_id \"/:/;s/\".*class_id \"/:/;s/\";.*$//' | awk -v OFS="	" '{print $4,$21}' > ${BASE}_TEAnnot.txt &
fi

wait

cat  ${BASE}_geneAnnot.txt ${BASE}_TEAnnot.txt | sort -k1,1 -S 70G -T $PWD --parallel=13 | groupBy -g 1 -c 2 -o distinct > ${BASE}_fullAnnot.txt

if [ $? -ne 0 ]; then
    echo "Error with combining annotations" >&2
    exit 1
else
    echo "No issue"
#    rm ${BASE}_geneAnnot.txt ${BASE}_TEAnnot.txt
fi

echo "Done annotating"
