#!/bin/bash
# ==============================================================================
# STAR Alignment — Unsorted BAM + samtools sort
# ==============================================================================
# STAR의 BAM SortedByCoordinate에서 메모리/디스크 에러 발생 시 사용
# ==============================================================================

GENOMEDIR="reference/star_genome"
THREADS=8

ulimit -n 65535

for i in $(cat list_dz.txt); do

    echo "=========================================="
    echo "STAR alignment: $i"
    echo "=========================================="

    # Step 1: STAR — unsorted BAM 출력
    STAR \
        --genomeDir ${GENOMEDIR} \
        --readFilesIn fastq/${i}_R1.fastq.gz fastq/${i}_R2.fastq.gz \
        --readFilesCommand gunzip -c \
        --runThreadN ${THREADS} \
        --twopassMode Basic \
        --alignEndsType Local \
        --alignSoftClipAtReferenceEnds Yes \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterMismatchNmax 10 \
        --outFilterMultimapNmax 20 \
        --outMultimapperOrder Random \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS nM XS MD \
        --outSAMunmapped Within \
        --outFilterScoreMinOverLread 0.5 \
        --outFilterMatchNminOverLread 0.5 \
        --outFileNamePrefix result/${i}

    # Step 2: samtools sort
    echo "  Sorting BAM..."
    samtools sort \
        -@ ${THREADS} \
        -m 2G \
        result/${i}Aligned.out.bam \
        -o result/${i}Aligned.sorted.out.bam

    # Step 3: samtools index
    echo "  Indexing BAM..."
    samtools index result/${i}Aligned.sorted.out.bam

    # Step 4: 원본 unsorted BAM 삭제 (디스크 절약)
    rm -f result/${i}Aligned.out.bam

    echo "Completed: $i"
    echo ""
done

echo "=========================================="
echo "All samples aligned!"
echo "=========================================="
