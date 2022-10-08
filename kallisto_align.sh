# conda deactivate
# conda activate py38

# for setting nextera transposase fragment length and sd
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08316-y

pjdir="/nas/homes/benyang/LiuLabCollab"
fastqdir="$pjdir/02_trimmed_fastqs"
outdir="$pjdir/03_Kallisto"

cd $fastqdir

# for i in 7016-BY-3_AGGCAGAA-CTCTCTAT_S3_R1_001 7016-BY-4_TCCTGAGC-CTCTCTAT_S4_R1_001 7016-BY-5_GGACTCCT-CTCTCTAT_S5_R1_001 7016-BY-6_TAGGCATG-CTCTCTAT_S6_R1_001
for i in `ls -1 *_trimmed.fq.gz | sed 's/_trimmed.fq.gz//'`
do
    echo @@@@@@ Processing "$i"

    kallisto quant \
    -b 100 \
    --single \
    -l 300 \
    -s 30 \
    -i "/nas/homes/benyang/Genome_References/Kallisto/Homo_sapiens.GRCh37.cdna.all.fa.idx" \
    -t 45 \
    -o $outdir/$i/ \
    $fastqdir/${i}_trimmed.fq.gz

    # --genomebam --gtf "/nas/homes/benyang/Genome_References/Kallisto/Homo_sapiens.GRCh37.87.gtf.gz" \
    # --chromosomes "/nas/homes/benyang/Genome_References/Kallisto/ensembl.hg19.chrom.sizes" \


    echo @@@@@@@ Saved to "$outdir/$i" 
done