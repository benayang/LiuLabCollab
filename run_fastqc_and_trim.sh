# conda deactivate
# conda activate py38

fastqc /nas/homes/benyang/LiuLabCollab/7016-BY/fastqs_7016-BY/* -o /nas/homes/benyang/LiuLabCollab/01_fastqc/ -t 16
cd /nas/homes/benyang/LiuLabCollab/01_fastqc/
multiqc .

nohup ~/tools/TrimGalore-0.6.7/trim_galore \
-q 30 \
--fastqc \
--fastqc_args "-t 16 --outdir /nas/homes/benyang/LiuLabCollab/01_fastqc/Trimmed" \
--nextera \
--gzip \
--output_dir /nas/homes/benyang/LiuLabCollab/02_trimmed_fastqs \
--cores 4 \
/nas/homes/benyang/LiuLabCollab/7016-BY/fastqs_7016-BY/* &

cd /nas/homes/benyang/LiuLabCollab/01_fastqc/Trimmed
multiqc .