- Assembly with MegaHit and Spades
```
module load megahit/1.1.3
module load spades/3.11.1
megahit -1 $reads/QUI85_MG-1.fastq -2 $reads/QUI85_MG-2.fastq --presets meta-large -t 56 -o $reads/MegaHit_QUI85_MG -m 150000000000000
megahit -1 $reads/SJ75_MG-1.fastq -2 $reads/SJ75_MG-2.fastq --presets meta-large -t 56 -o $reads/MegaHit_SJ75_MG -m 150000000000000
spades.py -o $reads/Spades_QUI85_MG --meta -1 $reads/QUI85_MG-1.fastq -2 $reads/QUI85_MG-2.fastq -t 56 -k 27,37,47,57,67,77,87,97,107,117,127 --only-assembler
spades.py -o $reads/Spades_SJ75_MG --meta -1 $reads/SJ75_MG-1.fastq -2 $reads/SJ75_MG-2.fastq-t 56 -k 27,37,47,57,67,77,87,97,107,117,127 --only-assembler

```
- Contig QC
```
module load bowtie2/2.2.5
module load samtools/1.8
module load bedtools2/2.24.0
bowtie2-build -f $rawcontig/MegaHit_QUI85_MG.fa $rawcontig/MegaHit_QUI85_MG.fa.build
bowtie2-build -f $rawcontig/MegaHit_SJ75_MG.fa $rawcontig/MegaHit_SJ75_MG.fa.build
bowtie2 --very-sensitive-local -x $rawcontig/MegaHit_QUI85_MG.fa.build -1 $reads/QUI85_MG-1.fastq -2 $reads/QUI85_MG-2.fastq -S $rawcontig/aln-QUI85 -p 28
bowtie2 --very-sensitive-local -x $rawcontig/MegaHit_SJ75_MG.fa.build -1 $reads/SJ75_MG-1.fastq -2 $reads/SJ75_MG-2.fastq -S $rawcontig/aln-SJ75 -p 28
samtools view -@ 28 -b -S -o $rawcontig/aln-QUI85.bam $rawcontig/aln-QUI85
samtools view -@ 28 -b -S -o $rawcontig/aln-SJ75.bam $rawcontig/aln-SJ75
samtools sort -@ 28 $rawcontig/aln-SJ75.bam -o $rawcontig/aln-SJ75.bam.sorted
samtools sort -@ 28 $rawcontig/aln-QUI85.bam -o $rawcontig/aln-QUI85.bam.sorted
for a in *fa; do ./seq_lengths $a | sed 's/ flag.*len.* / /' > $a.lenghts;done
genomeCoverageBed -d -split -ibam $rawcontig/aln-SJ75.bam.sorted -g $rawcontig/MegaHit_SJ75_MG.fa.lenghts | ~/02_Scripts/contig_cov_from_bed > $rawcontig/SJ75.fa.non_zero_cov_perc.txt
genomeCoverageBed -d -split -ibam $rawcontig/aln-QUI85.bam.sorted -g $rawcontig/MegaHit_QUI85_MG.fa.lenghts | ~/02_Scripts/contig_cov_from_bed > $rawcontig/QUI85.fa.non_zero_cov_perc.txt
awk '$3 >= 90 { print $2,$1,$3 }' $rawcontig/SJ75.fa.non_zero_cov_perc.txt | sort -nr | awk 'BEGIN { OFS=" " } { print $2,$1,$3 }' > $rawcontig/SJ75.fa.txt
awk '$3 >= 90 { print $2,$1,$3 }' $rawcontig/QUI85.fa.non_zero_cov_perc.txt| sort -nr | awk 'BEGIN { OFS=" " } { print $2,$1,$3 }' > $rawcontig/QUI85.fa.txt
~/02_Scripts/grep_ids $rawcontig/SJ75.fa.txt $rawcontig/MegaHit_SJ75_MG.fa > $rawcontig/SJ75_MG.QC.fa
~/02_Scripts/grep_ids $rawcontig/QUI85.fa.txt $rawcontig/MegaHit_QUI85_MG.fa > $rawcontig/QUI85_MG.QC.fa

```
- Assembly Quality
```
module load quast/5.0.0 
quast.py -o $contig/  --contig-thresholds 2000,5000,10000,25000,50000 QUI85_MG.QC.fa -t 28
quast.py -o $contig/  --contig-thresholds 2000,5000,10000,25000,50000 SJ75_MG.QC.fa -t 28
```
- Generate Coverage Profiles. Using differential coverage as a metric for clustering contings into bins requires specific file formats for each of the binning programs. These were built using the following commands
```
module load bowtie2/2.2.5
module load samtools/1.8
module load bedtools2/2.24.0
module load metabat/2.12.1
bowtie2-build -f $raw/QUI85_MG.QC.fa $bam/QUI85_MG.QC.fa.build
bowtie2-build -f $raw/SJ75_MG.QC.fa $bam/SJ75_MG.QC.fa.build
bowtie2 --very-sensitive-local -x $bam/QUI85_MG.QC.fa.build -1 $reads/QUI85_MG-1.fastq -2 $reads/QUI85_MG-2.fastq  -S $bam/aln-QUI85 -p 28
bowtie2 --very-sensitive-local -x $bam/SJ75_MG.QC.fa.build -1 $reads/SJ75_MG-1.fastq -2 $reads/SJ75_MG-2.fastq -S $bam/aln-SJ75 -p 28
samtools view -@ 28 -b -S -o $bam/aln-QUI85.bam $bam/aln-QUI85
samtools view -@ 28 -b -S -o $bam/aln-SJ75.bam $bam/aln-SJ75
samtools sort -@ 28 $bam/aln-SJ75.bam -o $bam/aln-SJ75.bam.sorted
samtools sort -@ 28 $bam/aln-QUI85.bam -o $bam/aln-QUI85.bam.sorted
jgi_summarize_bam_contig_depths --outputDepth $bam/SJ75.txt $bam/aln-SJ75.bam.sorted
jgi_summarize_bam_contig_depths --outputDepth $bam/QUI85.txt $bam/aln-QUI85.bam.sorted
metabat -i $raw/QUI85_MG.QC.fa -o $bam/QUI85/QUI85 -a $bam/QUI85.txt -m 2000 -t 28
metabat -i $raw/SJ75_MG.QC.fa -o $bam/SJ75/SJ75 -a $bam/SJ75.txt -m 2000 -t 28
```
- Bin with MaxBin with 107 and 40 SCG
```
module load bowtie2/2.2.5
module load perl/5.26.0
mkdir $output/SJ75_maxbin.107
run_MaxBin.pl -contig $fasta/SJ75_MG.QC.fa -out $output/SJ75_maxbin.107/SJ75.107 -abund $output/SJ75.txt  -min_contig_length 2000 -thread 28
mkdir $output/QUI85_maxbin.107
run_MaxBin.pl -contig $fasta/QUI85_MG.QC.fa -out $output/QUI85_maxbin.107/QUI85.107 -abund $output/QUI85.txt  -min_contig_length 2000 -thread 28
mkdir $output/SJ75_maxbin.40
run_MaxBin.pl -contig $fasta/SJ75_MG.QC.fa -out $output/SJ75_maxbin.40/SJ75.40 -abund $output/SJ75.txt  -min_contig_length 2000 -thread 28 -markerset 40
mkdir $output/QUI85_maxbin.40
run_MaxBin.pl -contig $fasta/QUI85_MG.QC.fa -out $output/QUI85_maxbin.40/QUI85.40 -abund $output/QUI85.txt  -min_contig_length 2000 -thread 28 -markerset 40
```
- CheckM and Anvi'o refinement
```
module load anaconda3/5.3.0
module load perl/5.26.0
module load samtools/1.8
module load checkm/1.0.11
module load hmmer/3.1b2 
conda activate anvio-6.1

checkm lineage_wf -x fa -t 28 --pplacer_threads 28 QUI_Metabat_Bins/ ../00_CheckM/QUI_Metabat_Bins
checkm lineage_wf -x fa -t 28 --pplacer_threads 28 SJ_Metabat_Bins ../00_CheckM/SJ_Metabat_Bins
checkm qa -o 2 --tab_table -t 28 -f QUI_Metabat_Bins.txt QUI_Metabat_Bins/lineage.ms QUI_Metabat_Bins
checkm qa -o 2 --tab_table -t 28 -f SJ_Metabat_Bins.txt SJ_Metabat_Bins/lineage.ms SJ_Metabat_Bins

checkm lineage_wf -x fa -t 28 --pplacer_threads 28 QUI_MaxBin_Bins/ ../00_CheckM/QUI_MaxBin_Bins
checkm lineage_wf -x fa -t 28 --pplacer_threads 28 SJ_MaxBin_Bins ../00_CheckM/SJ_MaxBin_Bins
checkm qa -o 2 --tab_table -t 28 -f QUI_MaxBin_Bins.txt QUI_MaxBin_Bins/lineage.ms QUI_MaxBin_Bins
checkm qa -o 2 --tab_table -t 28 -f SJ_MaxBin_Bins.txt SJ_MaxBin_Bins/lineage.ms SJ_MaxBin_Bins

anvi-script-reformat-fasta QUI85_MG.QC.fa -o QUI85-fixed.fa -l 2000 --simplify-names -r QUI85-fixed.txt
anvi-script-reformat-fasta SJ75_MG.QC.fa -o SJ75-fixed.fa -l 2000 --simplify-names -r SJ75-fixed.txt
anvi-gen-contigs-database -f $anvio/QUI85-fixed.fa -o $anvio/QUI85_contigs.db
anvi-gen-contigs-database -f $anvio/SJ75-fixed.fa -o $anvio/SJ75_contigs.db
anvi-run-hmms -c $anvio/QUI85_contigs.db -T 28
anvi-run-hmms -c $anvio/SJ75_contigs.db -T 28
anvi-init-bam $anvio/aln-QUI85.sorted.bam -o $anvio/aln-QUI85.sorted.anvio.bam
anvi-init-bam $anvio/aln-SJ75.sorted.bam -o $anvio/aln-SJ75.sorted.anvio.bam
anvi-profile -i $anvio/aln-SJ75.sorted.anvio.bam -c $anvio/SJ75_contigs.db -o $anvio/SJ75_PROFILE --skip-SNV-profiling --cluster-contigs -T 28
anvi-profile -i $anvio/aln-QUI85.sorted.anvio.bam -c $anvio/QUI85_contigs.db -o $anvio/QUI85_PROFILE --skip-SNV-profiling --cluster-contigs -T 28

```
