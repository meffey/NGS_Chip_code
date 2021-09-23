# NGS_Chip-seq_code_perl

## bed_overlap.pl:
   Description: This script find overlap peaks between two peaks. bed format is input type.
### perl bed_overlap.pl input_a.bed input_b.bed out.bed

## extract_seq.pl：
   Description: This script extract equally length sequence according to peaks.
### perl extract_seq.pl genome.fa peaks.bed length out.fa

    For example: perl extract_seq.pl A.fa A_peaks.bed 1000 A_1000.fa

## peaks_annotation.pl:
   Description: This script annotated peaks according to relationship between peaks and gene(include exon,intron,3'UTR,5'UTR and promote).
### perl peaks_annotation.pl genes.gff peaks.bed

## peaks_annotation_tss.pl:
   Description: This script annotated peaks according to relationship between peaks and transcription start sites(TSS-400bp).
### perl peaks_annotation_tss.pl genes_transposons.gff peaks.bed
