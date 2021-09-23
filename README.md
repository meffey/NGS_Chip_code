# NGS_Chip-seq_code_perl

## bed_overlap.pl:
perl bed_overlap.pl input_a.bed input_b.bed out.bed

## extract_seq.pl：
perl extract_seq.pl genome.fa peaks.bed > out.fa

## peaks_annotation.pl:
perl peaks_annotation.pl genes.gff peaks.bed

## peaks_annotation_tss.pl"
perl peaks_annotation_tss.pl genes_transposons.gff peaks.bed
