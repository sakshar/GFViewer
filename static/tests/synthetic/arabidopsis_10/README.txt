Synthetic dataset: 10 randomly generated gene families placed on the five nuclear chromosomes of Arabidopsis thaliana (TAIR10 lengths). Each family follows one of four spatial patterns (tandem array, sub-telomeric, single-chromosome-enriched, dispersed) so the analytics have something to find. Not real genes.

Files:
  genes.tsv    gene_id, gene_family, chromosome, start, end, strand (1-based, inclusive)
  genome.txt   chromosome,length
  colors.txt   gene_family,palette-index (optional)

Per-family spatial pattern used by the generator:
  AT_GF01      dispersed
  AT_GF02      tandem
  AT_GF03      subtelomeric
  AT_GF04      dispersed
  AT_GF05      dispersed
  AT_GF06      tandem
  AT_GF07      tandem
  AT_GF08      subtelomeric
  AT_GF09      dispersed
  AT_GF10      chromosome-enriched
