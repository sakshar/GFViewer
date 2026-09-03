The 6 Babesia sp. MO1 gene families (see data_test_2.csv) written in every input format GFViewer accepts.

  genes.bed         one BED file  -> column 4 (name) is the family
  per_family/*.bed  several BED files -> each file name is a family
  genes_named.bed   BED with gene ids in column 4 ...
  mapping.tsv       ... + this gene_id -> gene_family map
  genes.gff3        GFF3, family in the gene_family attribute
  genes.gtf         GTF, family in the gene_family attribute

Genome for all of the above: ../chrs_test_2.txt
