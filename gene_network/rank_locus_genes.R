#!/usr/bin/env Rscript
library(GenomicRanges)
library(pheatmap)
library(annotables)
library(liftOver)
library(tidyverse)

###############################################################################
# Get clear locus definition from merging OT locus2gene score from two
# IBD GWAS - Jimmy Liu's and Katrina de Lange's.
dir = "/Users/jeremys/work/opentargets/gene_network/ibd"
l2g.katie.df = read_tsv(file.path(dir, "l2g.IBD.deLange.tsv")) %>% rename(chr = chrom)
l2g.jimmy.df = read_tsv(file.path(dir, "l2g.IBD.Liu.tsv")) %>% rename(chr = chrom)

l2g.df = bind_rows(l2g.katie.df, l2g.jimmy.df) %>%
  arrange(chr, pos) %>%
  mutate(id = paste(study_id, chr, pos, ref, alt)) %>%
  dplyr::select(id, gene_id, y_proba_full_model)

l2g.df.spread = l2g.df %>%
  spread(key = "id", value = y_proba_full_model)
l2g.df.spread[is.na(l2g.df.spread)] = 0

# Look at a heatmap that shows the correlation between loci, based on the 
# overlap between their nearby genes. Some loci are unique to one study,
# and some are shared.
l2g.df.cor = cor(l2g.df.spread %>% dplyr::select(-gene_id))
pheatmap(l2g.df.cor, treeheight_row = 0, treeheight_col = 0, show_colnames = F, fontsize = 5)

# Make a table with all genes at each locus in a single line
grch38_nodup = grch38 %>% filter(!duplicated(ensgene))
l2g.df = bind_rows(l2g.katie.df, l2g.jimmy.df) %>%
  arrange(chr, pos) %>%
  mutate(id = paste(study_id, chr, pos, ref, alt, sep = "_")) %>%
  dplyr::select(study_id, id, chr, l2g_pos = pos, ref, alt, gene_id, l2g_score = y_proba_full_model) %>%
  left_join(grch38_nodup %>% dplyr::select(ensgene, symbol), by=c("gene_id" = "ensgene"))

l2g.locus.df = l2g.df %>%
  arrange(id, desc(l2g_score)) %>%
  group_by(id) %>%
  summarise(chr = dplyr::first(chr),
            l2g_pos = dplyr::first(l2g_pos),
            ref = dplyr::first(ref), 
            alt = dplyr::first(alt),
            top_gene = dplyr::first(symbol),
            top_gene_score = dplyr::first(l2g_score),
            locus_genes = paste(sprintf("%s_%g", symbol, l2g_score), collapse=(", ")),
            locus_genes_ens = paste(sprintf("%s_%g", gene_id, l2g_score), collapse=(", "))) %>%
  arrange(chr, l2g_pos)

# Save this file, which has all the l2g scores for each locus **for each study** (Jimmy's, Katie's)
write_tsv(l2g.locus.df, path = file.path(dir, "l2g.IBD_loci.by_study.tsv"))


# Merge together nearby loci (within 200 kb)
# To do this we first make a column with a "locus ID" that is the
# chr:pos of the first signal in a region
l2g.locus.mapping.df = l2g.locus.df %>%
  arrange(chr, l2g_pos) %>%
  mutate(l2g_locus = paste(chr, l2g_pos, sep = "_")) %>%
  ungroup() %>%
  mutate(prev_row_dist = l2g_pos - lag(l2g_pos)) %>%
  rowwise() %>%
  mutate(prev_row_dist = max(0, prev_row_dist)) %>%
  mutate(l2g_locus_merged = l2g_locus)
for (i in 2:nrow(l2g.locus.mapping.df)) {
  if (l2g.locus.mapping.df$prev_row_dist[i] < 2e5 & l2g.locus.mapping.df$chr[i] == l2g.locus.mapping.df$chr[i-1]) {
    l2g.locus.mapping.df$l2g_locus_merged[i] = l2g.locus.mapping.df$l2g_locus_merged[i-1]
  }
}
l2g.locus.mapping.df = l2g.locus.mapping.df %>%
  ungroup() %>%
  filter(!duplicated(l2g_locus))

# The unique loci remaining
nrow(l2g.locus.mapping.df %>% filter(!duplicated(l2g_locus_merged)))

# Now group by the new locus_id and get summary values for all genes at each locus.
# E.g. mean score, max score
l2g.df = l2g.df %>%
  mutate(l2g_locus = paste(chr, l2g_pos, sep = "_")) %>%
  left_join(l2g.locus.mapping.df %>% select(l2g_locus, l2g_locus_merged), by="l2g_locus")

l2g.merged.df = l2g.df %>%
  arrange(desc(l2g_score)) %>%
  group_by(l2g_locus_merged, gene_id) %>%
  summarise(symbol = first(symbol),
            id = first(id),
            chr = first(chr),
            l2g_pos = first(l2g_pos),
            ref = first(ref),
            alt = first(alt),
            l2g_score_mean = mean(l2g_score),
            l2g_score_max = max(l2g_score))
# This table has the genes associated with the merged loci. This is useful later so let's save it.
l2g.merged.df = l2g.merged.df %>%
  select(l2g_locus_merged, chr, l2g_pos, ref, alt, gene_id, symbol, l2g_score_mean, l2g_score_max)
write_tsv(l2g.merged.df, path = file.path(dir, "l2g.IBD_genes.merged.tsv"))

# This table now has the average (and max) L2G score for each gene at each of our
# merged loci.
l2g.locus.merged.df = l2g.merged.df %>%
  arrange(desc(l2g_score_mean)) %>%
  group_by(l2g_locus_merged) %>%
  summarise(chr = dplyr::first(chr),
            l2g_mean_pos = mean(l2g_pos),
            l2g_pos = dplyr::first(l2g_pos),
            ref = dplyr::first(ref), 
            alt = dplyr::first(alt),
            top_gene_mean = dplyr::first(symbol),
            top_gene_score_mean = dplyr::first(l2g_score_mean),
            locus_genes_mean = paste(sprintf("%s_%g", symbol, l2g_score_mean), collapse=(", ")),
            locus_genes_ens_mean = paste(sprintf("%s_%g", gene_id, l2g_score_mean), collapse=(", ")),
            locus_genes_max = paste(sprintf("%s_%g", symbol, l2g_score_max), collapse=(", ")),
            locus_genes_ens_max = paste(sprintf("%s_%g", gene_id, l2g_score_max), collapse=(", ")))


l2g.gr = makeGRangesFromDataFrame(l2g.locus.merged.df %>% dplyr::select(chr, pos = l2g_pos, l2g_locus_merged), keep.extra.columns = T, ignore.strand = T,
                                  seqinfo = NULL, seqnames.field = "chr", start.field = "pos", end.field = "pos")


###############################################################################
# Annotated loci from de Lange Supp Table 2
locus.df = read_tsv(file.path(dir, "deLange.locus_table.tsv"))
# ibd.locus.df = locus.df %>%
#   filter(grepl("IBD", Trait))
# Rethinking this - don't filter at all, so that we retain as many annotated loci as possible
ibd.locus.df = locus.df
ibd.locus.gr.hg19 = makeGRangesFromDataFrame(ibd.locus.df, keep.extra.columns = T, ignore.strand = T,
                                             seqinfo = NULL, seqnames.field = "chr", start.field = "pos", end.field = "pos")

# Liftover the SNP positions to hg38
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
chain = import.chain(path)
seqlevelsStyle(ibd.locus.gr.hg19) = "UCSC"
ibd.locus.gr.hg38 = liftOver(ibd.locus.gr.hg19, chain)
ibd.locus.gr.hg38 = unlist(ibd.locus.gr.hg38)
#as_tibble(ibd.locus.gr.hg38)
ibd.locus.df = ibd.locus.df %>%
  dplyr::rename(pos_hg19 = pos, deLange_locus_genes = locus_genes) %>%
  left_join(as_tibble(ibd.locus.gr.hg38) %>% dplyr::select(rsid, pos_hg38 = start), by="rsid")
ibd.locus.gr.hg38 = makeGRangesFromDataFrame(ibd.locus.df %>% select(chr, pos = pos_hg38, rsid), keep.extra.columns = T, ignore.strand = T,
                                             seqinfo = NULL, seqnames.field = "chr", start.field = "pos", end.field = "pos")


###############################################################################
# Get the nearest IBD locus from the Supp Table to each Locus2Gene region
hits.df = as_tibble( distanceToNearest(l2g.gr, ibd.locus.gr.hg38, ignore.strand = T) )
# This gives us the indices of hits in the l2g table and in the ibd.locus table (GRanges)
# Since the order is the same in the data.frames we've made, we just index into those.

nearest.df = bind_cols(l2g.locus.merged.df[hits.df$queryHits,],
                       ibd.locus.df[hits.df$subjectHits,] %>% dplyr::select(-chr),
                       distance = hits.df$distance)

merged.locus.df = nearest.df %>%
  mutate(rsid_group = rsid)
merged.locus.df$rsid_group[merged.locus.df$distance > 500000] = paste0("indep", 1:(sum(merged.locus.df$distance > 500000)))

# merged.locus.df = merged.locus.df %>%
#   group_by(rsid_group) %>%
#   mutate(consensus_gene0.5 = if_else(all(top_gene == first(top_gene)) & max(top_gene_score) > 0.5, first(top_gene), ""),
#          consensus_gene0.8 = if_else(all(top_gene == first(top_gene)) & max(top_gene_score) > 0.8, first(top_gene), "")) %>%
#   dplyr::select(id, chr, l2g_pos, pos_hg19, consensus_gene0.5, consensus_gene0.8, top_gene, top_gene_score, implicated_gene, rsid, distance, locus_genes, locus_genes_ens, everything()) %>%
#   arrange(chr, l2g_pos)

merged.locus.df = merged.locus.df %>%
  select(l2g_locus_merged, chr, l2g_pos, l2g_mean_pos, ref, alt, top_gene_mean, top_gene_score_mean, locus_genes_mean, locus_genes_max, lead_p, implicated_gene, rsid, rsid_group, distance, everything())
write_tsv(merged.locus.df, path = file.path(dir, "l2g.IBD_loci.merged.annotated.tsv"))

