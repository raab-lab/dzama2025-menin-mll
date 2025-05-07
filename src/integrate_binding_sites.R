## Use DMSO to define regions of menin binding
## Look at differential binding of H3K4 and ASH2L at those specific sites
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(nullranges)
library(ggplot2)
library(plyranges)
chrs <- paste0("chr", c(1:22, "X", "Y"))
fdr_cut <- 0.1
lfc_cut <- 0
menin_dmso_peak_files <- list.files("data/source_data/cnr_output/peaks/",
                                    pattern = 'menin_DMSO.*narrowPeak',
                                    full.names = T)
menin_dmso_peaks <- lapply(menin_dmso_peak_files, function(file) {
  peak <- read_narrowpeaks(file) %>%
    filter(seqnames %in% chrs)
  seqlevels(peak) <- seqlevelsInUse(peak)
  return(peak)
})
menin_binding_sites <- bind_ranges(menin_dmso_peaks[[1]],
                                   menin_dmso_peaks[[2]]) %>%
  csaw::mergeWindows(ranges = ., tol = 51, max.width = 2000)
menin_binding_sites <- menin_binding_sites$regions
menin_binding_sites$peak <- 1:NROW(menin_binding_sites)
k4_diff_bind <- readRDS("data/derived_data/binding/cnr_H3K4me3_threshold5_efficiency_SNDX_vs_DMSO.rds")
menin_k4 <- join_overlap_inner(x = k4_diff_bind,
                                   y = menin_binding_sites,
                                   maxgap = 50)
write.table(menin_k4 %>% select(FDR, best.logFC) %>% as.data.frame(),
            "data/derived_data/binding/menin_H3K4me3_overlap.bed",
            row.names = F, col.names = F, quote = F, sep = "\t")


k4_loss <- k4_diff_bind %>%
  filter(best.logFC < lfc_cut & FDR < fdr_cut)
k4_menin_loss_olap <- sum(count_overlaps(k4_loss, menin_binding_sites, maxgap = 50))

## Simple test, find the null distribution for the number of overlaps between
## k4 loss and menin binding sites
## Generate segmentation
exclusion <- read_bed("/proj/jraablab/users/pkuhlers/seq_resources/hg38_cnr_exclusion.bed")
g <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
g <- keepStandardChromosomes(g, pruning.mode = "coarse")
seqlevels(g, pruning.mode = "coarse") <-
  setdiff(seqlevels(g), "chrM")
g <- sortSeqlevels(g)
g <- sort(g)
seg_cbs <-
  segmentDensity(
    g,
    n = 3,
    L_s = 1e6,
    exclude = exclusion,
    type = "cbs"
  )

## Add seqlengths to menin and K4 GRanges
seqlengths(menin_binding_sites) <-
  seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)[names(seqlengths(menin_binding_sites))]
boots <- bootRanges(menin_binding_sites,
                    blockLength = 5e5,
                    R = 30,
                    seg = seg_cbs,
                    exclude = exclusion)

boot_olap <- join_overlap_inner(x = k4_loss %>% mutate(x_id = 1:NROW(k4_loss)), y = boots, maxgap = 50) %>%
  group_by(x_id, iter) %>%
  summarize(n_overlaps = n()) %>%
  tibble::as_tibble() %>%
  tidyr::complete(iter, fill = list(n_overlaps = 0)) %>%
  group_by(iter) %>%
  summarize(sumOverlaps = sum(n_overlaps))

ggplot(boot_olap, aes(sumOverlaps)) +
  geom_histogram(bins = 6, color = 'black', fill = 'blue') +
  theme_classic()
  # geom_vline(xintercept = k4_menin_loss_olap, linetype = 'dashed')

## Are there more menin peaks at regions of k4 loss than not?
k4_no_loss <- k4_diff_bind %>%
  filter(best.logFC > 0) %>%
  plyranges::select(best.logFC, FDR)
boot_no_loss <- replicate(10,
                          slice(k4_no_loss,
                                sample.int(n(), length(k4_loss))),
                          simplify = F) %>%
  bind_ranges(.id = 'resample')
k4_regions <- bind_ranges(loss = k4_loss %>% plyranges::select(best.logFC, FDR),
                          no_loss = boot_no_loss, .id = 'origin') %>%
  mutate(resample = ifelse(is.na(resample), 0L, as.integer(resample)))

k4_loss_menin <- join_overlap_left(x = k4_regions,
                                   y = menin_binding_sites)
k4_loss_menin %>%
  group_by(origin) %>%
  summarize(peak_count = sum(!is.na(peak)) / n_distinct(resample)) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = -origin) %>%
  tidyr::pivot_wider(names_from = origin, values_from = value) %>%
  mutate(enrichment = loss / no_loss)



