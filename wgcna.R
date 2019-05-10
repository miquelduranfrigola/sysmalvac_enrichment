# Import library

library(DCGL)

# Variables

setwd("/aloy/home/mduran/myscripts/sysmalvac/src_rebuttal/")

FILESPATH <- "data/individuals/"

comparisons <- list(
  genes_CSP_I7_C1 = c("malaria-protected_C1_CSP", "malaria-protected_I7_CSP"),
  genes_iRBC_RTS_CONTROL = c("malaria-protected_RTSS_iRBC", "malaria-protected_comparator_iRBC"),
  genes_malaria_iRBC_CONTROL = c("protected_comparator_iRBC", "malaria_comparator_iRBC"),
  genes_protected_CSP_I7 = c("protected_I7_CSP", "malaria_I7_CSP"),
  genes_CSP_RTS_CONTROL = c("malaria-protected_RTSS_CSP", "malaria-protected_comparator_CSP"),
  genes_malaria_CSP_CONTROL = c("protected_comparator_CSP", "malaria_comparator_CSP"),
  genes_malaria_iRBC_RTS = c("protected_RTSS_iRBC", "malaria_RTSS_iRBC"),
  genes_protected_iRBC_C1 = c("protected_C1_iRBC", "malaria_C1_iRBC"),
  genes_iRBC_I7_C1 = c("malaria-protected_C1_iRBC", "malaria-protected_I7_iRBC"),
  genes_malaria_CSP_RTS = c("protected_RTSS_CSP", "malaria_RTSS_CSP"),
  genes_protected_CSP_C1 = c("protected_C1_CSP", "malaria_C1_CSP"),
  genes_protected_iRBC_I7 = c("protected_I7_iRBC", "malaria_I7_iRBC")
  )


# Main function

do_wgcna <- function(comp, expr1, expr2) {

expr1 <- read.table(paste(FILESPATH, "/", expr1, ".tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)
expr2 <- read.table(paste(FILESPATH, "/", expr2, ".tsv", sep = ""), header = TRUE, sep = "\t", row.names = 1)

genes  <- sort(WGCNA(expr1, expr2, power = 6), decreasing = TRUE)

filename = paste("data/wgcna_cor/", comp, ".rnk", sep = "")
write.table(genes, file = filename, sep = "\t", quote = FALSE, col.names = FALSE)

}

for (comp in names(comparisons)) {
  expr1 = comparisons[[comp]][1]
  expr2 = comparisons[[comp]][2]
  do_wgcna(comp, expr1, expr2)
}


