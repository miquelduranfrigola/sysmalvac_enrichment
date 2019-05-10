# Enrichment analysis by Camera and Tmod, as suggested by Reviewer 1

# Load libraries

library(tmod)
library(limma)
library(GSA)

# Set working directory

wd <- "/aloy/home/mduran/myscripts/sysmalvac/src_rebuttal/"

condnames = c(
  "CSP_I7_C1",
  "CSP_RTS_CONTROL",
  "iRBC_I7_C1",
  "iRBC_RTS_CONTROL",
  "malaria_CSP_CONTROL",
  "malaria_CSP_RTS",
  "malaria_iRBC_CONTROL",
  "malaria_iRBC_RTS",
  "protected_CSP_C1",
  "protected_CSP_I7",
  "protected_iRBC_C1",
  "protected_iRBC_I7"
)

genesetnames = c(
  "btm",
  "modules",
  "biocarta",
  "kegg",
  "reactome"    
)

# Variables

for (condname in condnames) {
  for (genesetname in genesetnames) {
    
OUTPATH <- paste(wd, "/sysmalvac_results/", condname, "/", genesetname, "/", sep = "")

# Read the rnk file

rnk <- read.table(paste(wd, "data/diff_expr/genes_", condname, ".rnk", sep = ""), header = FALSE, col.names = c("gene", "fc", "pval", "adj.pval"))

print(head(rnk))

# Read modules

modules <- GSA.read.gmt(paste(wd, "/data/genesets/", genesetname, ".gmt", sep = ""))
mod.names <- paste("M", 1:length(modules$geneset.names), sep = "")
mod.titles <- modules$geneset.names
mod.genes <- modules$genesets
names(mod.genes) <- mod.names

# Universes

geneset_universe <- as.character(unlist(mod.genes))


### TMOD ###

mset <- makeTmod(
  modules = data.frame(
    ID = mod.names,
    Title = mod.titles),
  modules2genes = mod.genes
)

do_tmod <- function(rnk, full_universe = TRUE) {
    
  if (full_universe) {
    tt   <- rnk
    sur <- "fulluniv"
  } else {
    tt   <- rnk[rnk$gene %in% geneset_universe,]
    sur <- "reduniv"
  }
  
  outpath <- paste(OUTPATH, "/tmod_", sur, sep = "")
  dir.create(outpath, showWarnings = FALSE)
  
  ttpos  <- tt[order(tt$fc, decreasing = TRUE),]
  ttneg  <- tt[order(tt$fc, decreasing = FALSE),]
  
  respos <- tmodCERNOtest(ttpos$gene, mset = mset, qval = 2)
  respos$direction <- "pos"
  resneg <- tmodCERNOtest(ttneg$gene, mset = mset, qval = 2)
  resneg$direction <- "neg"
  
  res <- rbind(respos, resneg)
  
  res <- res[,c("Title", "N1", "direction", "cerno", "AUC", "cES", "P.Value", "adj.P.Val")]
  colnames(res) <- c("geneset", "size", "direction", "cerno", "auc", "ces", "pval", "adj_pval")
  
  res <- res[order(res$pval),]
  res <- res[!is.na(res$auc),]
  
  write.table(res, file = paste(outpath, "/output.tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
  
  return(res)
  
}

do_tmod(rnk, TRUE)
do_tmod(rnk, FALSE)


### CAMERA (preranked...) ###

do_camera <- function(rnk, full_universe = TRUE) {
    if (full_universe) {
        tt  <- rnk
        sur <- "fulluniv"
    } else {
        tt  <- rnk[rnk$gene %in% geneset_universe, ]
        sur <- "reduniv"
    }
    outpath <- paste(OUTPATH, "/camera_", sur, sep = "")
    dir.create(outpath, showWarnings = FALSE)
    l <- tt$fc
    names(l) <- tt$gene
    res <- cameraPR(l, mod.genes)
    colnames(res) <- c("size", "direction", "pval", "fdr")
    res <- res[order(row.names(res)),]
    res$geneset <- mod.titles[order(mod.names)]
    res$direction <- ifelse(res$direction == "Up", "pos", "neg")
    res <- res[order(res$pval),]
    res <- res[!is.na(res$pval),]
    res <- res[,c("geneset", "size", "direction", "pval", "fdr")]
    write.table(res, file = paste(outpath, "/output.tsv", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
    return(res)
}

do_camera(rnk, TRUE)
do_camera(rnk, FALSE)


  }

}


