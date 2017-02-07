library(stringr)
library(openxlsx) # exporting to excel

# NEW_IPA is where the finally used results are!

# DEFINE FUNCTIONS I USE: -----
GetGeneDirection <- function(gene, mydataframe, ... ){
  ifelse(gene %in% subset.data.frame(mydataframe, mydataframe$logFC >0)$mgi_symbol,
         1,
         ifelse(gene %in% subset.data.frame(mydataframe, mydataframe$logFC <0)$mgi_symbol,
                -1,
                0))
}

`%ni%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

cleanupIPA <- function(outprefix){
  tmp <- read.table(paste0("../raw/", paste0(outprefix, ".txt")), skip = 2, header = TRUE, sep = "\t",quote = "")
  tmp$Ntargets <- 1 + str_count(tmp$Target.molecules.in.dataset, ",")
  tmp.sig <- subset.data.frame(tmp, tmp$Ntargets >= 5 & tmp$p.value.of.overlap <= 0.01 & tmp$Molecule.Type %ni% c("chemical drug", "chemical toxicant", "chemical - endogenous non-mammalian", "chemical - kinase inhibitor", "chemical drug" ))
  write.table(unique(c(as.character(tmp.sig$Upstream.Regulator), unique(unlist(strsplit(as.character(tmp.sig$Target.molecules.in.dataset), ","))))), paste0(outprefix, "_genelist.txt"), sep = "\n", quote =FALSE, row.names = FALSE, col.names = FALSE)
  write.table(as.character(subset.data.frame(tmp.sig, tmp.sig$Predicted.Activation.State == "Activated")$Upstream.Regulator), paste0(outprefix, "_activated.txt"), sep = "\n", quote =FALSE, row.names = FALSE, col.names = FALSE)
  write.table(as.character(subset.data.frame(tmp.sig, tmp.sig$Predicted.Activation.State == "Inhibited")$Upstream.Regulator), paste0(outprefix, "_inhibited.txt"), sep = "\n", quote =FALSE, row.names = FALSE, col.names = FALSE)
  rm(tmp.sig, tmp)
}

listofobjectstoxls <- function(characterlistofobjectnames, filename = "outputworksheet", suffix = ""){
  # Test that all objects exist
  for (object in paste0(characterlistofobjectnames, suffix, sep = "")) {
    if (exists(object) == FALSE) {
      stop(paste0(paste0("The object ", object), " does not exist! Please edit your input list before continuing."))
    }
  }
  
  tmpworkbookname <- createWorkbook()
  for (object in characterlistofobjectnames) {
    addWorksheet(wb = tmpworkbookname, sheetName = object, gridLines = FALSE)
    writeData(wb = tmpworkbookname, sheet = object, x = eval(as.name(paste0(object, suffix))), borders = "none")
  }
  saveWorkbook(tmpworkbookname, file = filename, overwrite = TRUE)
  rm(tmpworkbookname, object)
}




## Load data -----


setwd("../raw/")
file_list <- substr(list.files(".",pattern = "\\.txt$"), 1, nchar(list.files(".",pattern = "\\.txt$")) - 4)
setwd("../out/")
sapply(file_list, cleanupIPA)

setwd("../tes/")
cleanupIPA("BPvsE_dw")

# QUESTION: do IPA results differ when you look at up or down individually vs up/down simultaneously?

setwd("./Final_IPA/")
 BP7vsE_logFC4 <- read.table("BP7vsE_logFC4.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
 BP7vsE_logFC4 <- subset.data.frame(BP7vsE_logFC4, BP7vsE_logFC4$p.value.of.overlap <= 0.01)
 EBvsS_2FC <- read.table("EBvsS_2FC.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
 EBvsS_2FC <- subset.data.frame(EBvsS_2FC, EBvsS_2FC$p.value.of.overlap <= 0.01)
# EBvsS_4FC_dw <- read.table("EBvsS_4FC_dw.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# EBvsS_4FC_dw <- subset.data.frame(EBvsS_4FC_dw, EBvsS_4FC_dw$p.value.of.overlap <= 0.01)
# EBvsS_4FC_up <- read.table("EBvsS_4FC_up.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# EBvsS_4FC_up <- subset.data.frame(EBvsS_4FC_up, EBvsS_4FC_up$p.value.of.overlap <= 0.01)
 #EBvsS_logFC4 <- read.table("EBvsS_logFC4.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
 #EBvsS_logFC4 <- subset.data.frame(EBvsS_logFC4, EBvsS_logFC4$p.value.of.overlap <= 0.01)
 PBvsS_2FC <- read.table("PBvsS_2FC.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
 PBvsS_2FC <- subset.data.frame(PBvsS_2FC, PBvsS_2FC$p.value.of.overlap <= 0.01)
# PBvsS_4FC_dw <- read.table("PBvsS_4FC_dw.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# PBvsS_4FC_dw <- subset.data.frame(PBvsS_4FC_dw, PBvsS_4FC_dw$p.value.of.overlap <= 0.01)
# PBvsS_4FC_up <- read.table("PBvsS_4FC_up.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# PBvsS_4FC_up <- subset.data.frame(PBvsS_4FC_up, PBvsS_4FC_up$p.value.of.overlap <= 0.01)
# PBvsS_logFC4 <- read.table("PBvsS_logFC4.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# PBvsS_logFC4 <- subset.data.frame(PBvsS_logFC4, PBvsS_logFC4$p.value.of.overlap <= 0.01)
 SP7vsE_logFC4 <- read.table("SP7vsE_logFC4.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
 SP7vsE_logFC4 <- subset.data.frame(SP7vsE_logFC4, SP7vsE_logFC4$p.value.of.overlap <= 0.01)
# P7BvsS <- read.table("P7BvsS.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# P7BvsS <- subset.data.frame(P7BvsS, P7BvsS$p.value.of.overlap <= 0.01)
# EBvsS <- read.table("EBvsS.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# EBvsS <- subset.data.frame(EBvsS, EBvsS$p.value.of.overlap <= 0.01)
# BP7vsE <- read.table("BP7vsE.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# BP7vsE <- subset.data.frame(BP7vsE, BP7vsE$p.value.of.overlap <= 0.01)
# SP7vsE <- read.table("SP7vsE.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# SP7vsE <- subset.data.frame(SP7vsE, SP7vsE$p.value.of.overlap <= 0.01)
# 
# 
# setwd("../raw/")
# BPvsE_dw <- read.table("BPvsE_dw.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# BPvsE_dw <- subset.data.frame(BPvsE_dw, BPvsE_dw$p.value.of.overlap <= 0.01)
# BPvsE_up <- read.table("BPvsE_up.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# BPvsE_up <- subset.data.frame(BPvsE_up, BPvsE_up$p.value.of.overlap <= 0.01)
# SPvsE_dw <- read.table("SPvsE_dw.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# SPvsE_dw <- subset.data.frame(SPvsE_dw, SPvsE_dw$p.value.of.overlap <= 0.01)
# SPvsE_up <- read.table("SPvsE_up.txt", sep = "\t", skip = 2, header = TRUE, quote = "")
# SPvsE_up <- subset.data.frame(SPvsE_up, SPvsE_up$p.value.of.overlap <= 0.01)
# 
# 
# table(EBvsS_logFC4$Upstream.Regulator %in% EBvsS_4FC_up$Upstream.Regulator)
# table(EBvsS_logFC4$Upstream.Regulator %in% EBvsS_4FC_dw$Upstream.Regulator)
# 
# 
# 
# library(limma)
# 
# drawvenn3D <- function(dataset1, dataset2, dataset3){
#   # What are the possible letters in the universe?
#   universe <- sort(unique(c(as.character(dataset1$Upstream.Regulator), as.character(dataset2$Upstream.Regulator), as.character(dataset3$Upstream.Regulator))))
#   
#   # Generate a matrix, with the sets in columns and possible letters on rows
#   Counts <- matrix(0, nrow=length(universe), ncol=3)
#   # Populate the said matrix
#   for (i in 1:length(universe)) {
#     Counts[i,1] <- universe[i] %in% as.character(dataset1$Upstream.Regulator)
#     Counts[i,2] <- universe[i] %in% as.character(dataset2$Upstream.Regulator)
#     Counts[i,3] <- universe[i] %in% as.character(dataset3$Upstream.Regulator)
#   }
#   
#   # Name the columns with the sample names
#   colnames(Counts) <- c(deparse(substitute(dataset1)), deparse(substitute(dataset2)), deparse(substitute(dataset3)))
#   
#   # Specify the colors for the sets
#   cols<-c("Red", "Green", "Blue")
#   vennDiagram(vennCounts(Counts), circle.col=cols)
#   rm(cols, Counts, universe)
# }
# 
# drawvenn2D <- function(dataset1, dataset2){
#   # What are the possible letters in the universe?
#   universe <- sort(unique(c(as.character(dataset1$Upstream.Regulator), as.character(dataset2$Upstream.Regulator))))
#   
#   # Generate a matrix, with the sets in columns and possible letters on rows
#   Counts <- matrix(0, nrow=length(universe), ncol=2)
#   # Populate the said matrix
#   for (i in 1:length(universe)) {
#     Counts[i,1] <- universe[i] %in% as.character(dataset1$Upstream.Regulator)
#     Counts[i,2] <- universe[i] %in% as.character(dataset2$Upstream.Regulator)
#   }
#   
#   # Name the columns with the sample names
#   colnames(Counts) <- c(deparse(substitute(dataset1)), deparse(substitute(dataset2)))
#   
#   # Specify the colors for the sets
#   cols<-c("Red", "Green")
#   vennDiagram(vennCounts(Counts), circle.col=cols)
#   rm(cols, Counts, universe)
# }
# 
# 
# drawvenn3D(EBvsS_logFC4, EBvsS_4FC_dw, EBvsS_4FC_up)
# drawvenn3D(PBvsS_logFC4, PBvsS_4FC_dw, PBvsS_4FC_up)
# drawvenn3D(BP7vsE_logFC4, BPvsE_dw, BPvsE_up)
# drawvenn3D(SP7vsE_logFC4, SPvsE_dw, SPvsE_up)
# 
# 
# 
# rm(EBvsS_4FC_dw, EBvsS_4FC_up, file_list, PBvsS_4FC_up, PBvsS_4FC_dw, SPvsE_up, SPvsE_dw, BPvsE_up, BPvsE_dw, EBvsS_logFC4, PBvsS_logFC4)
# 
# 
# 
# drawvenn2D(BP7vsE_logFC4, SP7vsE_logFC4)
# drawvenn2D(EBvsS_2FC, PBvsS_2FC)
# drawvenn2D(BP7vsE, SP7vsE)
#drawvenn2D(EBvsS, P7BvsS)

 
# Draw the venn diagram to show that E is more heterogeneous than P ----
library(VennDiagram)
# PvsE
venn.plot <- draw.pairwise.venn(209,229,
    cross.area = 104,
    category = c("BP7vsE", "SPvsE"),
  fill = c("darkblue", "darkgreen"),
  cex = 2.5,
  cat.cex = 2.5,
  euler.d =TRUE, 
);
grid.draw(venn.plot);
grid.newpage();

# BvsS
venn.plot <- draw.pairwise.venn(327,253,
                                cross.area = 125,
                                category = c("E13.5", "P7"),
                                fill = c("darkblue", "darkgreen"),
                                cex = 2.5,
                                cat.cex = 2.5,
                                euler.d =TRUE
);
grid.draw(venn.plot);
grid.newpage();

# Clean up extra column
BP7vsE_logFC4$X <- NULL
SP7vsE_logFC4$X <- NULL
EBvsS_2FC$X <- NULL
PBvsS_2FC$X <- NULL

# Now do the comparisons ----
PvsE_both <- merge.data.frame(BP7vsE_logFC4, SP7vsE_logFC4, by = "Upstream.Regulator", suffixes = c(".b", ".s"))
PvsE_onlybrain <- subset.data.frame(BP7vsE_logFC4, BP7vsE_logFC4$Upstream.Regulator %ni% PvsE_both$Upstream.Regulator)
PvsE_onlyspinalcord <- subset.data.frame(SP7vsE_logFC4, SP7vsE_logFC4$Upstream.Regulator %ni% PvsE_both$Upstream.Regulator)
BvsS_both <- merge.data.frame(EBvsS_2FC, PBvsS_2FC, by = "Upstream.Regulator", suffixes = c(".e", ".p"))
BvsS_onlyE <- subset.data.frame(EBvsS_2FC, EBvsS_2FC$Upstream.Regulator %ni% BvsS_both$Upstream.Regulator)
BvsS_onlyP <- subset.data.frame(PBvsS_2FC, PBvsS_2FC$Upstream.Regulator %ni% BvsS_both$Upstream.Regulator)

# and get the direction of the genes
# load the DE data
BPvsE <- readRDS("~/postdoc/03_OPC/160923_DGEAgo/BP7vsETotal.anno.Rds")
SPvsE <- readRDS("~/postdoc/03_OPC/160923_DGEAgo/SP7vsETotal.anno.Rds")
EBvsS <- readRDS("~/postdoc/03_OPC/160923_DGEAgo/EBvsSTotal.anno.Rds")
PBvsS <- readRDS("~/postdoc/03_OPC/160923_DGEAgo/P7BvsSTotal.anno.Rds")
# Manipulate to identify the directions 
PvsE_both$MGIgenes.b <- sapply(PvsE_both[,"Target.molecules.in.dataset.b"], function(x) str_to_title(x))
PvsE_both$MGIgenesDir.b <- sapply(sapply(PvsE_both$MGIgenes.b, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, BPvsE))
PvsE_both$MGIgenesSum.b <- sapply(PvsE_both$MGIgenesDir.b, sum)
PvsE_both$MGIgenesDir.b <- as.character(PvsE_both$MGIgenesDir.b)
PvsE_both$MGIgenes.s <- sapply(PvsE_both[,"Target.molecules.in.dataset.s"], function(x) str_to_title(x))
PvsE_both$MGIgenesDir.s <- sapply(sapply(PvsE_both$MGIgenes.s, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, SPvsE))
PvsE_both$MGIgenesSum.s <- sapply(PvsE_both$MGIgenesDir.s, sum)
PvsE_both$MGIgenesDir.s <- as.character(PvsE_both$MGIgenesDir.s)
PvsE_both$MGIgenesSum.multi <- PvsE_both$MGIgenesSum.b * PvsE_both$MGIgenesSum.s
table(PvsE_both$MGIgenesSum.multi >= 0)
#
BvsS_both$MGIgenes.e <- sapply(BvsS_both[,"Target.molecules.in.dataset.e"], function(x) str_to_title(x))
BvsS_both$MGIgenesDir.e <- sapply(sapply(BvsS_both$MGIgenes.e, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, EBvsS))
BvsS_both$MGIgenesSum.e <- sapply(BvsS_both$MGIgenesDir.e, sum)
BvsS_both$MGIgenesDir.e <- as.character(BvsS_both$MGIgenesDir.e)
BvsS_both$MGIgenes.p <- sapply(BvsS_both[,"Target.molecules.in.dataset.p"], function(x) str_to_title(x))
BvsS_both$MGIgenesDir.p <- sapply(sapply(BvsS_both$MGIgenes.p, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, PBvsS))
BvsS_both$MGIgenesSum.p <- sapply(BvsS_both$MGIgenesDir.p, sum)
BvsS_both$MGIgenesDir.p <- as.character(BvsS_both$MGIgenesDir.p)
BvsS_both$MGIgenesSum.multi <- BvsS_both$MGIgenesSum.p * BvsS_both$MGIgenesSum.e
table(BvsS_both$MGIgenesSum.multi >= 0)
# 
#
#  Now for the individual comparisons
PvsE_onlybrain$MGIgenes <- sapply(PvsE_onlybrain[,"Target.molecules.in.dataset"], function(x) str_to_title(x))
PvsE_onlybrain$MGIgenesDir <- sapply(sapply(PvsE_onlybrain$MGIgenes, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, BPvsE))
PvsE_onlybrain$MGIgenesSum <- sapply(PvsE_onlybrain$MGIgenesDir, sum)
PvsE_onlybrain$MGIgenesDir <- as.character(PvsE_onlybrain$MGIgenesDir)
#
PvsE_onlyspinalcord$MGIgenes <- sapply(PvsE_onlyspinalcord[,"Target.molecules.in.dataset"], function(x) str_to_title(x))
PvsE_onlyspinalcord$MGIgenesDir <- sapply(sapply(PvsE_onlyspinalcord$MGIgenes, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, SPvsE))
PvsE_onlyspinalcord$MGIgenesSum <- sapply(PvsE_onlyspinalcord$MGIgenesDir, sum)
PvsE_onlyspinalcord$MGIgenesDir <- as.character(PvsE_onlyspinalcord$MGIgenesDir)
#
BvsS_onlyE$MGIgenes <- sapply(BvsS_onlyE[,"Target.molecules.in.dataset"], function(x) str_to_title(x))
BvsS_onlyE$MGIgenesDir <- sapply(sapply(BvsS_onlyE$MGIgenes, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, EBvsS))
BvsS_onlyE$MGIgenesSum <- sapply(BvsS_onlyE$MGIgenesDir, sum)
BvsS_onlyE$MGIgenesDir <- as.character(BvsS_onlyE$MGIgenesDir)
#
BvsS_onlyP$MGIgenes <- sapply(BvsS_onlyP[,"Target.molecules.in.dataset"], function(x) str_to_title(x))
BvsS_onlyP$MGIgenesDir <- sapply(sapply(BvsS_onlyP$MGIgenes, function(x) strsplit(x, ",") ), function (x) GetGeneDirection(x, PBvsS))
BvsS_onlyP$MGIgenesSum <- sapply(BvsS_onlyP$MGIgenesDir, sum)
BvsS_onlyP$MGIgenesDir <- as.character(BvsS_onlyP$MGIgenesDir)

# Now export the data
listofobjectstoxls(characterlistofobjectnames = c("BvsS_both", "BvsS_onlyE", "BvsS_onlyP", "PvsE_both", "PvsE_onlybrain", "PvsE_onlyspinalcord"), filename = "IPA_Upstream.xlsx", suffix = "")


subset.data.frame(SPvsE, SPvsE$mgi_symbol == "Hoxb9")
subset.data.frame(SPvsE, SPvsE$mgi_symbol == "Lhx1")
subset.data.frame(PBvsS, PBvsS$genes == "ENSMUSG00000024985")#TCF7L2
subset.data.frame(EBvsS, EBvsS$genes == "ENSMUSG00000024985")#TCF7L2
subset.data.frame(SPvsE, SPvsE$genes == "ENSMUSG00000024985")#TCF7L2
subset.data.frame(BPvsE, BPvsE$genes == "ENSMUSG00000024985")#TCF7L2
