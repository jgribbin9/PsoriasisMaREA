##Project- Psoriasis RNAseq- Metabolic Profiles

##Data obtained from GEO: GSE54456
read.delim("C:\\Users\\jgrib\\Desktop\\Derm Conference Ideas\\GSE54456_RPKM_samples.txt.gz") -> RKPMCountsPsoriasis
write.csv(RKPMCountsPsoriasis, "RPKM.Psoriasis.csv")
read.delim("C:\\Users\\jgrib\\Desktop\\Derm Conference Ideas\\SraRunTable-PsoriasisPatients.txt", sep = ",") ->PsoriasisPatients
read.csv("C:\\Users\\jgrib\\Desktop\\Derm Conference Ideas\\ConversionSheetGSM.csv") -> ConversionSheetDerm
library(dplyr)
left_join(PsoriasisPatients, ConversionSheetDerm, by = c("GEO_Accession..exp." = "GSM.ID")) ->PsoriasisPatientsMeta
PsoriasisPatientsMeta%>%
  rename("GEO_Accession" = "GEO_Accession..exp.")
##In preparation for MaREA- We need seperate data for psoriasis and controls
PsoriasisPatientsMeta%>%
  filter(source_name=="Psoriasis_skin") -> PsoriasisOnlyMeta

PsoriasisPatients%>%
  filter(source_name == "normal_skin") -> NormalOnlyMeta

##MaREA
library(data.table)
marea_genes <- fread('marea_genes.csv')
marea_genes <- marea_genes %>%
  select(Symbol,
         ENS)
names(marea_genes)[names(marea_genes) == "ENS"] <- "Gene"

## Pull from data sets

rna.1.RPKM <- fread('RPKM.Psoriasis.csv')
rna.1.RPKM%>%
  rename("Symbol"="X") ->rna.1.RPKM
rna.1.RPKM <- rna.1.RPKM %>% select(!V1)
rna.1.marea <- left_join(marea_genes, rna.1.RPKM, by = 'Symbol')
write.csv(rna.1.marea, 'rna.PSORIASIS.marea.csv')


## save combined list
write.csv(rna.1.marea, 'rna.PSORIASIS.marea.csv')



## for write_TSV function
library(readr)

## remove ensembl id column to use marea tool in Galaxy
rna.1.marea.trim <- rna.1.marea %>% select(!Gene)

## write TSV file
write_tsv(rna.1.marea.trim, "rna.PSORIASIS.marea.TSV")

##Input this TSV into MaREA to get the Reaction Activity Scores

ras <- fread("C:\\Users\\jgrib\\Desktop\\Derm Conference Ideas\\RAS.Psoriasis.tabular")


ras.clus.prep <- ras %>%
  select(!Reactions)

ras.clus.prep[ras.clus.prep == 'None'] <- NA
rownames(ras.clus.prep) <- ras$Reactions
ras.clus.prep1 <- ras.clus.prep[!rowSums(is.na(ras.clus.prep)) > 0,]
ras.clustering <- t(as.matrix(ras.clus.prep1))

ras.clustering.num <- as.data.frame(sapply(as.data.frame(ras.clustering), as.numeric))
ras.rownames <- rownames(ras.clustering)

## Figuring out RAS Names

ras.ids <- ras
ras.ids[ras.ids == 'None'] <- NA
ras.ids1 <- ras.ids %>% select(!Reactions)
ras.ids1 <- as.matrix(ras.ids1)
rownames(ras.ids1) <- ras$Reactions
ras.ids2 <- ras.ids1[!rowSums(is.na(ras.ids1)) > 0,]
ras.names <- rownames(ras.ids2)

##UMAP
set.seed(123)

library(umap)
## Removing all the NA's
ras.clustering.num.no.na <- ras.clustering.num %>%
  select_if(~ !any(is.na(.)))

ras.umap <- umap(ras.clustering.num.no.na, preserve.seed = T, method = 'naive')
ras.umap.df <- as.data.frame(ras.umap[["layout"]])
colnames(ras.umap.df) <- c('X', 'Y')

library(ggplot2)
ggplot(ras.umap.df, aes(x = X, y = Y)) +
  geom_point()
## Clustering

## Determine optimal number of clusters

install.packages('factoextra')
library(factoextra)

## silhouette method
fviz_nbclust(ras.umap.df,
             FUNcluster = cluster::pam,
             method = 'silhouette',
             k.max = 6,
             print.summary = T) ##BUT USE GAP STAT INSTEAD


## kmeans method
install.packages('NbClust')
library(NbClust)
NbClust(ras.umap.df, distance = "euclidean", min.nc = 2, max.nc = 10, method = 'kmeans')


## gap_stat method
fviz_nbclust(ras.umap.df,
             FUNcluster = cluster::pam,
             method = 'gap_stat',
             k.max = 6,
             print.summary = T,
             nboot = 50)


## set up clustering

## truncate submitter id's to match the other data frame
##ras.umap.df$submitter_id <- stringr::str_trunc(ras.umap.df$submitter_id, 12, side = "right", ellipsis = '')

## add metabolic profiles
ras.umap.df$metabolic.profile <- ras.pam$clustering

## PAM

library(cluster)
ras.pam <- pam(ras.umap.df, 4, metric = "euclidean", stand = FALSE)

fviz_cluster(ras.pam, repel = T, geom = 'point', pointsize = .75)

library(dplyr)
## Add patient id's- ##I MOVED THIS STEP- Clustering is different based on if sumitter_id is in this df or not. 

patient.ids <- colnames(ras)
patient.ids <- patient.ids[2:length(patient.ids)]

ras.umap.df$submitter_id <- patient.ids


ras.umap.df$metabolic.profile <- ras.pam$clustering

ras.umap.df.clinical <- right_join(ras.umap.df, PsoriasisPatientsMeta, by = c("submitter_id" = "Patient.ID"))

ras.umap.df.clinical$metabolic.profile <- as.character(ras.umap.df.clinical$metabolic.profile)
class(ras.umap.df.clinical$metabolic.profile) ##Should say "character"

ggplot(ras.umap.df.clinical, aes(x = X, y = Y, color =metabolic.profile)) +
  geom_point()

##Eyeball investigation- Metabolic Profiles
ras.umap.df.clinical%>%
  filter(metabolic.profile == 1)

ras.umap.df.clinical%>%
  filter(metabolic.profile == 2)

ras.umap.df.clinical%>%
  filter(metabolic.profile == 3)

ras.umap.df.clinical%>%
  filter(metabolic.profile == 4)


##From eyeballing it, it looks like clusters 1,2, and 3 have the psoriatic patients (with some normal) and clusters 4 and 5 are mostly normal. Let's compare 1,2,and 3 vs. 4 and 5

ras.umap.df.clinical%>%
  filter(tissue_type == "lesional psoriatic skin") -> PsO
as.vector(PsO$submitter_id) ->PsOnames

ras.umap.df.clinical%>%
  filter(tissue_type == "normal skin") -> NormalSkin
as.vector(NormalSkin$submitter_id) ->NormalSkinnames


ras
ras$Reactions ->rasreactions


##In this step - bind columns with names of patiens to the rasreactions df. This lets you compare one subset of patients to another using the MaREA GSEA tool.
ras%>%
  select(PsO123names) -> PsO123RAS
bind_cols(rasreactions, PsO123RAS) -> PsO123RAS
PsO123RAS%>%
  rename("Reactions" = "...1") -> PsO123RAS

ras%>%
  select(PsO45names) -> PsO45RAS
bind_cols(rasreactions, PsO45RAS) -> PsO45RAS
PsO45RAS%>%
  rename("Reactions" = "...1") -> PsO45RAS

ras%>%
  select(PsOnames) -> PsORAS
bind_cols(rasreactions, PsORAS) -> PsORAS
PsORAS%>%
  rename("Reactions" = "...1") -> PsORAS

ras%>%
  select(NormalSkinnames) -> NormalSkinRAS
bind_cols(rasreactions, NormalSkinRAS) -> NormalSkinRAS
NormalSkinRAS%>%
  rename("Reactions" = "...1") -> NormalSkinRAS

##Need to export as TSV for MaREA
library(readr)
write_tsv(PsO123RAS, "PsO.123.RAS")
write_tsv(PsO45RAS, "PsO.45.RAS")
write_tsv(PsORAS, "PsO.RAS")
write_tsv(NormalSkinRAS, "NormalSkin.RAS")

ras.umap.df.clinical


{
##GEO2R Script for Microarray pateints
# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
read.delim("C:\\Users\\jgrib\\Desktop\\Derm Conference Ideas\\GSE54456_MAoverlappedsamples.txt.gz") ->Microarray.Patients
Microarray.Patients[nrow(Microarray.Patients) + 1,] <-c("M3059")
Microarray.Patients
colnames(Microarray.Patients) <-c("Patient.ID")
right_join(Microarray.Patients, ras.umap.df.clinical, by=c("Patient.ID" = "submitter_id")) ->MicroPtsClinical

MicroPtsClinical%>%
  filter(tissue_type=="lesional psoriatic skin") ->PsoMicro
as.vector(PsoMicro$Patient.ID) ->PsoMicronames


MicroPtsClinical%>%
  filter(tissue_type=="normal skin") ->NSMicro
as.vector(NSMicro$Patient.ID) ->NSMicronames

MicroPtsClinical%>%
  rename("0" = "lesional psoriatic skin")%>%
  rename("1"="normal skin") ->MicroPtsClinical.01

df <- data.frame(x = c("Bad", "Good", "Bad", "Bad", "Good"))
df$x <- as.factor(df$x)

library(tidyverse)
df <- df %>% 
  mutate(x = recode(x, 
                    "Bad" = "0", 
                    "Good" = "1"))

MicroPtsClinical$tissue_type <-as.factor(MicroPtsClinical$tissue_type)
MicroPtsClinical%>%
  mutate(tissue_type = recode(tissue_type, "lesional psoriatic skin" = "0", "normal skin"="1")) -> MicroPtsClinical.01
MicroPtsClinical.01%>%
  select(GEO_Accession..exp., tissue_type, metabolic.profile, Patient.ID)
##The above may not be necessary

################################################################
#   Differential expression analysis with limma
# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE13355", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000111111111111111111111111111111111111",
               "11111111111111111111112222222222222222222222222222",
               "222222222222222222222222222222")
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Normal No PSO","Normal Skin From PSO Patient","Psoriasis Skin from Psoriasis Patient"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE13355", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE13355", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("bottomright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE13355")
}##GEO2R

