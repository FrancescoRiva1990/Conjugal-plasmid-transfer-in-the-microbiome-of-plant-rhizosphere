#Data Francesco DTU

setwd("C:/Users/Dell/Desktop/CNR/Francesco DTU/")

library("phyloseq")
library("dada2") #caricare il pacchetto


 path <- "C:/Users/Dell/Desktop/CNR/Francesco DTU/Trimmed/" # CHANGE ME to the directory containing the fastq files after unzipping, per dirgli che deve lavorare in quella cartella
# list.files(path) # digli di farti vedere tutti i file in quella cartella
# #
# # # Forward and reverse fastq filenames have format: SAMPLENAME_22_001.fastq and SAMPLENAME_R2_001.fastq
# fnFs <- sort(list.files(path, pattern="1_001.fastq.gz", full.names = TRUE)) # gli dico di caricare tutti i file  con il nome per esempio 1.fastq.gz
# fnRs <- sort(list.files(path, pattern="2_001.fastq.gz", full.names = TRUE))
# # # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`,1)
# sample.names
# 
plotQualityProfile(fnFs[6:10])
# #
plotQualityProfile(fnRs[6:10])
# 
# #
# # # Place filtered files in filtered/ subdirectory nuovo file Assegnare i nomi dei file per i file fastq.gz filtrati.
# filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names
# 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,  trimRight = c(20,25),
                     maxN=0, maxEE=c(2,10), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#  head(out) # ti fa vedere i primi elementi della tabella
#  View(out) # per vedere tutto ciò che c'è nella tabella
# #
errF <- learnErrors(filtFs, multithread=TRUE)
# #
errR <- learnErrors(filtRs, multithread=TRUE)
# #
save.image("dada.RData")
# #
# # plotErrors(errF, nominalQ=TRUE)
# #
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
# #
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# #
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# # # Inspect the merger data.frame from the first sample
head(mergers[[1]])
#
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#
# # Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#
seqtab.nochim2 <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 402:431]

#
sum(seqtab.nochim2)/sum(seqtab)
#
 getN <- function(x) sum(getUniques(x))
 track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim),rowSums(seqtab.nochim2))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "lenghtreduced")
# rownames(track) <- sample.names
 head(track)
# #
write.csv(seqtab.nochim2, "sequence_tab.csv")
# #
taxa <- assignTaxonomy(seqtab.nochim2, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
# #

 write.csv(taxa, "taxonomy.csv")
# # taxa2 <- addSpecies(taxa, "C:/Users/GIULIA/Desktop/schiume/silva_species_assignment_v138.1.fa.gz")
# #

   taxa.print <- taxa # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL
  head(taxa.print)
#
##load("dada.RData")

  

rm(dadaFs,dadaRs,dna,errF,errR,mergers, asvs,seqtab,seqtab.nochim)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings") 
library("ggplot2")

theme_set(theme_bw())

#seqtab.nochim2<-read.csv2("sequence_tab.csv")

rownames(seqtab.nochim2)<-DTU_phyloseq@sam_data$Sample_id
samples.out <- rownames(seqtab.nochim2)
vari<-as.data.frame(samples.out)
vari$Treatment<-c("Ecoli","Ecoli","Ecoli","Ecoli","Klebsiella","Klebsiella","Klebsiella","Klebsiella","none","none" )
vari$Replicate<-c("1", "2", "3","4","1", "2", "3","4","1","2")

rownames(vari) <- samples.out
vari<-vari[,-1]
ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE), 
               sample_data(vari), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))

names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # 
ps1 <- subset_taxa(ps0, !is.na(Genus) & !Genus %in% c("Klebsiella", "Escherichia-Shigella")) # 
ps2 <- subset_taxa(ps1, !Order %in% c("Chloroplast")) #
ps<-ps2
table(phyloseq::tax_table(ps)[, "Order"], exclude = NULL) # check orders  ##

table(phyloseq::tax_table(ps)[, "Phylum"], exclude = NULL) # check phyla
table(phyloseq::tax_table(ps)[, "Class"], exclude = NULL) # check classes
table(phyloseq::tax_table(ps)[, "Order"], exclude = NULL) # check orders, chloroplasts

table(phyloseq::tax_table(ps)[, "Genus"], exclude = NULL) # check orders, chloroplasts
table(phyloseq::tax_table(ps)[, "Phylum"], exclude = NULL) # check phyla


plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"))

library("vegan")
betabray<-vegdist(ps@otu_table,method="bray")
plot(hclust(betabray, method="average"),hang=-1, main='Bray-Curtis Bacteria', sub='', xlab='', cex=1) #plot cluster analysis of betapair

ado<-adonis2(betabray~vari$Treatment)
ado
raOTU <-  ps@otu_table
traOTU<-t(raOTU)
paOTU<-raOTU
paOTU[paOTU>0]=1


#Alpha diversity

alpha<-rowSums(paOTU) #Number of OTUs
x<-vari$Treatmet

library(GUniFrac)


library("ggplot2")
alphaV<-as.data.frame(cbind(rowSums(paOTU), vari))
colnames(alphaV)<-c("sum", "Treatment", "Replicate")
pR<-ggplot(alphaV, aes(x=factor(Treatment), y=sum))  + geom_dotplot(binaxis = "y", stackdir = "center")
pRR<-pR+labs(x = "", y=" 16S richness")+theme(text = element_text(size=15), legend.position = "bottom")+ scale_fill_manual(name = NULL, values=c('lightskyblue3', 'black', 'gold3'))
pRR
# library(GUniFrac)
# #################Rarefaction##################
# 
# rOTU<-Rarefy(raOTU,dept=min(rowSums(raOTU)))
# 
# raOTU<-as.data.frame(rOTU$otu.tab.rff)
# traOTU<-t(rOTU)
# paOTU<-raOTU
# paOTU[paOTU>0]=1



saveRDS(ps, file="ps.RDS")

table(phyloseq::tax_table(ps)[, "Phylum"], exclude = NULL) # check phyla
table(phyloseq::tax_table(ps)[, "Class"], exclude = NULL) # check classes
table(phyloseq::tax_table(ps)[, "Order"], exclude = NULL) # check orders, chloroplasts
table(phyloseq::tax_table(ps)[, "Genus"], exclude = NULL) # check orders, chloroplasts

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
# 
colNew2<-c("turquoise4","orange","gray5","lightblue4", "turquoise","gray65","gray5", "darkblue")   
ordi<-plot_ordination(ps, ord.nmds.bray, color="Treatment", title="A) nmds")+
  geom_point(size=7, alpha=0.4) +
  scale_color_manual(values=colNew2)
ordi
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Treatment", fill="Family") 


#############SUBSET only eco and kleb

psEK<-subset_samples(ps, Treatment!="none")
table(phyloseq::tax_table(psEK)[, "Phylum"], exclude = NULL) # check phyla
table(phyloseq::tax_table(psEK)[, "Class"], exclude = NULL) # check classes
table(phyloseq::tax_table(psEK)[, "Order"], exclude = NULL) # check orders, chloroplasts
table(phyloseq::tax_table(psEK)[, "Genus"], exclude = NULL) # check orders, chloroplasts

#outEK = ancombc2(data = psEK,fix_formula="Treatment")


res = outEK$res

write.csv2(res, "ancombc_res.csv")
res_global = outEK$res_global
resDB<-as.data.frame(res)
TrueF<-subset(resDB, resDB$diff_TreatmentKlebsiella=="TRUE")
tax_split<-stringr::str_split_fixed(TrueF$taxon, ":", 2)


library(randomForest)

library("randomForest")
set.seed(2014)
ntree<-100000
group<-psEK@sam_data$Treatment

sample.size.per.treatment=table(group)
predictorASV<-randomForest(y=as.factor(group), x=psEK@otu_table, ntree=ntree, importance=TRUE)
predictorASV
write.csv(predictorASV$importance, "meandecreaseASV.csv")
impASV<-as.data.frame(importance(predictorASV))
varImpPlot(predictorASV, n.var=50)

impS<-imp[order(-impASV$MeanDecreaseAccuracy),]
imp50OTU<-impS[1:50,]
#hist(imp$MeanDecreaseAccuracy)

tOTUEK<-as.data.frame(t(psEK@otu_table))
library(dplyr)
imp2<-tibble::rownames_to_column(imp50OTU)
tOTUEK<-tibble::rownames_to_column(tOTUEK)
taxaEK<-as.data.frame(psEK@tax_table)
taxaEK<-tibble::rownames_to_column(taxaEK)

impOTUTaxa<-left_join(imp2, tOTUEK, by="rowname")
impOTUTaxa2<-left_join(impOTUTaxa, taxaEK, by="rowname")

write.csv2(impOTUTaxa2, "outcome_RandomForest.csv")

dev.off() 
library("microbiome")

genus_data = aggregate_taxa(psEK, "Genus")

group<-genus_data@sam_data$Treatment
genusEK<-t(genus_data@otu_table)
predictor<-randomForest(y=as.factor(group), x=genusEK, ntree=ntree, importance=TRUE)
predictor
write.csv(predictor$importance, "meandecrease.csv")
imp<-as.data.frame(importance(predictor))
varImpPlot(predictor, n.var=50)

impS<-imp[order(-imp$MeanDecreaseAccuracy),]
imp50OTU<-impS[1:50,]
hist(imp$MeanDecreaseAccuracy)


library(dplyr)
imp2<-tibble::rownames_to_column(imp50OTU)
genusEK<-tibble::rownames_to_column(as.data.frame(genusEK))
taxaEK<-as.data.frame(genus_data@tax_table)
taxaEK<-tibble::rownames_to_column(taxaEK)

impOTUTaxa<-left_join(imp2, genusEK, by="rowname")
impOTUTaxa2<-left_join(impOTUTaxa, taxaEK, by="rowname")











