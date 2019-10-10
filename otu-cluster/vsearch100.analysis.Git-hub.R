# Summarizing the OTU result for Cas-16S-seq data
# Kabin Xie
# 2019.9.5

library(ggplot2)
library(cowplot)
library(gridExtra)

# remove the # label in the first header line in all.otutab100.txt
# In this dataset, OTUs (clustered by 100%) matched  rice mito rDNA with 2 mismatches are in all.otu100.mt2.txt
# Pay attention to Col names if it start with digit.

########################################################################################
# A Input otu table and sample information
########################################################################################
otutab<-read.table(file="../vsearch.out/all.otutab100.txt", head=T, check.names=FALSE)
head(otutab)
rownames(otutab)<-otutab[,1]
OTU.mt<-read.table(file="all.otu100.mt2.txt", head=T)
table(!otutab$OTU_ID %in% OTU.mt$OTU_ID)
otu.matrix<-as.matrix(otutab[!otutab$OTU_ID %in% OTU.mt$OTU_ID,-1])  #-1, remove first name column

head(otu.matrix)
dim(otu.matrix)
dim(otutab)

mito.no<-colSums(otutab[otutab$OTU_ID %in% OTU.mt$OTU_ID, -1])
mito.ratio<-mito.no/apply(otutab[,-1],2, sum)
mito.ratio
mito_bac<-rbind(mito.no, apply(otu.matrix,2, sum))
rownames(mito_bac)<-c("Mito", "Bact")
barplot(as.matrix(mito_bac[,1:10]))
mito_bac<-rbind(mito.ratio, 1-mito.ratio)
row.names(mito_bac)<-c("Mito", "bact")

##read sample meta info
sample.info<-read.table("../vsearch_sample_info.csv", head=T, sep=",")
rownames(sample.info)<-sample.info$Seq_ID
mito_bac<-t(mito_bac)
as.character(sample.info$Seq_ID)==as.character(rownames(mito_bac)) # check the order
sample.info<-cbind(sample.info, mito_bac[rownames(sample.info),])

## read vsearch-pipline summary##
vs.summary<-read.table(file="../vs.pipline.summary100.csv", head=TRUE, sep=",")
head(vs.summary)
rownames(vs.summary)<-vs.summary$sample
sample.info<-cbind(sample.info, vs.summary[rownames(sample.info), c(2,3,4,5)])
OTU_No<-data.frame(OTU_No=colSums(otutab[-1,]>0))
 

sample.info<-cbind(sample.info, OTU_No=OTU_No[rownames(sample.info),])
sample.info<-sample.info[with(sample.info, order(gRNA, Compartment,Cas9, Cycle_No)),]
write.table(sample.info, file="sample.info.combined100.txt", sep="\t")
rm(vs.summary, OTU_No)

##########################################################
## Plot OTU numbers and Mito ratio Cas+ vs Cas-
##########################################################
#remove the negative control
sample.info=sample.info[sample.info$Seq_ID!="bS1_4",]  #remove the empty control
qplot(interaction(Cas9, Compartment, Cycle_No), derep/merged, data=sample.info, geom=c("boxplot", "jitter"))
p.reads<-qplot(interaction(Cas9, Compartment, Cycle_No), 100*derep/filtered, 
               fill = Cas9, data=sample.info)

#p.reads<-ggplot(aes(y = derep/filtered, x = Compartment, fill = Cas9),ylim(0,1), 
#                data = rbind(sample.info.root, sample.info.soil)) 
p.reads<-p.reads+geom_boxplot();
#p.reads<-p.reads+geom_jitter();
p.reads<-p.reads+theme(text=element_text(size=16))
p.reads<-p.reads+theme_grey();
p.reads
ggsave(filename = "derep_ratio.pdf", p.reads,
       width = 10, height = 10, 
       units = "cm")

sample.info.root<-sample.info[sample.info$Compartment == "Root", ]
p.OTU.No<-qplot(interaction(Cas9, gRNA,Cycle_No), OTU_No, fill = Cas9, data=sample.info.root)
p.OTU.No<-p.OTU.No+geom_boxplot()
p.OTU.No

p.ratio<-ggplot(aes(y=Mito*100, x=interaction(gRNA, Cycle_No), fill=Cas9),ylim(0,100), data=sample.info.root)
p.ratio<-p.ratio+geom_boxplot()
p.ratio<-p.ratio+theme_grey()
p.ratio
ggsave(filename="Mito_vs_Cas9.pdf", p.ratio,
       width = 10, height = 10, 
       units = "cm")
ggsave(filename="Figure5A.pdf", plot_grid(p.ratio, p.OTU.No, rel_widths = c(1, 1)),
       width = 18, height = 6, 
       units = "cm")

rm(p.ratio, p.reads, p.OTU.No)


#############################################################
#   Another boxplot to show mito-rDNA ratio
#########################################################

sample.info.root<-sample.info[sample.info$Compartment == "Root", ]
sample.info.root<-sample.info.root[with(sample.info.root, order(Cycle_No,gRNA, Cas9)),]
barplot(sample.info.root$Mito, ylim=c(-1.5,1.5))
barplot(-sample.info.root$bact, add=T, col="black")
legend(18,0.8, "Mito ratio")
legend(18,-1, "Bacterial")

###############################################################
##  B. analyze OTUs in soil sample (Fig. 4)
##
###########################################################################

sum.threshold=10  #each OTU have more then 10 reads
sample.info.soil<-sample.info[sample.info$Compartment == "Soil", ]

#split sample info according gRNA #####
lst.sample<-split(sample.info.soil, as.character(sample.info.soil$gRNA))
gRNA_id<-names(lst.sample)

# calculate the corealtion coefficient value
dim(otu.matrix)  #check the data
cor.df<-cor(otu.matrix[,colnames(otu.matrix) %in% sample.info.soil$Seq_ID])
write.table(cor.df, file="corr-R.txt")
rm(cor.df)

p<-ggplot()
otu.matrix<-as.matrix(otutab[,-1])  #Mt OTUs are included here

for(j in 1:length(gRNA_id)) { 
  subset.sample<-lst.sample[[gRNA_id[j]]]
  subset.seq<-colnames(otu.matrix) %in% subset.sample$Seq_ID  #Get sequence_ID
  subset.otutab<-subset(otu.matrix,select=subset.seq)

## filter low adundant OTUs ##
  subset.otutab<-subset.otutab[apply(subset.otutab, 1, sum) > 0, ]
  filter.abund<-apply(subset.otutab, 1,sum)>sum.threshold
  filter.freq<-apply(subset.otutab>0, 1, min)>0  # at least 1 read in 6 samples
  subset.otutab.filtered<-subset.otutab[filter.abund & filter.freq, ]
  dim(subset.otutab.filtered)
  dim(subset.otutab)
  cor(subset.otutab.filtered)
  cor(subset.otutab)
  head(subset.otutab)
## compare OTU abundances between Cas9 +/- samples ###

  subset.sample$Plant_sample<-as.character(subset.sample$Plant_sample)

# get each pair of Cas9+/- samples
  split_sample_info<-split(subset.sample, subset.sample$Plant_sample)
  t.otu.all<-data.frame(plus=as.numeric(), minus=as.numeric(),rep=as.character())
  for (i in 1:length(split_sample_info)) {
    t<-split_sample_info[[i]]
    t.otu<-subset.otutab.filtered[,colnames(subset.otutab.filtered) %in% t$Seq_ID]
    t.otu<-data.frame(t.otu)
    colnames(t.otu)<-t[colnames(t.otu),5]
    t.otu$rep<-names(split_sample_info)[i]
    t.otu.all<-rbind(t.otu.all, t.otu)
    #ggplot(x=log(plus), y=log(minus), data=data.frame(t.otu), geom=c("point"), size=I(0.5),asp = 1)
  }
  p[[j]]<-qplot(x=log(Yes), y=log(No), data=t.otu.all, 
        geom=c("point"), 
        size=I(2),
        alpha=0.5,
        asp = 1, 
        shape= factor(rep), 
        col=factor(rep)) + geom_smooth(method="lm", se=F)
  p[[j]]<-p[[j]]+labs(title =gRNA_id[j], x="log(Cas9+)", y= "log(Cas9-)")
  p[[j]]<-p[[j]]+theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5))

  p[[j]]<-p[[j]]+theme_gray()

}

#plot with bar chart
p.all<-plot_grid(p[[1]],p[[2]],p[[3]], nrow=1, labels=c("A", "B","C"),align = 'v', label_size = 20)
#ggsave(filename = "soil_otu-paired.pdf", p.all,
#       width = 15, height = 5, 
#       units = "in")

########## Get OTU observed in Cas9+ and Cas9- samples ###################

#OTUs less the threshold labeled as no detected.
count.overlap <- function(t, threshold=10) { #make sure col 1 a is Cas9 plus
  a<-length(t[t[,1]>threshold,1])
  b<-length(t[t[,2]>threshold,1])
  a_b<-length(t[t[,2]>threshold&t[,1]>threshold,1]);
  total<-a+b-a_b
  count.ov<-c(total, a-a_b, b-a_b, a_b);
  names(count.ov)<-c("Total", colnames(t),"OV")
  count.ov
}

overlap.info.all<-data.frame(
  gRNA  = as.character(),
  Sample= as.character(), 
  Cat   = as.character(),
  Num   = as.numeric()
  #Minus_sp = as.numeric(),
  #Ov    = as.numeric()
)
p.ov<-ggplot()
p.ov.all<-ggplot()
sum.threshold=10
for(j in 1:length(gRNA_id)) {
  subset.sample<-lst.sample[[gRNA_id[j]]]
  subset.seq<-colnames(otu.matrix) %in% subset.sample$Seq_ID
  subset.otutab<-subset(otu.matrix,select=subset.seq)
  ##remove zero OTUs ##
  subset.otutab<-subset.otutab[apply(subset.otutab, 1, sum) > 0, ]
  filter.abund<-apply(subset.otutab, 1,sum)>sum.threshold
  #total frequencies big than thredshold in both samples
  subset.otutab.filtered<-subset.otutab[filter.abund, ]

  subset.sample$Plant_sample<-as.character(subset.sample$Plant_sample)
  split_sample_info<-split(subset.sample, subset.sample$Plant_sample)
  overlap.info<-data.frame(
    gRNA  = as.character(),
    Sample= as.character(), 
    Cat   = as.character(),
    Num   = as.numeric()
    #Minus_sp = as.numeric(),
    #Ov    = as.numeric()
  )
  for (i in 1:length(split_sample_info)) {
    t<-split_sample_info[[i]]
    t.otu<-subset.otutab.filtered[,colnames(subset.otutab.filtered) %in% t$Seq_ID]
    t.otu<-data.frame(t.otu)
    colnames(t.otu)<-t[colnames(t.otu),5]  # sample_info rownames = t.otu colnames
    overlap<-count.overlap(t.otu, threshold = 0)
    overlap.df<-data.frame(gRNA=gRNA_id[j],Sample= names(split_sample_info)[i], Cat=names(overlap)[c(2,4,3)], Num=overlap[c(2,4,3)])
    overlap.info<-rbind(overlap.info, overlap.df,make.row.names=F)
    overlap.info.all<-rbind(overlap.info.all, overlap.df,make.row.names=F)
  }
  p.ov[[j]]<-ggplot(data=overlap.info, aes(x=Sample, y=Num, fill=Cat))+geom_bar(stat="identity", width = 0.6)
  p.ov[[j]]<-p.ov[[j]]+labs(title =gRNA_id[j], x=paste(split_sample_info[[1]][1,1],"samples"), y= "OTU count")
  p.ov[[j]]<-p.ov[[j]]+geom_text(aes(label=Num),vjust=1.0, color="black")
  p.ov[[j]]<-p.ov[[j]]+theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5))
  p.ov[[j]]<-p.ov[[j]]+theme_gray()
}  

p.all<-plot_grid(p.ov[[1]],p.ov[[2]],p.ov[[3]], 
                 p[[1]], p[[2]], p[[3]],
                 nrow=2)


##########Output the result, change the filename##########
# Fig 4
ggsave(filename = "Soil_otu100-paired-new-OV-thred-10.pdf", p.all,
       width = 15, height = 10, 
       units = "in")

write.table(overlap.info.all, file="Soil_OTU100_overlap_thred_10.txt")
rm(p.all,p.ov, p)


##########################################################
#  C. Taxanomic analysis                                    #
##########################################################
library(phyloseq)
library(gridExtra)
library(DECIPHER)
library(dada2)
library(Biostrings)
library(ShortRead)
#get OTU sequences and assign taxanomy
otu.seq<-readDNAStringSet("../vsearch.out/all.otus100.fasta", format="fasta")
table((width(otu.seq)))

#taxonomy assign taken more than 4 hours with 7 cores.
set.seed(7)
taxa <- assignTaxonomy(otu.seq, "~/mystation/microbiome/silva_nr_v132_train_set.fa.gz", multithread=TRUE, minBoot=60)
#save the workspace here, RStudio may crash here. (set n to read n seq each time, save memory)
taxa.specis<-addSpecies(taxa, "~/mystation/microbiome/silva_species_assignment_v132.fa.gz", n=1000, verbose=T)

## rowname of taxa is the OTU sequence which is named by fasta header ##
taxa.rownames.matrix<-matrix(unlist(strsplit(names(rownames(taxa)), ";")), ncol=2, byrow=T)
rownames(taxa)<-taxa.rownames.matrix[,1]

ps <- phyloseq(otu_table(otutab[,-1], taxa_are_rows=T), 
               sample_data(sample.info), 
               tax_table(taxa))  # the control bs1-4 was not read

# subset root and soil samples #
ps.root<-phyloseq(otu_table(otutab[, as.character(sample.info.root$Seq_ID)], taxa_are_rows = T),
                  sample_data(sample.info.root),
                  tax_table(taxa))
ps.root<-prune_taxa(rowSums(otu_table(ps.root))>0, ps.root )
ps.soil<-phyloseq(otu_table(otutab[,as.character(sample.info.soil$Seq_ID)],taxa_are_rows = T), 
                  sample_data(sample.info.soil),
                  tax_table(taxa))
ps.soil<-prune_taxa(rowSums(otu_table(ps.soil))>0, ps.soil )

sample_sums(prune_taxa(rowSums(otu_table(ps.soil))>10,ps.soil))

#ps<-prune_taxa(rownames(otu_table(ps)[-1,]), ps)

ps_richness<-estimate_richness(ps, measures =c("Shannon","Simpson"))
sample.info<-cbind(sample.info, estimate_richness(ps, measures =c("Shannon","Simpson"))[as.character(sample.info$Seq_ID),])

#removing Mt sequences
ps.noMt<-prune_taxa(!rownames(otu_table(ps)) %in% OTU.mt$OTU_ID, ps)
ps_richness<-estimate_richness(ps.noMt, measures =c("Shannon","Simpson"))
colnames(ps_richness)<-c("Shannon_NoMt","Simpson_NoMt")
sample.info<-cbind(sample.info, ps_richness[as.character(sample.info$Seq_ID),])

write.table(sample.info, file="sample.info.with.alpha.txt", sep="\t")

#############  ANOVA  test diversity using the all datasets.################
# 
library(agricolae)
sample.info.root<-sample.info[sample.info$Compartment == "Root",]
sample.info.root$gRNA<-relevel(sample.info.root$gRNA,ref="No")
sample.info.root$Cycle_No<-as.factor(sample.info.root$Cycle_No)
aov.test<-aov(filtered~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
outHSD<-HSD.test(aov.test, "gRNA")
outHSD
outHSD<-HSD.test(aov.test, "Cycle_No")
outHSD
TukeyHSD(aov.test)
plot(TukeyHSD(aov.test, "Cycle_No"))

# Mito relative abundance
aov.test<-aov(Mito~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD
outHSD<-HSD.test(aov.test, "Cycle_No")
outHSD

sample.info.root.8<-sample.info.root[sample.info.root$Cycle_No == 8,]
aov.test<-aov(Mito~gRNA, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

#OTU number
aov.test<-aov(OTU_No~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

# Shannon
#aov.test<-aov(Shannon~gRNA+Cycle_No, data=sample.info.root)
aov.test<-aov(Shannon_NoMt~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD
# Simpson
aov.test<-aov(Simpson~gRNA+Cycle_No, data=sample.info.root)
summary(aov.test)
TukeyHSD(aov.test)
outHSD<-HSD.test(aov.test,"gRNA")
outHSD

estimate_richness(ps.root)
estimate_richness(ps.soil)

prune_taxa(data.frame(tax_table(ps.root))$Family != "Mitochondria", ps.root)
ps.mito<-prune_taxa(data.frame(tax_table(ps.root))$Family=="Mitochondria",ps.root)
table(rowSums(otu_table(ps.mito))>2)

####################################################
# D. microbiome diversities 
#    Fig 5C
####################################################

## Root samples
sample.info.root<-sample.info[sample.info$Compartment == "Root", ]
t<-colSums(otu_table(prune_taxa(rowSums(otu_table(ps.root))>10,ps.root))>0)
sample.info.root<-cbind(sample.info.root, Effective_OTU=t[rownames(sample.info.root)])
p.OTU.No<-qplot(interaction(Cas9, gRNA,Cycle_No), Effective_OTU, fill = Cas9, data=sample.info.root)
p.OTU.No<-p.OTU.No+geom_boxplot()
p.OTU.No
rm(p.OTU.No)
sample.info.root<-sample.info.root[with(sample.info.root, order(Cycle_No,gRNA, Cas9)),]

#p.richness<-qplot(interaction(Cas9, gRNA, Cycle_No),Shannon, fill = Cas9, data=sample.info.root)
p.richness<-qplot(interaction(Cas9, gRNA, Cycle_No),Shannon_NoMt, fill = Cas9, data=sample.info.root)
p.richness<-p.richness+geom_boxplot()
p.richness<-qplot(interaction(Cas9, gRNA, Cycle_No),Simpson, fill = Cas9, data=sample.info.root)
p.richness

# generate rarefy plot #
lst.sample<-split(sample.info.root, sample.info.root$gRNA)
gRNA_id<-names(lst.sample)

root.otu.sum<-sample_sums(ps.root)

######### Alpha diversity  ################

 #rarefy to even depth and repeat xx times. 
 # makesure the mito OTU is included in analyze.
 ps.root.noMt<-prune_taxa(!rownames(otu_table(ps.root)) %in% OTU.mt$OTU_ID, ps.root)
 sample_size<-min(sample_sums(ps.root));
 set.seed(7)
 
 rarefy.otu<-otu_table(rarefy_even_depth(ps.root, sample.size=sample_size, replace=T))
 rarefy.otu
 rarefy.rich<-estimate_richness(rarefy.otu, measures = c("Shannon","Simpson"))
 rarefy.otu.NoMt<-otu_table(rarefy_even_depth(ps.root.noMt, sample.size = min(sample_sums(ps.root.noMt)),replace =T))
 rarefy.rich.NoMt<-estimate_richness(rarefy.otu.NoMt, measures = c("Shannon","Simpson"))
 
 sample.info.root$gRNA<-relevel(sample.info.root$gRNA, ref="No")
 ggplot(data=sample.info.root, aes(x=gRNA, y=Shannon, color=gRNA))+geom_boxplot()
 qplot(interaction(gRNA, Cycle_No), Shannon_NoMt, data=sample.info.root, color=gRNA)+geom_boxplot()

 # Generate rarefy curve 

 #generate ps.root and ps.root.Mt rarefy plot indepdently
 sample_size<-min(sample_sums(ps.root.noMt));  #change to ps.root/ps.root.noMt
 rarefy.observed<-data.frame()
 rarefy.otu<-data.frame()
 for( i in seq(1, sample_size, by = 500)) {
  rarefy.otu<-otu_table(rarefy_even_depth(prune_taxa(rowSums(otu_table(ps.root))>10,ps.root), sample.size=i, replace=T))  
  # change to ps.root, prune_taxa(rowSums(otu_table(ps.root))>10,ps.root), ps.root.noMt
  
  t<-data.frame(Rarefy=i,OTU.No=colSums(rarefy.otu>0))
  t<-cbind(sample.info.root[rownames(t), c("Plant_sample","gRNA", "Cas9", "Cycle_No")], t)
  rarefy.observed<-rbind(rarefy.observed, t)
  rm(t)
 }
 
 tail(rarefy.observed)
 
 rarefy.observed$gRNA<-relevel(rarefy.observed$gRNA, ref="No")
 t<-rarefy.observed[rarefy.observed$Cycle_No=="8",]
 p.rarefy<-ggplot(aes(x=Rarefy/1000, y=OTU.No, colour=gRNA, Cycle_No),data=t)
 p.rarefy<-p.rarefy+geom_smooth(fill="white",weight=0.5)+facet_wrap(~Plant_sample)
 p.rarefy<-p.rarefy+labs(x="Sequences (x 1000)", y="OTU Number")
 p.rarefy
 ggsave(filename = "Rare_Curve_8cycles_OTU100_threshold10.pdf", p.rarefy,
        width = 15, height = 5, 
        units = "in")
rm(t, p.rarefy, rarefy.ovserved, rarefy.otu)

ps.soil<-prune_taxa(rowSums(otu_table(ps.soil))>0, ps.soil)
soil_richness<-estimate_richness(ps.soil)
sample.info.soil<-cbind(sample.info.soil, soil_richness[rownames(sample.info.soil),])
p.richness<-qplot(interaction(Cas9, gRNA), Chao1, fill = Cas9, data=sample.info.soil)
p.richness<-p.richness+geom_boxplot()

## Analyze beta diversities
ps.root.noMt<-prune_taxa(rowSums(otu_table(ps.root.noMt))>10,ps.root.noMt)
pslog<-transform_sample_counts(ps.root.noMt, function(x) log(1+x))
out.log<-ordinate(pslog, method="NMDS", distance="bray")
p.ord<-plot_ordination(pslog, out.log, color="gRNA", shape="Plant_sample")+geom_point(size=I(2))
ggsave(filename = "plot_Ordination_NMD_OTU100_10.pdf", p.ord, width=3, height=3, units= "in")

pslog.soil<-transform_sample_counts(prune_taxa(rowSums(otu_table(ps.soil))>10, ps.soil), function(x) log(1+x))
out.log.soil<-ordinate(pslog.soil, method="NMDS", distance="bray")
p.ord.soil<-plot_ordination(pslog.soil, out.log.soil, color="Cas9", shape="gRNA")+geom_point(size=I(2))

ps.root.phyla<-tax_glom(ps.root.noMt, taxrank = "Phylum")
ps.root.phyla<-transform_sample_counts(ps.root.phyla, function(x) {x/sum(x)} )
ps.root.phyla<-psmelt(ps.root.phyla)

ps.root.phyla<-ps.root.phyla[ps.root.phyla$Abundance>0.01,]
table(as.character(ps.root.phyla$Phylum))
summary(ps.root.phyla$Abundance)

aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Proteobacteria",])
summary(aov.test)
aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Bacteroidetes",])
summary(aov.test)
aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Firmicutes",])
summary(aov.test)
aov.test<-aov(Abundance~gRNA+Cycle_No, data=ps.root.phyla[ps.root.phyla$Phylum=="Actinobacteria",])
summary(aov.test)

write.table(file="Root.Phylum.comparison.txt",ps.root.phyla, sep="\t")

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

#Figure 5e
p.phylum<-ggplot(data=ps.root.phyla,aes(x=gRNA, y=Abundance, colour=Phylum, fill=Phylum)) + 
  geom_bar(stat = "identity")+facet_wrap(~Cycle_No + Plant_sample)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="Cas9",y="Abundance (>1%)")

ggsave(filename = "Relative_abundance_Phyl_otu100_10.pdf", p.phylum,
       width = 5, height = 10, 
       units = "in")
rm(ps.root.phyla, p.phylum,phylum_colors)


###########################################################################
#  
#    E. DESeq2 analysis 2019.8.13
###########################################################################
library(DESeq2)

ps.soil.deseq2<-phyloseq_to_deseq2(ps.soil, ~ Cas9+gRNA)
ps.soil.deseq2
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

Deseq.analyze= function(dds, sig.alpha=0.05) {
  geoMeans<- apply(counts(dds),1, gm_mean)
  dds<- estimateSizeFactors(dds, geoMeans = geoMeans)
  dds<- DESeq(dds, fitType="local")
  res<-results(dds)
  #summary(res)
  sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
  sigtab = cbind(gRNA=i, OTU=rownames(sigtab), as(sigtab, "data.frame"), as(tax_table(ps.soil)[rownames(sigtab), ], "matrix"))
  sigtab
}

# Analyzing the differential OTUs in soil samples
# For each gRNA, three biological repeats were analyzed.
   sig.alpha=0.01  # set significant p value
   ps.soil
   prune_taxa(rowSums(otu_table(ps.soil))>0,ps.soil)
   sample_data(ps.soil)
   gRNA_id=levels(sample_data(ps.soil)$gRNA)
   Dif.OTU<-data.frame();
for(i in gRNA_id) {
   dds<-prune_samples(sample_data(ps.soil)$gRNA == i,ps.soil)
   keep <- rowSums(otu_table(dds)) >= 10  #set the threshold to filter OTUs
   table(keep)
   dds <- prune_taxa(keep,dds)
   sample_data(dds)$Cas9<-relevel(sample_data(dds)$Cas9,ref="No")
   ps.dds<-phyloseq_to_deseq2(dds, ~Cas9)
   geoMeans<- apply(counts(ps.dds),1, gm_mean)
   dds<- estimateSizeFactors(ps.dds, geoMeans = geoMeans)
   dds<- DESeq(dds, test="Wald",fitType="parametric")
   res<-results(dds)
   #summary(res)
   #plotMA(res)
   #sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
   sigtab = res[(res$pvalue < sig.alpha & !is.na(res$pvalue)), ]
   if(dim(DataFrame(sigtab))[1]>0) {
      sigtab = cbind(gRNA=i, OTU=rownames(sigtab), as(sigtab, "data.frame"), as(tax_table(ps.soil)[rownames(sigtab), ], "matrix"))
      sigtab
      #sigtab = Deseq.analyze(dds, sig.alpha)
      Dif.OTU<-rbind(Dif.OTU, sigtab, make.row.names=F)
   }
}
Dif.OTU
OTU.id<-matrix(unlist(strsplit(names(otu.seq),";size=",fixed=TRUE)), ncol=2, byrow=T)
names(otu.seq)<-OTU.id[,1]
Dif.OTU$Sequence<-as.character(otu.seq[Dif.OTU$OTU])

table(Dif.OTU$OTU)
write.table(file="Dif.OTU100.soil.threshold10.txt", Dif.OTU, sep="\t")
p<-ggplot(aes(y=Mito*100, x=interaction(gRNA, Cas9), fill=Cas9),ylim(0,100), data=sample.info.soil)
p<-p+geom_boxplot()
p<-p+theme_grey()
p
ggsave(filename = "Mito_soil_frequency.pdf", p,
       width = 15, height = 10, 
       units = "in")

rm(p,Dif.OTU, OTU.id)


##Analyzing differential OTU in root samples

# 1. Analyze differential OTU between No vs gRNA + cycles = 8.   ps.root
   sample_data(ps.root)
Deseq.analyze= function(dds, sig.alpha=0.05) {
     geoMeans<- apply(counts(dds),1, gm_mean)
     dds<- estimateSizeFactors(dds, geoMeans = geoMeans)
     dds<- DESeq(dds, fitType="local")
     res<-results(dds)
     #summary(res)
     sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
     sigtab = cbind(gRNA=i, OTU=rownames(sigtab), as(sigtab, "data.frame"), as(tax_table(ps.root)[rownames(sigtab), ], "matrix"))
     sigtab
}   

   dds<-prune_samples(sample_data(ps.root.noMt)$Cycle_No == 8,ps.root.noMt)
   dds<-prune_samples(sample_data(dds)$gRNA %in% c("No", "gRNA1196"), dds)
   dds<-prune_taxa(rownames(otu_table(dds)[-1,]), dds) # remove mito sequences
   dds
   keep <- rowSums(otu_table(dds)) > 10
   table(keep)
   dds <- prune_taxa(keep,dds)
   dds<-phyloseq_to_deseq2(dds, ~gRNA)
   sigtab<-Deseq.analyze(dds,sig.alpha = 0.05)

   rownames(sigtab) %in% OTU.mt$OTU_ID
   
   otu_table(prune_taxa(rownames(otu_table(ps.root)) %in%  rownames(sigtab),ps.root))
   tax_table(ps.root)[rownames(otu_table(ps.root)) %in%  rownames(sigtab),]
   writeXStringSet(otu.seq[names(otu.seq) %in% rownames(sigtab)], "Root.1196.dif.otu100.fasta")
   
 
# 2. Analyze differential OTUs between Cas9+ vs - per gRNA for each sample
#  results from 8 and 15 cycles were used as technicalrepeat.
# Fig. 5d

sample_id<-levels(sample_data(ps.root)$Plant_sample)
p_FCvsP<-ggplot();
p_FCvsPadj<-ggplot();
Dif.OTU<-data.frame() # output the result in this dataframe.
n=0; sig.alpha=0.01 # set to 0.05 if use padj
for(plant_sample in sample_id){  
   ps.sel<-prune_samples(sample_data(ps.root.noMt)$Plant_sample == plant_sample,ps.root.noMt)
   for (i in 1:(length(gRNA_id))) {
     ps.sel.gRNA<-prune_samples(sample_data(ps.sel)$gRNA %in% c("No", gRNA_id[i]), ps.sel)
      # choose one prune method
     ps.sel.gRNA<-prune_taxa(data.frame(tax_table(ps.sel.gRNA))$Family != "Mitochondria", ps.sel.gRNA)
     ps.sel.gRNA
     keep <- rowSums(otu_table(ps.sel.gRNA)) > 0
     table(keep)
     ps.sel.gRNA <- prune_taxa(keep,ps.sel.gRNA)
     sample_data(ps.sel.gRNA)$gRNA<-relevel(sample_data(ps.sel.gRNA)$gRNA,ref="No")
     dds<-phyloseq_to_deseq2(ps.sel.gRNA, ~gRNA)
     geoMeans<- apply(counts(dds),1, gm_mean)
     dds<- estimateSizeFactors(dds, geoMeans = geoMeans)
     dds<- DESeq(dds, fitType="parametric", test="Wald")
     res<-results(dds)
     #res["OTU_446",]
     #sigtab = res[(res$padj < sig.alpha & !is.na(res$padj)), ]
     sigtab = res[(res$pvalue<sig.alpha & ! is.na(res$pvalue)),]
     sigtab = cbind(gRNA=gRNA_id[i], OTU=rownames(sigtab), as(sigtab, "data.frame"), as(tax_table(ps.root)[rownames(sigtab), ], "matrix"))
     sigtab$Plant_sample=plant_sample     
     Dif.OTU<-rbind(Dif.OTU,sigtab)
   #summary(res)
   #plotMA(res)
   #col=res$padj>0
   p_FCvsP[[i+n]]<-ggplot(data=data.frame(res), aes(x=log2FoldChange, y=-log10(pvalue)))  
   p_FCvsP[[i+n]]<-p_FCvsP[[i+n]]+geom_point(aes(colour=pvalue>0.01),size=I(2),alpha=0.7)
   p_FCvsP[[i+n]]<-p_FCvsP[[i+n]]+geom_vline(xintercept=0, color= "grey")
   p_FCvsP[[i+n]]<-p_FCvsP[[i+n]]+labs(title = paste(plant_sample,"\n",  gRNA_id[i], " vs CK"))
   p_FCvsPadj[[i+n]]<-ggplot(data=data.frame(res), aes(x=log2FoldChange,y=-log10(padj)))
   p_FCvsPadj[[i+n]]<- p_FCvsPadj[[i+n]]+geom_point(aes(colour=padj>0.01),size=I(2),alpha=0.7)
   p_FCvsPadj[[i+n]]<- p_FCvsPadj[[i+n]]+geom_vline(xintercept=0, color= "grey")
   p_FCvsPadj[[i+n]]<- p_FCvsPadj[[i+n]]+labs(title = paste(plant_sample,"\n",gRNA_id[i], " vs CK "),hjust=0.5)
   }
   n=n+3
}    

Dif.OTU$Sequence<-as.character(otu.seq[Dif.OTU$OTU])
Dif.OTU[Dif.OTU$log2FoldChange< 0,]
#write.table(Dif.OTU, file="Dif_OTU_Root.txt", sep = "\t")   
#write.table(Dif.OTU, file="Dif_OTU_Root_pval0.01.txt", sep = "\t")  

# Get the abundance of differential OTUs LFC<0
unique(Dif.OTU$OTU)
unique(as.character(Dif.OTU$OTU[Dif.OTU$Family=="Mitochondria"]))
writeXStringSet(otu.seq[names(otu.seq) %in% unique(as.character(Dif.OTU$OTU[Dif.OTU$Family=="Mitochondria"]))], "Root.dif.Decrease.otu100.fasta")
writeXStringSet(otu.seq["OTU_1"], "otu100.otu_1.fasta")
otu_table(ps.root)[unique(as.character(Dif.OTU$OTU))]
write.table(otu_table(ps.root)[unique(as.character(Dif.OTU$OTU))], file="Dif_OTU_NoMt_rawReads.txt", sep="\t")
otu_table(ps.root)[unique(as.character(Dif.OTU$OTU))]/sample_sums(ps.root.noMt)
summary(otu_table(ps.root)[,1]/sample_sums(ps.root.noMt)[1])
write.table(otu_table(ps.root)[unique(as.character(Dif.OTU$OTU))]/sample_sums(ps.root.noMt), file="Dif_OTU_NoMt_RA.txt", sep="\t")

t<-as.character(unique(Dif.OTU$OTU[Dif.OTU$log2FoldChange>0]))
t<-otu.matrix[t,]
colSums(t)/colSums(otu.matrix)  #relative abundance of all differential OTU logFC>1.
cbind(sample.info[rownames(t),1:6],t)

write.table(Dif.OTU, file="Dif_OTU_Root_pval0.01-NoMt.txt", sep = "\t")    #-w-ChrMr to store resultes including Mt OTU in analysis


p.all<-plot_grid(p_FCvsP[[1]],p_FCvsP[[4]],p_FCvsP[[7]],
                 p_FCvsP[[2]],p_FCvsP[[5]],p_FCvsP[[8]],
                 p_FCvsP[[3]],p_FCvsP[[6]],p_FCvsP[[9]],nrow=3)

ggsave(filename = "Dif_OTU100_root_LFCvsPvalue_NoMt.pdf", p.all,
       width = 15, height = 10, 
       units = "in")

p.all<-plot_grid(p_FCvsPadj[[1]],p_FCvsPadj[[4]],p_FCvsPadj[[7]],
                 p_FCvsPadj[[2]],p_FCvsPadj[[5]],p_FCvsPadj[[8]],
                 p_FCvsPadj[[3]],p_FCvsPadj[[6]],p_FCvsPadj[[9]],nrow=3)

ggsave(filename = "Dif_OTU_root_LFCvsPadj_NoMt.pdf", p.all,
       width = 15, height = 10, 
       units = "in")


   
