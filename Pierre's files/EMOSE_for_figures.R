setwd ("~/Documents/EMOSE2017/R")
#--------------------------------------------------------------
#Box plot shannon

divtotal16S <- read.table("16S_div_boxplot_day1.txt", header=T, sep="\t", check.names=F,row.names=1)
divtotal16S$method_size_f <- factor(divtotal16S$method_size,levels = c("metab_S", "mitag_S", "metab_L", "mitag_L"))
allGroupsColors<- c( "#F8766D", "#00BFC4", "#7CAE00", "#C77CFF")
ggplot(divtotal16S, aes(x=method_size_f, y=Chao_1))  + theme_bw()+
  geom_boxplot( alpha = 0.5)+
  scale_colour_hue(l=50) +
  scale_color_manual(values = allGroupsColors)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_jitter(aes(color = Size),position=position_jitter(0.2), size=5,alpha=0.7)+
  labs(x = "Method", y="Shannon index")+
  scale_x_discrete(breaks=c("metab_S", "mitag_S", "metab_L", "mitag_L"),
                     labels=c("metab", "mitag", "metab", "mitag"))

divtotal16S <- read.table("mitag_eukvsprok_div_boxplot_day1.txt", header=T, sep="\t", check.names=F)
divtotal16S <- read.table("metab_eukvsprok_div_boxplot_day1.txt", header=T, sep="\t", check.names=F)
divtotal16S$method_size_f <- factor(divtotal16S$method_size,levels = c( "prok_S","prok_L", "euk_S","euk_L" ))
allGroupsColors<- c( "#F8766D", "#00BFC4", "#7CAE00", "#C77CFF")
ggplot(divtotal16S, aes(x=method_size_f, y=Chao_1))  + theme_bw()+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))+
  geom_boxplot( alpha = 0.5)+
  scale_colour_hue(l=50) +
  scale_color_manual(values = allGroupsColors)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_jitter(aes(color = Size),position=position_jitter(0.2), size=7,alpha=0.7)+
  labs( x="",y="Richness (chao-1)")+
  scale_x_discrete(breaks=c("prok_S","prok_L", "euk_S","euk_L" ),
                   labels=c( "Prokaryote small",  "Prokaryote large","Eukaryote small","Eukaryote large"))



--------------------------------------------------------------
  #MDS
  
library("phyloseq")

#dada
otutabledada <- read.table("16S_dada_formds_day1.txt", header=T, sep="\t", check.names=F,row.names=1)
metadada <- read.table("16Sdada_div_day1.txt", header=T, sep="\t", check.names=F,row.names=1)
otutabledada<-t(otutabledada)
otutabledada<-decostand(otutabledada, method="hellinger")
otutabledada<-t(otutabledada)
OTUTabdada<-as.matrix(otutabledada)  
OTUdada = otu_table(OTUTabdada, taxa_are_rows = TRUE)
sampledatadada=sample_data(metadada)
physeqdadan = phyloseq(OTUdada, sampledatadada)
p4title = "metab"
allGroupsColors<- c( "#F8766D", "#00BFC4", "#7CAE00", "#C77CFF")
dada.ord = ordinate(physeqdadan, method = "MDS", distance = "bray")
pwdada=plot_ordination(physeqdadan, dada.ord, "samples", color = "Size",  title = p4title)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_colour_hue(l=50)+
  scale_color_manual(values = allGroupsColors)+
  geom_point(shape=19, size=5,alpha=1/2)
pwdada$layers <- pwdada$layers[-1]
pwdada
--------------------------------------------------------------
  
  #betadisper 

otutabledada<-read.table("16S_dada_formds_day1.txt", header=T, sep="\t", check.names=F,row.names=1)
otutabledadabis <- t(otutabledada)
otutabledadabis<-decostand(otutabledadabis, method="hellinger")
comm.bc.datadadabis<-vegdist(otutabledadabis, method="bray")

data16S_day1 <- read.table("16Sdada_div_day1.txt", header=T, sep="\t", check.names=F,row.names=1)
moddada <- betadisper(comm.bc.datadadabis, data16S_day1$Size2)
moddada
permutest(moddada, control = permControl(nperm = 100))
anova(moddada)
plot(moddada)
boxplot(moddada,xlab="Size fractions",ylab="Betadispersion",cex.lab=1.4,cex.axis=1.4)
plot(TukeyHSD(moddada))

#--------------------------------------------------------------
  #-----dada 16 S cluster
  #annotated dendrogram
otutabledada<-read.table("16S_dada_formds_day1_02.txt", header=T, sep="\t", check.names=F,row.names=1)
otutabledada<-read.table("16dada_02_formds.txt", header=T, sep="\t", check.names=F,row.names=1)
otutabledadabis <- t(otutabledada)
otutabledadabis<-decostand(otutabledadabis, method="hellinger")
comm.bc.databis<-vegdist(otutabledadabis, method="bray")
#Communaut??s du cluster en utilisant un algorithme de lien moyen
comm.cluster.databis<-hclust(comm.bc.databis,method="average")
# load packages
library("ggtree")
library("ape")
dend <- as.phylo(comm.cluster.databis)
info <- read.csv(paste0("tip_data_02.csv"))
genotype <- read.table("genotype2_02_prc-euk.txt", header=T, sep="\t", check.names=F,row.names=1)
genotype <- read.table("genotype2_02_vol.txt", header=T, sep="\t", check.names=F,row.names=1)
#genotype <- read.table("genotype2_02_prc-euk.txt", header=T, sep="\t", check.names=F,row.names=1)
p <- ggtree(dend) %<+% info + xlim(0, 0.5)+ theme_tree2()
p2 <- p + geom_tiplab(offset = .005,size=3, color="black") +
  geom_tippoint(aes(color=Size),size=4,alpha=0.5) 

p2

gheatmap(p2, genotype, offset=0.1, width=0.2, low="white", high="black", colnames_position = "top",  color="black",font.size=3)
par(new=TRUE)
gheatmap(p2, genotype, offset=0.1, width=0.2,  colnames_position = "top",  color="black",font.size=3)+
  scale_fill_manual(breaks=c("V_10-100", "V_100-120", "V_2.5-10", "V_500+"), 
                    values=c("steelblue1", "steelblue4", "slategray1", "black"), name="Filter")

#-----mitag 16 S cluster------------
otutabledada<-read.table("mitag_otutableformds_day1-02.txt", header=T, sep="\t", check.names=F,row.names=1)
otutabledadabis <- t(otutabledada)
comm.bc.databis<-vegdist(otutabledadabis, method="bray")
#Communaut??s du cluster en utilisant un algorithme de lien moyen
comm.cluster.databis<-hclust(comm.bc.databis,method="average")
dend <- as.phylo(comm.cluster.databis)
info <- read.csv(paste0("tip_data_02_mitag.csv"))
genotype <- read.table("genotype2_02_metag_vol.txt", header=T, sep="\t", check.names=F,row.names=1)
p <- ggtree(dend) %<+% info + xlim(0, 0.5)+ theme_tree2()
p2 <- p + geom_tiplab(offset = .005,size=3, color="black") +
  geom_tippoint(aes(color=size),size=4,alpha=0.5) 
gheatmap(p2, genotype, offset=0.04, width=0.1, low="white", high="black", colnames_position = "top", color="black", font.size=3)+
scale_fill_manual(breaks=c("1-2.5","10", "60-100",  "500+"), 
                  values=c("slategray1","steelblue1", "black" ,"steelblue4"), name="Filter")

#----metab18S----
otutable18S <- read.table("metab18S_formds_eukonly_02.txt", header=T, sep="\t", check.names=F,row.names=1)
otutable18Sbis <- t(otutable18S)
#otutablemitagbis<-decostand(otutablemitag, method="hellinger")
otutable18Sbis<-decostand(otutable18Sbis, method="hellinger")
comm.bc.otutable18Sbis<-vegdist(otutable18Sbis, method="bray")
#Communaut??s du cluster en utilisant un algorithme de lien moyen
comm.cluster.otutable18Sbis<-hclust(comm.bc.otutable18Sbis,method="average")
dend <- as.phylo(comm.cluster.otutable18Sbis)
info <- read.csv(paste0("tip_data_02_18S.csv"))
genotype <- read.table("genotype_02_18S_artro.txt", header=T, sep="\t", check.names=F,row.names=1)
p <- ggtree(dend) %<+% info + xlim(0, 0.5)+ theme_tree2()
p2 <- p + geom_tiplab(offset = .005,size=3, color="black") +
  geom_tippoint(aes(color=size),size=4,alpha=0.5)
gheatmap(p2, genotype, offset=0.1, width=0.1, low="white", high="black", colnames_position = "top",  color="black",font.size=3)

comm.cluster.otutable18Sbis<-hclust(comm.bc.otutable18Sbis,method="average")
dend <- as.phylo(comm.cluster.otutable18Sbis)
info <- read.csv(paste0("tip_data_02_18S.csv"))
genotype <- read.table("genotype_02_18S_vol.txt", header=T, sep="\t", check.names=F,row.names=1)
p <- ggtree(dend) %<+% info + xlim(0, 0.5)+ theme_tree2()
p2 <- p + geom_tiplab(offset = .005,size=3, color="black") +
  geom_tippoint(aes(color=size),size=4,alpha=0.5)
gheatmap(p2, genotype, offset=0.1, width=0.1, low="white", high="black", colnames_position = "top",  color="black",font.size=3)+
  scale_fill_manual(breaks=c("V_2.5","V_10", "V_60-100",  "V_500+"), 
                    values=c("steelblue1","slategray1", "black" ,"steelblue4"), name="Filter")



#----mitagb18S----
otutablemitag <- read.table("mitag_euk_02_formds.txt", header=T, sep="\t", check.names=F,row.names=1)
otutablemitagbis <- t(otutablemitag)
otutablemitagbis<-decostand(otutablemitagbis, method="hellinger")
comm.bc.mitagdatabis<-vegdist(otutablemitagbis, method="bray")
#Communaut??s du cluster en utilisant un algorithme de lien moyen
comm.cluster.mitagdatabis<-hclust(comm.bc.mitagdatabis,method="average")
dend <- as.phylo(comm.cluster.mitagdatabis)
info <- read.csv(paste0("tip_data_02_mitag.csv"))
genotype <- read.table("genotype2_02_metag_artro.txt", header=T, sep="\t", check.names=F,row.names=1)
p <- ggtree(dend) %<+% info + xlim(0, 0.5)+ theme_tree2()
p2 <- p + geom_tiplab(offset = .005,size=3, color="black") +
  geom_tippoint(aes(color=size),size=4,alpha=0.5)
gheatmap(p2, genotype, offset=0.1, width=0.1, low="white", high="black", colnames_position = "top",  color="black",font.size=3)

#----shannon vs seq per volume 16S----
data16S_day1 <- read.table("16Sdada_div_day1.txt", header=T, sep="\t", check.names=F,row.names=1)
#Define ggplot2 custom colors:
allGroupsColors<- c( "#F8766D", "#00BFC4", "#C77CFF","#7CAE00")
# New facet label names for size variable
size.labs <- c("0.22", "0.22-3", "20", "3-20")
names(size.labs) <- c("A02", "B023", "D20", "C320")
#define mean line
mean_H <- data.frame(Size = c("A02", "B023", "D20", "C320"), H = c(5.20, 5.18, 5.40,5.29))

ggplot(data16S_day1, aes(x=seqpervol, y=Shannon_H , color=Size)) + geom_point(shape=19, size=5,alpha=1/2)+ theme_bw()+
  scale_colour_hue(l=50) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Sequence per L (log scale)",y='Shannon index')+
  theme(legend.position = "none")+
  scale_color_manual(values = allGroupsColors)+
  scale_x_continuous(trans='log2',breaks=c(5000, 50000, 500000)) +
  geom_smooth(method=lm,se=FALSE, size=0.3,data = subset(data16S_day1, Size =="C320"))+# Don't add shaded confidence region
  facet_wrap( ~Size, ncol=2,labeller = labeller(Size = size.labs))+# Divide by Size, going horizontally and wrapping with 2 columns
  geom_hline(aes(yintercept = H), mean_H,size=0.2)




#----taxa per seq vs volume 16S----
#Define ggplot2 custom colors:
allGroupsColors<- c( "#F8766D", "#00BFC4", "#C77CFF","#7CAE00")
# New facet label names for size variable
size.labs <- c("0.22", "0.22-3", "20", "3-20")
names(size.labs) <- c("A02", "B023", "D20", "C320")
#define mean line
mean_H <- data.frame(Size = c("A02", "B023", "D20", "C320"), H = c(8.29, 8.52, 8.49,8.48))

ggplot(data16S_day1, aes(x=volume, y=taxaperseq , color=Size)) + geom_point(shape=19, size=5,alpha=1/2)+ theme_bw()+
  scale_colour_hue(l=50) +
  expand_limits(x=c(0,1000))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Volume (L)",y='Number of ASVs (normalized per 10 000 sequences)')+
  theme(legend.position = "none")+
  scale_color_manual(values = allGroupsColors)+
  scale_x_continuous(trans='log2',breaks=c(1,10,100, 1000)) + # Add log scale
  facet_wrap( ~ Size, ncol=2,labeller = labeller(Size = size.labs))+
  geom_hline(aes(yintercept = H), mean_H,size=0.2)

#----shannon vs seq per volume 18S----
data18S_day1<- read.table("18Sdada_div_day1.txt", header=T, sep="\t", check.names=F,row.names=1)

#Define ggplot2 custom colors:
allGroupsColors<- c( "#F8766D", "#00BFC4", "#C77CFF","#7CAE00")
# New facet label names for size variable
size.labs <- c("0.22", "0.22-3", "20", "3-20")
names(size.labs) <- c("A0.22", "B0.22-3", "D20", "C3-20")
#define mean line
mean_H <- data.frame(Size = c("A0.22", "B0.22-3", "D20", "C3-20"), H = c(4.86, 5.58, 2.80,5.69))


ggplot(data18S_day1, aes(x=seqpervol, y=Shannon_H , color=Size)) + geom_point(shape=19, size=5,alpha=1/2)+ theme_bw()+
  scale_colour_hue(l=50) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Sequence per L (log scale)",y='Shannon index')+
  theme(legend.position = "none")+
  scale_color_manual(values = allGroupsColors)+
  scale_x_continuous(trans='log2',breaks=c(5000, 50000, 500000)) +
  facet_wrap( ~Size, ncol=2,labeller = labeller(Size = size.labs))+# Divide by Size, going horizontally and wrapping with 2 columns
  geom_hline(aes(yintercept = H), mean_H,size=0.2)

#----taxa per seq vs volume 18S----
#Define ggplot2 custom colors:
allGroupsColors<- c( "#F8766D", "#00BFC4", "#C77CFF","#7CAE00")
# New facet label names for size variable
size.labs <- c("0.22", "0.22-3", "20", "3-20")
names(size.labs) <- c("A0.22", "B0.22-3", "D20", "C3-20")
#define mean line
mean_H <- data.frame(Size = c("A0.22", "B0.22-3", "D20", "C3-20"), H = c(10.68, 9.67, 6.43,14.95))

ggplot(data18S_day1, aes(x=vol, y=taxperseq , color=Size)) + geom_point(shape=19, size=5,alpha=1/2)+ theme_bw()+
  scale_colour_hue(l=50) +
  expand_limits(x=c(0,1000))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Volume (L)",y='Number of ASVs (normalized per 10 000 sequences)')+
  theme(legend.position = "none")+
  scale_color_manual(values = allGroupsColors)+
  scale_x_continuous(trans='log2',breaks=c( 10, 100,1000)) +
  facet_wrap( ~Size, ncol=2,labeller = labeller(Size = size.labs))+# Divide by Size, going horizontally and wrapping with 2 columns
geom_hline(aes(yintercept = H), mean_H,size=0.2)

###----taxa per seq vs volume mitag16S
datamitag <- read.table("mitag_div_day1_minus1.txt", header=T, sep="\t", check.names=F,row.names=1)
#Define ggplot2 custom colors:
allGroupsColors<- c( "#F8766D", "#00BFC4", "#C77CFF","#7CAE00")
# New facet label names for size variable
size.labs <- c("0.22", "0.22-3", "20", "3-20")
names(size.labs) <- c("A02", "B023", "D20", "C320")

ggplot(datamitag, aes(x=volume, y=taxperseq , color=Size)) + geom_point(shape=19, size=5,alpha=1/2)+ theme_bw()+
  scale_colour_hue(l=50) +
  labs(x = "Volume (L)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(legend.position = "none")+
  scale_color_manual(values = allGroupsColors)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Volume (L)",y='Number of mitags (normalized per 100 000 sequences)')+
  scale_x_continuous(trans='log2')+# Add log scale
  geom_smooth(method=lm,size=0.3, data = subset(datamitag,  Size =="A02"),   # Add linear regression lines
              se=FALSE)+
  geom_smooth(method=lm,size=0.3, data = subset(datamitag,  Size =="B023"),   # Add linear regression lines
              se=FALSE)+
  facet_wrap( ~ Size, ncol=2,labeller = labeller(Size = size.labs))# Divide by Size, going horizontally and wrapping with 2 columns

###----Shannon vs seq per volume mitag16S
datamitag <- read.table("mitag_div_day1_minus1.txt", header=T, sep="\t", check.names=F,row.names=1)
#Define ggplot2 custom colors:
allGroupsColors<- c( "#F8766D", "#00BFC4", "#C77CFF","#7CAE00")
# New facet label names for size variable
size.labs <- c("0.22", "0.22-3", "20", "3-20")
names(size.labs) <- c("A02", "B023", "D20", "C320")
#calculate mean
meanmitag16s<-ddply(datamitag, "Size", summarise, mean=mean(Shannon_H), sd=sd(Shannon_H), se = sd(Shannon_H) / sqrt((length(Size))))
meanmitag16s
#define mean line
mean_H <- data.frame(Size = c("A02", "B023", "D20", "C320"), H = c(7.107692, 7.086364, 6.652500,6.858455))

ggplot(datamitag, aes(x=seqpervol, y=Shannon_H , color=Size)) + geom_point(shape=19, size=5,alpha=1/2)+ theme_bw()+
  scale_colour_hue(l=50) +
  labs(x = "Volume (L)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(legend.position = "none")+
  scale_color_manual(values = allGroupsColors)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Sequence per L (log scale)",y='Shannon index')+
  scale_x_continuous(trans='log2',breaks=c(500000,5000000, 50000000))+# Add log 
  facet_wrap( ~ Size, ncol=2,labeller = labeller(Size = size.labs))+# Divide by Size, going horizontally and wrapping with 2 columns
geom_hline(aes(yintercept = H), mean_H,size=0.2)

#-------SIMKA----MDS
allGroupsColors<- c( "#F8766D" ,"#00BFC4", "#C77CFF","#7CAE00" )
dist_matrix<-read.csv(file="mat_abundance_braycurtis.csv", header=TRUE, row.names=1,sep=";")
info <- read.csv(paste0("tip_data_all_mitag.csv"), sep=";")
cr3 = as.matrix(dist_matrix)
#
mds<-metaMDS(cr3)
plot(mds)
MDS_xy <- data.frame(mds$points)
MDS_xy$size <- info$size
MDS_xy$size <- factor(MDS_xy$size, levels = c("0.22", "0.22-3", "3-20","20"))
ggplot(MDS_xy, aes(MDS1, MDS2, color = size)) + geom_point(size=7,alpha=0.5) + theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  expand_limits(x=c(0,0.75))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12), legend.title=element_text(size=14))+
  scale_color_manual(values = allGroupsColors)+
  labs(color="Size (\u03BCm)") +
  theme(legend.position = c(0.926, 0.864))+
  theme( legend.background = element_rect(color = "black", linetype = "solid",size=0.2))

#-------SIMKA----cluster
dist_matrix<-read.csv(file="mat_abundance_braycurtis.csv", header=TRUE, row.names=1,sep=";")
inv_cr3 = matrix(100, ncol=dim(cr3)[1], nrow=dim(cr3)[1]) - cr3
simka_distance = as.dist(cr3)
dendo_cr3 = hclust(simka_distance, method="average")
plot(dendo_cr3, main="Simka normalized analysis", sub = NA)

dend <- as.phylo(dendo_cr3)
info <- read.csv(paste0("tip_data_all_mitag.csv"), sep=";")
info$size <- factor(info$size, levels = c("0.22", "0.22-3", "3-20","20"))
genotype <- read.table("genotype2_all_metag_vol.txt", header=T, sep="\t", check.names=F,row.names=1)
p <- ggtree(dend)%<+% info +xlim(0,1)
p2 <- p + geom_tiplab(offset = .005,size=3, color="black") +
  geom_tippoint(aes(color=size),size=4,alpha=0.5) +
  scale_color_manual(values = allGroupsColors)+
  geom_treescale(x=0)
p2

gheatmap(p2, genotype, offset=0.15, width=0.1, low="white", high="black", colnames_position = "top", color="black", font.size=3)+
  scale_fill_manual(breaks=c("1-2.5","10-60", "100-120",  "500+"), 
                    values=c("slategray1","steelblue1" ,"steelblue4", "black"), name="Filter")

