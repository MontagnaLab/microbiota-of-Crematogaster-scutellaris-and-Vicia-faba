setwd("D:/Progetti attuali/FormicheOscar/_ToPaper/HillComputation")

##### load data #####
bact <- read.csv("16S-CleanTable-SRSnorm.tsv", header=TRUE, row.names=1, sep="\t", stringsAsFactors = F)
fung <- read.csv("ITS-CleanTable-SRSnorm.tsv", header=TRUE, row.names=1, sep="\t", stringsAsFactors = F)


##### hilldiv analyses #####
library("hilldiv")
library("ggplot2")

vec=list(bact,fung)
names(vec)=c("bacteria","fungi")

# remove singletons
for (i in 1:length(vec)) {
  vec[[i]][vec[[i]]<2]=0
  vec[[i]]=vec[[i]][rowSums(vec[[i]])!=0,]
}

## alpha Diversity
for (i in 1:length(vec)) {
  
  xx=data.frame(hill_div(vec[[i]],0),
                hill_div(vec[[i]],1),
                hill_div(vec[[i]],2))
  names(xx)=c(paste0("Hill-q0_",names(vec)[i]),
              paste0("Hill-q1_",names(vec)[i]),
              paste0("Hill-q2_",names(vec)[i]))
  write.table(xx,paste0("HillAlphaDiversity_",names(vec)[i],".tsv"),
              quote=F,sep="\t",col.names = NA)
}


## alpha Diversity (single tables)
for (i in 1:length(vec)) {
  
  xx=data.frame(hill_div(vec[[i]],0))
  names(xx)="Hill_q0"
  write.table(xx,paste0("HillAlphaDiversity_singleTab_",names(vec)[i],"_q0.tsv"),
              quote=F,sep="\t",col.names = NA)
  
  xx=data.frame(hill_div(vec[[i]],1))
  names(xx)="Hill_q1"
  write.table(xx,paste0("HillAlphaDiversity_singleTab_",names(vec)[i],"_q1.tsv"),
              quote=F,sep="\t",col.names = NA)
  
  xx=data.frame(hill_div(vec[[i]],2))
  names(xx)="Hill_q2"
  write.table(xx,paste0("HillAlphaDiversity_singleTab_",names(vec)[i],"_q2.tsv"),
              quote=F,sep="\t",col.names = NA)
  
}



### Beta diversity

library(Matrix)

for (i in 1:length(vec)) {
  print(paste0("starting with ",names(vec)[i]))
  
  xx=list(pair_dis(vec[[i]],0),
          pair_dis(vec[[i]],1),
          pair_dis(vec[[i]],2))
  
  names(xx)=c(paste0("Hill-q0-betaDis_",names(vec)[i]),
              paste0("Hill-q1-betaDis_",names(vec)[i]),
              paste0("Hill-q2-betaDis_",names(vec)[i]))
  
  saveRDS(xx,paste0("HillBetaDissmilarities-ALL_",names(vec)[i],".rds"))
  
  
  yy=Matrix::forceSymmetric(xx[[1]]$L1_CqN,uplo="L")
  yy[is.na(yy)]=0
  write.table(as.matrix(yy),paste0("HillBetaDissmilarities-q0-Cqn_",names(vec)[i],".tsv"),
              quote=F,sep="\t",col.names = NA)
  
  yy=Matrix::forceSymmetric(xx[[1]]$L1_UqN,uplo="L")
  yy[is.na(yy)]=0
  write.table(as.matrix(yy),paste0("HillBetaDissmilarities-q0-Uqn_",names(vec)[i],".tsv"),
              quote=F,sep="\t",col.names = NA)

  yy=Matrix::forceSymmetric(xx[[2]]$L1_CqN,uplo="L")
  yy[is.na(yy)]=0
  write.table(as.matrix(yy),paste0("HillBetaDissmilarities-q1-Cqn_",names(vec)[i],".tsv"),
              quote=F,sep="\t",col.names = NA)

  yy=Matrix::forceSymmetric(xx[[2]]$L1_UqN,uplo="L")
  yy[is.na(yy)]=0
  write.table(as.matrix(yy),paste0("HillBetaDissmilarities-q1-Uqn_",names(vec)[i],".tsv"),
              quote=F,sep="\t",col.names = NA)
  
  yy=Matrix::forceSymmetric(xx[[3]]$L1_CqN,uplo="L")
  yy[is.na(yy)]=0
  write.table(as.matrix(yy),paste0("HillBetaDissmilarities-q2-Cqn_",names(vec)[i],".tsv"),
              quote=F,sep="\t",col.names = NA)
  
  yy=Matrix::forceSymmetric(xx[[3]]$L1_UqN,uplo="L")
  yy[is.na(yy)]=0
  write.table(as.matrix(yy),paste0("HillBetaDissmilarities-q2-Uqn_",names(vec)[i],".tsv"),
              quote=F,sep="\t",col.names = NA)
}





### Plots BetaDiversity
meta_bact <- read.csv("16S_hierarchy.tsv", header=T, sep="\t", stringsAsFactors = F)
meta_fung <- read.csv("ITS_hierarchy.tsv", header=T, sep="\t", stringsAsFactors = F)

beta_bact=list("q0"=pair_dis(vec$bacteria,0,hierarchy = meta_bact,level=2),
               "q1"=pair_dis(vec$bacteria,1,hierarchy = meta_bact,level=2),
               "q2"=pair_dis(vec$bacteria,2,hierarchy = meta_bact,level=2))

beta_fung=list("q0"=pair_dis(vec$fungi,0,hierarchy = meta_fung,level=2),
               "q1"=pair_dis(vec$fungi,1,hierarchy = meta_fung,level=2),
               "q2"=pair_dis(vec$fungi,2,hierarchy = meta_fung,level=2))

beta_bact_all=readRDS("HillBetaDissmilarities-ALL_bacteria.rds")
beta_fung_all=readRDS("HillBetaDissmilarities-ALL_fungi.rds")

col_bact=c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4","#FB9A99","#E31A1C","#CAB2D6","#6A3D9A")
col_fung=c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A")
colori=list(col_bact,col_fung)

#colori<-c("Foglie_C"="#B2DF8A","Foglie_T"="#33A02C",
#          "Formiche_a_C"="#A6CEE3","Formiche_a_T"="#1F78B4",
#          "Formiche_ct_C"="#FB9A99","Formiche_ct_T"="#E31A1C",
#          "Nettari_C"="#CAB2D6","Nettari_T"="#6A3D9A",
#          "Formiche_za_C"="#FDBF6F","Formiche_za_T"="#FF7F00")



dati=list("bacteria"=list("q0"=list("L1"=beta_bact_all$`Hill-q0-betaDis_bacteria`$L1_CqN,
                                    "L2"=beta_bact$q0$L2_CqN),
                          "q1"=list("L1"=beta_bact_all$`Hill-q1-betaDis_bacteria`$L1_CqN,
                                    "L2"=beta_bact$q1$L2_CqN),
                          "q2"=list("L1"=beta_bact_all$`Hill-q2-betaDis_bacteria`$L1_CqN,
                                    "L2"=beta_bact$q2$L2_CqN)),
          "fungi"=list("q0"=list("L1"=beta_fung_all$`Hill-q0-betaDis_fungi`$L1_CqN,
                                 "L2"=beta_fung$q0$L2_CqN),
                       "q1"=list("L1"=beta_fung_all$`Hill-q1-betaDis_fungi`$L1_CqN,
                                 "L2"=beta_fung$q1$L2_CqN),
                       "q2"=list("L1"=beta_fung_all$`Hill-q2-betaDis_fungi`$L1_CqN,
                                 "L2"=beta_fung$q2$L2_CqN)))
metadati=list("meta_bact"=meta_bact,
              "meta_fung"=meta_fung)


set.seed(1492)
for (tt in 1:length(dati)) {
  for (qq in 1:length(dati[[tt]])) {
    
    png(paste0("BetaDissimilarity_CqN_",names(dati)[tt],"_",names(dati[[tt]])[qq],"_qgplot.png"),
        width=30,height=20,units="cm",bg="white",res=800)
    qgplot=pair_dis_plot(dati[[tt]][[qq]]$L2,hierarchy = metadati[[tt]],
                         type="qgraph",level=2, magnify=T,colour = colori[[tt]])
    dev.off()
    
    png(paste0("BetaDissimilarity_CqN_",names(dati)[tt],"_",names(dati[[tt]])[qq],"_nmds.png"),
        width=21,height=14,units="cm",bg="white",res=800)
    nmdsplot=pair_dis_plot(dati[[tt]][[qq]]$L1,hierarchy = metadati[[tt]],
                           type="NMDS",level=1,colour = colori[[tt]])
    dev.off()


  }
}





################################ prove




library("RColorBrewer")

brewer.pal(10, "Paired")

# [1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C"
# [7] "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

# to display that palette:
display.brewer.pal(10, "Paired")

c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A")

col_bact=c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4","#FB9A99","#E31A1C","#CAB2D6","#6A3D9A")
col_fung=c("#B2DF8A","#33A02C","#A6CEE3","#1F78B4","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A")


xx<-c("bact","fung")

names(paste0(""))


## Qgraphs

meta_bact <- read.csv("16S_hierarchy.tsv", header=T, sep="\t", stringsAsFactors = F)
meta_fung <- read.csv("ITS_hierarchy.tsv", header=T, sep="\t", stringsAsFactors = F)

beta_bact=list("q0"=pair_dis(vec$bacteria,0,hierarchy = meta_bact,level=2),
               "q1"=pair_dis(vec$bacteria,1,hierarchy = meta_bact,level=2),
               "q2"=pair_dis(vec$bacteria,2,hierarchy = meta_bact,level=2))

beta_fung=list("q0"=pair_dis(vec$fungi,0,hierarchy = meta_fung,level=2),
               "q1"=pair_dis(vec$fungi,1,hierarchy = meta_fung,level=2),
               "q2"=pair_dis(vec$fungi,2,hierarchy = meta_fung,level=2))
  
pair_dis_plot(beta_bact$q2$L2_CqN,hierarchy = meta_bact,type="qgraph",level=2, magnify=T,colour = col_bact)
pair_dis_plot(beta_fung$q2$L2_CqN,hierarchy = meta_fung,type="qgraph",level=2, magnify=T,colour = col_fung)



## NMDS
library(vegan)

beta_bact_all=readRDS("HillBetaDissmilarities-ALL_bacteria.rds")
beta_fung_all=readRDS("HillBetaDissmilarities-ALL_fungi.rds")

pair_dis_plot(beta_bact_all$`Hill-q2-betaDis_bacteria`$L1_CqN,hierarchy = meta_bact,type="NMDS",level=1,colour = col_bact)
pair_dis_plot(beta_fung_all$`Hill-q0-betaDis_fungi`$L1_CqN,hierarchy = meta_fung,type="NMDS",level=1,colour = col_fung)



values.NMDS=metaMDS(as.dist(beta_bact_all$`Hill-q0-betaDis_bacteria`$L1_CqN),k=2,trymax=400)

NMDS = data.frame(x = values.NMDS$point[, 1], y = values.NMDS$point[,2], 
                  Sample = as.factor(rownames(values.NMDS$point)))
NMDS = merge(NMDS, meta_bact, by = "Sample")
colnames(NMDS) <- c("Sample", "x", "y", 
                    "Group")

nmds.plot <- ggplot() + geom_point(data = NMDS, aes_(x = ~x, 
                                                     y = ~y, colour = ~Group), size = 2, alpha = 0.5) + 
  scale_colour_manual(values = colour) + scale_shape_manual(values = 16) + 
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "grey"))
print(nmds.plot)
