#micropsitta Fst and fixed differences investigation
library(vcfR)
library(StAMPP)
library(adegenet)
library(ggplot2)
library(gridExtra)

#read in vcf
vcf<-read.vcfR("/Users/devder/Desktop/micropsitta/micro.populations.snps.vcf")
vcf
vcfR::write.vcf(vcf.100, "~/Downloads/micropsitta.100.vcf.gz")


#convert vcfr to genlight
gen<- vcfR2genlight(vcf)

#read in popmap
pop<-read.table("/Users/devder/Desktop/micropsitta/pops.popmap.txt", sep="\t", header=FALSE)
colnames(pop)<-c("id","pop")
island.pop<-read.table("/Users/devder/Desktop/micropsitta/clean.popmap.txt", sep="\t", header=FALSE)
colnames(island.pop)<-c("id","pop")

#make pca
pca<-glPca(gen, nf=6)

#pull pca scores out of df
pca.scores<-as.data.frame(pca$scores)
rownames(pca.scores)
pop$sample
pca.scores$pop<-pop$pop

#ggplot color by species
ggplot(pca.scores, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(cex = 4, alpha=.75)+
  theme_classic()+
  theme(legend.position = "bottom", legend.key.size = unit(0, 'cm'),
        legend.title = element_text(size=8), legend.text = element_text(size=6.5))

#ggplot color by species
ggplot(pca.scores, aes(x=PC1, y=PC2, color=pop)) +
  geom_point(cex = 4, alpha=.75)+
  xlim(c(-12,-3))+
  ylim(c(1.5,6))+
  theme_classic()+
  theme(legend.position = "bottom", legend.key.size = unit(0, 'cm'),
        legend.title = element_text(size=8), legend.text = element_text(size=6.5))


#make splitstree
#remove invalid splitstree characters
gen@ind.names<-gsub(gen@ind.names, pattern = ".*_", replacement = "")
pop(gen)<-pop$pop
sample.div <- stamppNeisD(gen, pop = FALSE)
plot(nj(sample.div), type = "unrooted", cex = .65)
#export for splitstree
stamppPhylip(distance.mat=sample.div, file="~/Downloads/micro.splits.txt") 
#all samples splitstree
knitr::include_graphics(c("/Users/devder/Desktop/micropsitta/splitstree.png"))

#SNAPP tree with each of the main phylogenetic groupings designated as terminals
knitr::include_graphics(c("/Users/devder/Desktop/micropsitta/snapp.early.png"))

#dsuite test for introgression between these tips in a phylogenetic context
#heatmap of D statistic (ABBA/BABA statistic)
knitr::include_graphics(c("/Users/devder/Desktop/micropsitta/dsuite/7spec_BBAA_D.png"))
#heatmap of f4 statistic
knitr::include_graphics(c("/Users/devder/Desktop/micropsitta/dsuite/7spec_BBAA_f4ratio.png"))
#nothing approaching signficance

#calculate Fst and fixed differences
#calc pairwise Fst
gen@pop<-as.factor(pop$pop)
gen@ploidy<-as.integer(rep(2, times=24))
di.heat<-stamppFst(gen)
m<-di.heat$Fsts
#fill in upper triangle of matrix
m[upper.tri(m)] <- t(m)[upper.tri(m)]

#melt for plotting
heat <- reshape::melt(m)

#plot with labels
ggplot(data = heat, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()+
  geom_text(data=heat,aes(label=round(value, 2)))+
  theme_minimal()+
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1))

#identify the number of pairwise fixed differences between pops
mat<-extract.gt(vcf)
conv.mat<-mat
conv.mat[conv.mat == "0/0"]<-0
conv.mat[conv.mat == "0/1"]<-1
conv.mat[conv.mat == "1/1"]<-2
conv.mat<-as.data.frame(conv.mat)
#convert to numeric
for (i in 1:ncol(conv.mat)){
  conv.mat[,i]<-as.numeric(as.character(conv.mat[,i]))
}

#show colnames to verify you're subsetting correctly
colnames(conv.mat) == pop$id

#make vector to fill with fixed diff values
f<-c()

#write for loop to calc number of fixed diffs between each pop
for (i in 1:nrow(heat)){
  #calc af of pop1 and pop2
  pop1.af<-(rowSums(conv.mat[,pop$pop == heat$X1[i]], na.rm=T)/(rowSums(is.na(conv.mat[,pop$pop == heat$X1[i]]) == FALSE)))/2
  pop2.af<-(rowSums(conv.mat[,pop$pop == heat$X2[i]], na.rm=T)/(rowSums(is.na(conv.mat[,pop$pop == heat$X2[i]]) == FALSE)))/2
  #store number of fixed differences
  f[i]<-sum(is.na(abs(pop1.af - pop2.af)) == FALSE & abs(pop1.af - pop2.af) == 1) #find fixed SNPs and add to vector
}

#add number of fixed diffs to df
heat$fixed<-f

#fix the assignments
heat$mixed<-heat$value
heat$mixed[c(2:7,10:14,18:21,26:28,34:35,42)]<-heat$fixed[c(2:7,10:14,18:21,26:28,34:35,42)]

#plot with labels
fst.plot<-ggplot(data = heat, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()+
  #scale_x_discrete(limits=levels(heat$X2)[c(1,7,2,3,4,5,6)], labels=levels(heat$X1))+ 
  #scale_y_discrete(limits=levels(heat$X2)[c(1,7,2,3,4,5,6)], labels=levels(heat$X1))+
  geom_text(data=heat,aes(label=round(mixed, 2)), size=2.5)+
  theme_minimal()+
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="Fst") +
  theme(axis.text.x = element_text(angle = 45, vjust=.9, hjust = .9, size=12),
        axis.text.y = element_text(angle = 45, hjust = 1, size=12),
        axis.title.x = element_blank(), axis.title.y = element_blank())


#visualize heterozygosity and pi
pi<-read.table("~/Desktop/micropsitta/micro.populations.snps.p.sumstats_summary.tsv", header=T, sep='\t')
gen.mat<-as.matrix(gen)
loci<-rowSums(is.na(gen.mat) == FALSE)
hets<-rowSums(gen.mat == 1, na.rm = TRUE)/loci
het.df<-data.frame(id=row.names(gen.mat),
                   het=hets,
                   pop=pop$pop,
                   island=island.pop$pop)

#plot heterozygosity violin plots
het.plot<-ggplot(het.df, aes(x=pop, y=het)) + 
  #geom_violin(trim = FALSE)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 3, alpha=.75, aes(fill=pop, color=pop))+
  theme_classic()+
  #scale_fill_manual(values=c("blue","pink","red","purple","orange","green"))+
  #scale_color_manual(values=c("blue","pink","red","purple","orange","green"))+
  #scale_x_discrete(labels=c("California","Florida","Island","Sumichrast","Texas","Woodhouse"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "none")+
  geom_point(pi, mapping=aes(x=c(2,4,7,6,1,5,3), y=Pi), pch=8, cex=3)+
  labs(x="",y="heterozygosity")+
  scale_y_continuous(sec.axis = sec_axis(trans = (~.*1), name="Pi"))+
  coord_flip()

#calculate fixed, shared, and private SNPs
#logic for calculating shared, private, and fixed SNPs
#if target AF == 1, and rest AF==0, fixed non-ref
#if target AF == 0, and rest AF==1, fixed ref
#if target AF > 0 & < 1, and rest AF==0, private non-ref
#if target AF > 0 & < 1, and rest AF==1, private ref
#if rest AF > 0 & < 1, shared

#generate popmap file. 
#Two column popmap with the same format as stacks, and the columns must be named 'id' and 'pop'
###NOTE: popmap$pop MUST be a factor with levels
levels(pop$pop)
pop$pop<-as.factor(pop$pop)
levels(pop$pop)
#function for making pies
fix.private.share<-function(vcfR, popmap){
  
  #make vcfR into gt matrix
  gt.matrix<-extract.gt(vcfR)
  #convert values for easier computation
  gt.matrix[gt.matrix=="0/0"]<-0
  gt.matrix[gt.matrix=="0/1"]<-1
  gt.matrix[gt.matrix=="1/1"]<-2
  gt.matrix<-as.data.frame(gt.matrix)
  for (i in 1:ncol(gt.matrix)){
    gt.matrix[,i]<-as.numeric(as.character(gt.matrix[,i]))
  }
  
  #open df to fill up with information by pop
  df<-data.frame(pop=character(), class=character(), snps=numeric())
  #separate gt.matrix based on the population of interest
  for (i in 1:length(levels(popmap$pop))){
    target<-rowSums(gt.matrix[,as.character(popmap$id[popmap$pop==levels(popmap$pop)[i]])], na.rm = T)/
      (rowSums(!is.na(gt.matrix[,as.character(popmap$id[popmap$pop==levels(popmap$pop)[i]])]))*2)
    rest<-rowSums(gt.matrix[,as.character(popmap$id[popmap$pop!=levels(popmap$pop)[i]])], na.rm = T)/
      (rowSums(!is.na(gt.matrix[,as.character(popmap$id[popmap$pop!=levels(popmap$pop)[i]])]))*2)
    #remove sites where AF can't be calculated due to missing data
    missing.sites<-c(!is.na(target) & !is.na(rest))
    rest<-rest[missing.sites]
    target<-target[missing.sites]
    
    ##calculate # of fixed SNPs
    fix<-sum(abs(target-rest) == 1) #488
    
    #calc # of private non ref 
    y<-0
    for (j in 1:length(target)){
      if(target[j] > 0 & target[j] < 1 & rest[j] == 0){
        y<-y+1
      }
    }
    #calc # of private ref 
    for (k in 1:length(target)){
      if(target[k] > 0 & target[k] < 1 & rest[k] == 1){
        y<-y+1
      }
    }
    
    #calc # of shared SNPs
    z<-0
    for (m in 1:length(target)){
      if(rest[m] > 0 & rest[m] < 1){
        z<-z+1
      }
    }
    
    #combine this pop's info
    x<-cbind(rep(levels(popmap$pop)[i], times=3),
             c("fixed","private","shared"),
             c(fix,y,z))
    
    #add it to communal tidy df
    df<-as.data.frame(rbind(df,x))
  }
  #fix df colnames
  colnames(df)<-c("pop","class","snps")
  return(df)
}

#calc pie df
pie.df<-fix.private.share(vcf, pop)

#plot pies
#open empty list
plots<-list()
for (i in levels(as.factor(pie.df$pop))){
#plot
plots[[i]]<-ggplot(pie.df[pie.df$pop == i,], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
  geom_bar(stat="identity", width=1, color="black")+
  coord_polar("y", start=0)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("black","grey","white"))
}

#plot legend
legend <- cowplot::get_legend(ggplot(pie.df[pie.df$pop == "Makira",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
                                geom_bar(stat="identity", width=1, color="black")+
                                coord_polar("y", start=0)+
                                theme(legend.title = element_blank())+
                                scale_fill_manual(values=c("black","grey","white")))

#plot pies (order matches order of the 'plots' list)
pies<-grid.arrange(legend, plots[[7]],
             plots[[6]],
             plots[[5]],
             plots[[4]],
             plots[[3]],
             plots[[2]],
             plots[[1]], ncol=1)


#visualize this as a coherent figure
#make same figure and add compoplot
grid.arrange(fst.plot, het.plot, pies, ncol=3, widths=c(1,.7,.4))



