# pemba ml trees
# pemba + all AF1b
library(treeio)
library(cowplot)
library(phytools)
library(strap)
library(ggtree)
library(ggtreeExtra)
library(wesanderson)
library(ggplot2)  #need this to change ggtree text
#remotes::install_github("JLSteenwyk/ggpubfigs")
#remotes::install_github("fmichonneau/phyloch")
library(ggpubfigs)
library(ggpubr)
library(ggnewscale)
library(viridis)
library(ggprism)
library(ggpubr)
library(ggtreeExtra)
library(Biostrings)
library(RColorBrewer)
library(ggnewscale)
library(ggmsa)

## DATA#################
### trees
gtree=read.tree("data/genetic/ALL.aln_g_pruned.fasta.treefile")
ntree=read.tree("data/genetic/ALL.aln_n_pruned.fasta.treefile")

## metadata
g.annot=read.csv("data/genetic//ALL.aln_g_meta.csv")
g.annot$sequence.gb_place_sampled[grep("RAB16|PM0|RAB12|RV2776|RV2777|RV2778|RV2817", g.annot$sequence.isolate, ignore.case=T)] = "Pemba Island"
g.annot$sequence.gb_place_sampled[grep("RV2782", g.annot$sequence.isolate, ignore.case=T)] = "Zanzibar"
g.annot$Pemba=NA
g.annot$Pemba[grep("RAB16|PM0",g.annot$sequence.isolate, ignore.case=T)]="Outbreak"
g.annot$Pemba[grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",g.annot$sequence.isolate, ignore.case=T)]="Historical"
n.annot=read.csv("data/genetic/ALL.aln_n_meta.csv")
n.annot$sequence.gb_place_sampled[grep("RAB16|PM0|RAB12|RV2776|RV2777|RV2778|RV2817", n.annot$sequence.isolate, ignore.case=T)] = "Pemba Island"
n.annot$sequence.gb_place_sampled[grep("RV2782", n.annot$sequence.isolate, ignore.case=T)] = "Zanzibar"
n.annot$Pemba=NA
n.annot$Pemba[grep("RAB16|PM0",n.annot$sequence.isolate, ignore.case=T)]="Outbreak"
n.annot$Pemba[grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",n.annot$sequence.isolate, ignore.case=T)]="Historical"


### alignments
n.aln=read.fasta("data/genetic/ALL.aln_n_pruned.fasta")
g.aln=read.fasta("data/genetic/ALL.aln_g_pruned.fasta")

## G gene ####################
### tree manipulations
#remove extreme outliers (long branches)
#a=ggtree(gtree)
#tail(sort(gtree$edge.length),50)
#remove really long branches from tree:
quantile=quantile(gtree$edge.length,probs=0.99)
drop=gtree$edge[which(gtree$edge.length>quantile),2]
#drop=gtree$edge[which(gtree$edge.length>0.04),2]
gtree.mod=drop.tip(gtree,tip=drop,2)
#and remove from metadata
g.annot.mod=g.annot[g.annot$sequence.sequenceID %in% gtree.mod$tip.label,]

# tree plots:  aesthetics, midpoint root, ladderize
gplot=midpoint.root(gtree.mod)
gplot <- ggtree(gtree.mod, ladderize = TRUE,col="gray", size=0.1) %<+% g.annot.mod

# highlight pemba seq- new and old
pemba.new=grep("RAB16|PM0",gplot$data$sequence.isolate,ignore.case=T)
pemba.old=grep("RAB12|RV2782|RV2776|RV2777|RV2778|RV2817",gplot$data$sequence.isolate, ignore.case=T)

# tree plots:  aesthetics, midpoint root, ladderize
gplot+ geom_tippoint(aes(color=sequence.m49_country.display_name), size=0.4)+ theme(legend.position = "left", legend.key.size = unit(0.2, 'cm'),legend.text = element_text(size=8),legend.title=element_blank())

## identify colour associations
# function that replicates gnplot col palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# based on length n
n = length(unique(g.annot.mod$sequence.m49_country.display_name))
# col palette 
cols = gg_color_hue(n)
#check col:
plot(1:n, pch = 16, cex = 2, col = cols)
#find which number col is Tza (+1 to account for NA):
tza=which(sort(unique(g.annot.mod$sequence.m49_country.display_name))=="Tanzania")+1

# change labels by adding new cols to attached dataframe
gplot$data$new=""
gplot$data$new[pemba.old]="+"
gplot$data$new[pemba.new]="*"
gplot$data$length=gplot$data$sequence.gb_length
# anything specified as length>1575 has come from wgs
gplot$data$length[gplot$data$length>1575]=1575

pm.mrca=getMRCA(gtree.mod, c(pemba.new,pemba.old))

# add bar tile to show country associations
# geom_tippoint(data=td_filter(node %in% pemba.old), color="blue", size=3,shape=8)
rect=gplot+
  geom_tippoint(data=td_filter(node %in% pemba.new), fill="darkred",col="black", size=2, shape=24, stroke=0.01, alpha=0.4)+
  geom_tippoint(data=td_filter(node %in% pemba.old), fill="black",col="black", size=2, shape=21, stroke=0.01, alpha=0.7)+
   # geom_tiplab(data=td_filter(node %in% c(pemba.old,pemba.new)),aes(label=sequence.isolate),size=0.4,align=F, linetype="blank")+
  geom_fruit(geom=geom_tile,mapping=aes(y=sequence.sequenceID,fill=sequence.m49_country.display_name), width=0.008,offset=0.07) +
  geom_treescale(fontsize = 6)+
  theme(legend.position = "left", legend.key.size = unit(0.2, 'cm'),legend.text = element_text(size=11),legend.title=element_blank())+
  guides(fill=guide_legend(ncol =4))+
  #layout_circular()+
  labs(title="A")+
  theme(plot.title = element_text(face="bold"))+
  geom_hilight(node=pm.mrca,, fill="gray", alpha=.4); rect
 # geom_hilight(node=pm.mrca, fill="gray", alpha=.6,type = "gradient", gradient.direction = 'rt'); rect
# geom_fruit(geom=geom_bar,mapping=aes(y=label, x=length, fill=sequence.m49_country.display_name), pwidth=0.38, orientation="y",  stat="identity",axis.params = list(axis = "x"))  ; rect

#rect+layout_fan(angle=270)
#rect+layout_circular()


###########
# N gene tree

#trees
#remove extreme outliers
#remove really long branches from tree:
quantile=quantile(ntree$edge.length,probs=0.99)
drop=ntree$edge[which(ntree$edge.length>quantile),2]
ntree.mod=drop.tip(ntree,tip=drop)
#and remove from metadata
drop_names=ntree$tip.label[drop]
n.annot.mod=n.annot[n.annot$sequence.sequenceID %in% ntree.mod$tip.label,]

#
# tree aesthetics, midpoint root, ladderize
nplot=midpoint.root(ntree.mod)
nplot <- ggtree(ntree.mod, ladderize = TRUE,col="gray", size=0.4) %<+% n.annot.mod
#nplot+ geom_tippoint(aes(color=sequence.m49_country.display_name), size=0.4)+ theme(legend.position = "left", legend.key.size = unit(0.2, 'cm'),legend.text = element_text(size=8),legend.title=element_blank())

#pemba seq- new and old
pemba.new.n=grep("RAB16|PM0",nplot$data$sequence.isolate, ignore.case=T)
pemba.old.n=grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817",nplot$data$sequence.isolate, ignore.case=T)
pemba.old.n2=grep("RAB12001|RV2782|RV2776",nplot$data$sequence.isolate, ignore.case=T)

## identify colour associationsL
# function that replicates gnplot col palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# based on length n
n = length(unique(n.annot.mod$sequence.m49_country.display_name))
# col palette 
cols = gg_color_hue(n)
#check col:
plot(1:n, pch = 16, cex = 2, col = cols)
#find which number col is Tza:
tza.n=which(sort(unique(n.annot.mod$sequence.m49_country.display_name))=="Tanzania") #same as tza

# change labels by adding new cols to attached dataframe
nplot$data$new=""
nplot$data$new[c(pemba.old.n, pemba.old.n2)]="+"
nplot$data$new[pemba.new.n]="*"
nplot$data$length=nplot$data$sequence.gb_length
# anything specified as length>1575 has come from wgs
nplot$data$length[nplot$data$length>1353]=1353

pm.mrca.n=getMRCA(ntree.mod, c(pemba.new.n,pemba.old.n))
pm.mrca.n1=getMRCA(ntree.mod, pemba.old.n)
pm.mrca.n2=getMRCA(ntree.mod, pemba.old.n2)
# add bar tile to show country associations
# geom_tippoint(data=td_filter(node %in% pemba.old), color="blue", size=3,shape=8)
rect2=nplot+
  geom_tippoint(data=td_filter(node %in% pemba.old.n2), fill="black",col="black", size=2, shape=21, stroke=0.01, alpha=0.4)+
  geom_tippoint(data=td_filter(node %in% pemba.new.n), fill="darkred",col="black", size=2, shape=24, stroke=0.01, alpha=0.7)+
  geom_tippoint(data=td_filter(node %in% pemba.old.n), fill="black",col="black", size=2, shape=21, stroke=0.01, alpha=0.7)+
  #geom_tiplab(data=td_filter(node %in% c(pemba.old.n,pemba.new.n, pemba.old.n2)),aes(label=sequence.isolate),size=0.4,align=F, linetype="blank")+
  geom_tiplab(data=td_filter(node %in% pemba.old.n),aes(label=sequence.isolate),size=0.4,colour=cols[13],align=F, linetype="blank")+
  geom_fruit(geom=geom_tile,mapping=aes(y=sequence.sequenceID,fill=sequence.m49_country.display_name), width=0.006,offset=0.05) +
  geom_treescale(fontsize = 6)+
  theme(legend.position = "none", legend.key.size = unit(0.2, 'cm'),legend.text = element_text(size=11),legend.title=element_blank())+
  guides(fill=guide_legend(ncol =4))+
  #layout_circular()+
  geom_hilight(node=pm.mrca.n, fill="gray", alpha=.6)+
  # geom_hilight(node=pm.mrca.n2, fill="gray", alpha=.6)+
  labs(title="B")+
  theme(plot.title = element_text(face="bold"))+
  geom_hilight(node=pm.mrca.n2, fill="gray", alpha=.4);rect2
#geom_hilight(node=pm.mrca.n2, fill="gray", alpha=.6); rect; rect2
# geom_fruit(geom=geom_bar,mapping=aes(y=label, x=length, fill=sequence.m49_country.display_name), pwidth=0.38, orientation="y",  stat="identity",axis.params = list(axis = "x"))  ; rect


ggarrange(rect,rect2, ncol=2, common.legend = TRUE, legend="bottom")              

### tree subsets
sub.g=getMRCA(gtree.mod, c(pemba.new, pemba.old))
sub.n=getMRCA(ntree.mod, c(pemba.new.n, pemba.old.n))
sub.n.2=getMRCA(ntree.mod, pemba.old.n2)

# g gene closest relatives excluding tz5 lineag
bi_subset <- tree_subset(gtree.mod, sub.g,levels_back = 0)
# n gene closest relatives excluding tz5 lineage
bi_subset_n <- tree_subset(ntree.mod, sub.n,levels_back = 1)
#loner tz5 lineage- n gene
bi_subset_n2 <- tree_subset(ntree.mod, sub.n.2,levels_back = 4)
#sub=g.annot.mod[g.annot.mod$sequence.sequenceID %in% bi_subset$tip.label,]

ggtree(bi_subset)  %<+% g.annot.mod+ 
  geom_tiplab(aes(label=sequence.isolate), size=0.4) +
  geom_tippoint(data=td_filter(node %in% grep("RAB16|PM0|RAB12",sequence.isolate,ignore.case=T)), fill=cols[14],col="gray", size=2, shape=21, stroke=0.01)

ggtree(bi_subset_n, size=0.4)  %<+% n.annot.mod+ 
  geom_tiplab(aes(label=sequence.isolate), size=0.4) +
  geom_tippoint(data=td_filter(node %in% pemba.old), fill=cols[14],col="gray", size=2, shape=21, stroke=0.01)+
  geom_tippoint(data=td_filter(node %in% pemba.new), fill=cols[14],col="gray", size=2, shape=21, stroke=0.01)


g.sub=g.aln[names(g.aln) %in% bi_subset$tip.label]

ggtree(bi_subset_n2)
n.sub=n.aln[names(n.aln) %in% c(bi_subset_n$tip.label,bi_subset_n2$tip.label)]
n.sub2=n.aln[names(n.aln) %in% bi_subset_n$tip.label]

##### bootstrapped sub trees
g_bs=read.iqtree("data/genetic/g_af1a_subset_v2.contree")
n_bs=read.iqtree("data/genetic/n_af1a_subset.contree")

## plot n tree
#drop the outgroup from plot
outgp=grep("Af1a",n_bs@phylo$tip.label,ignore.case=T)
n_bs2 <- root(n_bs, outgp)
n_bs3=treeio::drop.tip(n_bs2,tip=outgp)

n_bs_tree=ggtree(n_bs3,size=0.4, col="gray")  %<+% n.annot.mod

pemba.new.n.bs=grep("RAB16|PM0",n_bs_tree$data$sequence.isolate, ignore.case=T)
pemba.old.n.bs=grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817",n_bs_tree$data$sequence.isolate, ignore.case=T)
pemba.old.n2.bs=grep("RAB12001|RV2782|RV2776",n_bs_tree$data$sequence.isolate, ignore.case=T)

n_bs_tree$data$new=""
n_bs_tree$data$new[c(pemba.old.n.bs, pemba.old.n2.bs)]="+"
n_bs_tree$data$new[pemba.new.n.bs]="!"
n_bs_tree$data$length=n_bs_tree$data$sequence.gb_length
# anything specified as length>1575 has come from wgs
n_bs_tree$data$length[n_bs_tree$data$length>1353]=1353
n_bs_tree$data$Pemba=NA
n_bs_tree$data$Pemba[grep("RAB16|PM0",n_bs_tree$data$sequence.isolate, ignore.case=T)]="Outbreak"
n_bs_tree$data$Pemba[grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",n_bs_tree$data$sequence.isolate, ignore.case=T)]="Historical"
n_bs_tree$data$NGS=""
n_bs_tree$data$NGS[grep("RAB16|RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",n_bs_tree$data$sequence.isolate, ignore.case=T)]="Illumina_metagenomic"
n_bs_tree$data$NGS[grep("PM0|RAB16031|RAB16033",n_bs_tree$data$sequence.isolate, ignore.case=T)]="Minion_PCR"

dna= readDNAStringSet("data/genetic/ALL.aln_n_subset.fasta")
dna2=tidy_msa(dna)

nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)

n_sub_final=n_bs_tree+ geom_tippoint(data=subset(n_bs_tree$data, !is.na(Pemba)),aes(fill=Pemba,shape=Pemba), col="gray", size=2, stroke=0.01, alpha=0.7)+
  scale_fill_manual(values=c("black","darkred")) +
  scale_shape_manual(values=c(21,24)) +
  geom_nodepoint(color="darkgray", shape=18, alpha=0.6, size=1.5 , aes(subset=as.numeric(UFboot) >90))+
  new_scale_fill() +
  #geom_tippoint(data=td_filter(node %in% grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",sequence.isolate, ignore.case=T)), fill="black",col="gray", size=2, shape=24, stroke=0.01, alpha=0.6)+
   geom_tiplab(data=td_filter(node %in% grep("RAB16|PM0",sequence.isolate, ignore.case=T)),aes(label=sequence.isolate),size=1,align=F, linetype="blank")+
  #  geom_fruit(geom=geom_tile, mapping=aes(y=sequence.sequenceID,fill=NGS) ,width=0.001,offset=0.05)+
  #scale_fill_manual(values=c("white", "black","red")) +
  #geom_fruit(geom=geom_tile, mapping=aes(y=sequence.sequenceID,fill=Pemba) ,width=0.001,offset=0.05) +
  # scale_fill_manual(values=c("white", "black","red")) +
  geom_fruit(geom=geom_tile, mapping=aes(y=sequence.sequenceID,fill=sequence.m49_country.display_name), width=0.002,offset=0.05) +
  geom_treescale(fontsize = 6)+
  theme(legend.position = "bottom", legend.key.size = unit(0.2, 'cm'),legend.text = element_text(size=10),legend.title=element_blank());n_sub_final


#geom_hilight(node=pm.mrca.n2, fill="gray", alpha=.6); rect; rect2
# geom_fruit(geom=geom_bar,mapping=aes(y=label, x=length, fill=sequence.m49_country.display_name), pwidth=0.38, orientation="y",  stat="identity",axis.params = list(axis = "x"))  ; rect

msaplot(p=n_sub_final, fasta="data/genetic/ALL.aln_n_subset.fasta")

ggarrange(rect2+labs(title="")+layout_circular()+theme(plot.margin = margin(t = 0,r = 0,b = 0,l = 0)),n_sub_final, ncol=2, common.legend = TRUE, legend="bottom", labels="AUTO", heights=c(1,2), widths=c(1,2))  

#getMRCA(n_bs3@phylo,grep("PM015|PM018|PM019|PM024",n_bs_tree$data$sequence.isolate, ignore.case=T))
#p2=viewClade(n_bs_tree, node=660) + geom_tiplab(aes(label=sequence.isolate))


# closest relatives - N

n.annot.mod$newname=gsub("_.*","",n.annot.mod$sequence.sequenceID)
for (i in 1:length(n.annot.mod$sequence.sequenceID)){
  if(!is.na(n.annot.mod$sequence.isolate[i]) & n.annot.mod$sequence.isolate[i]==gsub("_.*","",n.annot.mod$sequence.sequenceID[i])) {n.annot.mod$newname[i]=paste0("-")}}
n.annot.mod$sequence.gb_place_sampled[grep("SD",n.annot.mod$sequence.isolate)]="Serengeti District"

ln1.mrca=getMRCA(n_bs3@phylo,grep("Rab1602|Rab16031|PM023",n_bs_tree$data$sequence.isolate, ignore.case=T))

ln1<- tree_subset(n_bs3, ln1.mrca,levels_back = 2)
a=ggtree(ln1, col="black", size=0.4)%<+% n.annot.mod
#ggtree(ln1)%<+% n.annot.mod+geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate,sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1)+geom_fruit(geom=geom_msa,data=dna2, panel = 'msa', font=NULL,border=NA,color = "LETTER")

ln2.mrca=getMRCA(n_bs3@phylo,grep("PM015|PM018|PM019|PM024",n_bs_tree$data$sequence.isolate, ignore.case=T))

ln2<- tree_subset(n_bs3, ln2.mrca,levels_back = 1)
b=ggtree(ln2, col="black", size=0.4)%<+% n.annot.mod

ln3.mrca=grep("PM025",n_bs_tree$data$sequence.isolate, ignore.case=T)

ln3<- tree_subset(n_bs3, ln3.mrca,levels_back = 1)
c=ggtree(ln3, col="black", size=0.4)%<+% n.annot.mod

ln4.mrca=grep("RAB16033",n_bs_tree$data$sequence.isolate, ignore.case=T)

ln4<- tree_subset(n_bs3, ln4.mrca,levels_back = 2)
d=ggtree(ln4, col="black", size=0.4)%<+% n.annot.mod

ln5.mrca=getMRCA(n_bs3@phylo,grep("PM025|Rab1602|Rab16031|PM023",n_bs_tree$data$sequence.isolate, ignore.case=T))
#ln5.mrca=385

ln5<- tree_subset(n_bs3, ln5.mrca,levels_back = 1)
#ln5.mrca=385
#ln5.2<- tree_subset(ln5, ln5.mrca,levels_back = 1)
e=ggtree(ln5, col="black", size=0.4)%<+% n.annot.mod
#e=as.polytomy(n_bs3@phylo, feature="node.label", fun=function(x) as.numeric(x) < 70)
#e=ggtree(e, col="black", size=0.4)%<+% n.annot.mod

a1=a+
  geom_rootedge(rootedge=0.001)+
  geom_tippoint(size=1)+#geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate, sequence.m49_country.id, sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1,align=T, linesize=0.1)+
  geom_tippoint(data=td_filter(node %in% grep("Rab1602|Rab16031|PM023",sequence.isolate, ignore.case=T)),col="red")+
  geom_treescale(fontsize=4)+
  geom_tiplab(aes(color=(grepl("Rab1602|Rab16031|PM023",sequence.isolate,ignore.case=T)),label=paste(sequence.isolate,newname, sequence.m49_country.display_name, sequence.gb_place_sampled,sequence.collection_year,sep=" | ")),size=2,align=T, linesize=0.1)+
  scale_color_manual(values=c("black","red"),labels=c("Closest relatives","Pemba Outbreak"))+ guides(color=guide_legend(title=NULL))+
  # geom_fruit(geom=geom_msa,data=dna2, panel = 'msa', font=NULL,border=NA,use_dot=T,color = "Chemistry_NT", offset=0.5, ignore_gaps=T)+
  geom_nodelab(aes(label=label, node= getMRCA(ln1@phylo,tip=ln1@phylo$tip.label)),color="dark gray", alpha=1,size=2.5,hjust=0)+
  #subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color="dark gray", alpha=1,size=2.5,hjust=1)
  theme(legend.position="top")+
  geom_rootpoint(pch=0)+
  xlim(-0.001,0.02); a1


b1=b+
  geom_rootedge(rootedge=0.001)+
  geom_tippoint(size=1)+#geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate, sequence.m49_country.id, sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1,align=T, linesize=0.1)+
  geom_tippoint(data=td_filter(node %in% grep("PM015|PM018|PM019|PM024",sequence.isolate, ignore.case=T)),col="red")+
  geom_treescale(fontsize=4,offset=0.1)+
  geom_tiplab(aes(color=(grepl("PM015|PM018|PM019|PM024",sequence.isolate,ignore.case=T)),label=paste(sequence.isolate,newname, sequence.m49_country.display_name, sequence.gb_place_sampled,sequence.collection_year,sep=" | ")),size=1.5,align=T, linesize=0.1)+
  scale_color_manual(values=c("black","red"),labels=c("Closest relatives","Pemba Outbreak-cluster 1"))+ guides(color=guide_legend(title=NULL))+
 # geom_fruit(geom=geom_msa,data=dna2, panel = 'msa', font=NULL,border=NA,use_dot=T,color = "Chemistry_NT", offset=0.5, ignore_gaps=T)+
  geom_nodelab(aes(label=UFboot, subset=!is.na(as.numeric(UFboot))& node==c(26:30)),color="dark gray")+
 # geom_nodelab(aes(label=label, node= getMRCA(ln2@phylo,tip=ln2@phylo$tip.label)),color="dark gray", alpha=1,size=2.5,hjust=0)+
  #subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color="dark gray", alpha=1,size=2.5,hjust=1)
  theme(legend.position="top")+
 # geom_rootpoint(pch=0)+
  xlim(-0.001,0.02); b1

c1=c+
  geom_rootedge(rootedge=0.001)+
  geom_tippoint(size=1)+#geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate, sequence.m49_country.id, sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1,align=T, linesize=0.1)+
  geom_tippoint(data=td_filter(node %in% grep("PM025",sequence.isolate, ignore.case=T)),col="red")+
  geom_treescale(fontsize=4)+
  geom_tiplab(aes(color=(grepl("PM025",sequence.isolate,ignore.case=T)),label=paste(sequence.isolate,newname, sequence.m49_country.display_name, sequence.gb_place_sampled,sequence.collection_year,sep=" | ")),size=2,align=T, linesize=0.1)+
  scale_color_manual(values=c("black","red"),labels=c("Closest relatives","Pemba Outbreak"))+ guides(color=guide_legend(title=NULL))+
  # geom_fruit(geom=geom_msa,data=dna2, panel = 'msa', font=NULL,border=NA,use_dot=T,color = "Chemistry_NT", offset=0.5, ignore_gaps=T)+
  geom_nodelab(aes(label=label, node= getMRCA(ln3@phylo,tip=ln3@phylo$tip.label)),color="dark gray", alpha=1,size=2.5,hjust=0)+
  #subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color="dark gray", alpha=1,size=2.5,hjust=1)
  theme(legend.position="top")+
  geom_rootpoint(pch=0)+
  xlim(-0.001,0.02); c1

d1=d+
  geom_rootedge(rootedge=0.001)+
  geom_tippoint(size=1)+#geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate, sequence.m49_country.id, sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1,align=T, linesize=0.1)+
  geom_tippoint(data=td_filter(node %in% grep("RAB16033",sequence.isolate, ignore.case=T)),col="red")+
  geom_treescale(fontsize=4)+
  geom_tiplab(aes(color=(grepl("RAB16033",sequence.isolate,ignore.case=T)),label=paste(sequence.isolate,newname, sequence.m49_country.display_name, sequence.gb_place_sampled,sequence.collection_year,sep=" | ")),size=2,align=T, linesize=0.1)+
  scale_color_manual(values=c("black","red"),labels=c("Closest relatives","Pemba Outbreak"))+ guides(color=guide_legend(title=NULL))+
  # geom_fruit(geom=geom_msa,data=dna2, panel = 'msa', font=NULL,border=NA,use_dot=T,color = "Chemistry_NT", offset=0.5, ignore_gaps=T)+
  geom_nodelab(aes(label=label, node= getMRCA(ln4@phylo,tip=ln4@phylo$tip.label)),color="dark gray", alpha=1,size=2.5,hjust=0)+
  #subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color="dark gray", alpha=1,size=2.5,hjust=1)
  theme(legend.position="top")+
  geom_rootpoint(pch=0)+
  xlim(-0.001,0.04); d1

e1=e+
  geom_rootedge(rootedge=0.001)+
  geom_tippoint(size=1)+#geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate, sequence.m49_country.id, sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1,align=T, linesize=0.1)+
  geom_tippoint(data=td_filter(node %in% grep("PM025|Rab1602|Rab16031|PM023",sequence.isolate, ignore.case=T)),col="blue")+
  geom_treescale(fontsize=4,offset=0.1)+
  geom_tiplab(aes(color=(grepl("PM025|Rab1602|Rab16031|PM023",sequence.isolate,ignore.case=T)),label=paste(sequence.isolate,newname, sequence.m49_country.display_name, sequence.gb_place_sampled,sequence.collection_year,sep=" | ")),size=1.5,align=T, linesize=0.1)+
  scale_color_manual(values=c("black","blue"),labels=c("Closest relatives","Pemba Outbreak-cluster 2"))+ guides(color=guide_legend(title=NULL))+
  # geom_fruit(geom=geom_msa,data=dna2, panel = 'msa', font=NULL,border=NA,use_dot=T,color = "Chemistry_NT", offset=0.5, ignore_gaps=T)+
  #geom_nodelab(aes(label=UFboot))
geom_nodelab(aes(label=UFboot, subset=!is.na(as.numeric(UFboot))& node==c(49:54,87,88)),color="dark gray")+
 # geom_nodelab(aes(label=UFboot, node=getMRCA(ln5@phylo,tip=ln5@phylo$tip.label)),color="dark gray", alpha=1,size=2.5,hjust=0)+
  #subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color="dark gray", alpha=1,size=2.5,hjust=1)
  theme(legend.position="top")+geom_rootedge(rootedge=0.0001)+
  #geom_rootpoint(pch=0)+
  xlim(-0.001,0.04); e1 


ggarrange(b1,e1, common.legend = T, labels = "AUTO", widths=c(0.9,1))  


# closest relatives- G

g.annot.mod$newname=gsub("_.*","",g.annot.mod$sequence.sequenceID)
for (i in 1:length(g.annot.mod$sequence.sequenceID)){
  if(!is.na(g.annot.mod$sequence.isolate[i]) & g.annot.mod$sequence.isolate[i]==gsub("_.*","",g.annot.mod$sequence.sequenceID[i])) {g.annot.mod$newname[i]=paste0("-")}}
g.annot.mod$sequence.gb_place_sampled[grep("SD",g.annot.mod$sequence.isolate)]="Serengeti District"

outgp=grep("Af1a",g_bs@phylo$tip.label,ignore.case=T)
g_bs2 <- root(g_bs, outgp)
g_bs3=treeio::drop.tip(g_bs2,tip=outgp)

g_bs_tree=ggtree(g_bs3,size=0.4, col="gray")  %<+% g.annot.mod

pemba.new.g.bs=grep("RAB16|PM0",g_bs_tree$data$sequence.isolate, ignore.case=T)
pemba.old.g.bs=grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817",g_bs_tree$data$sequence.isolate, ignore.case=T)
pemba.old.g2.bs=grep("RAB12001|RV2782|RV2776",g_bs_tree$data$sequence.isolate, ignore.case=T)

g_bs_tree$data$new=""
g_bs_tree$data$new[c(pemba.old.g.bs, pemba.old.g2.bs)]="+"
g_bs_tree$data$new[pemba.new.g.bs]="!"
g_bs_tree$data$length=g_bs_tree$data$sequence.gb_length
# anything specified as length>1575 has come from wgs
g_bs_tree$data$length[g_bs_tree$data$length>1353]=1353
g_bs_tree$data$Pemba=NA
g_bs_tree$data$Pemba[grep("RAB16|PM0",g_bs_tree$data$sequence.isolate, ignore.case=T)]="Outbreak"
g_bs_tree$data$Pemba[grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",g_bs_tree$data$sequence.isolate, ignore.case=T)]="Historical"
g_bs_tree$data$NGS=""
g_bs_tree$data$NGS[grep("RAB16|RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",g_bs_tree$data$sequence.isolate, ignore.case=T)]="Illumina_metagenomic"
g_bs_tree$data$NGS[grep("PM0|RAB16031|RAB16033",g_bs_tree$data$sequence.isolate, ignore.case=T)]="Minion_PCR"

dna.g= readDNAStringSet("data/genetic/ALL.aln_g_subset.fasta")
dna2.g=tidy_msa(dna)

nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)

g_sub_final=g_bs_tree+ geom_tippoint(data=subset(g_bs_tree$data, !is.na(Pemba)),aes(fill=Pemba,shape=Pemba), col="gray", size=2, stroke=0.01, alpha=0.7)+
  scale_fill_manual(values=c("black","darkred")) +
  scale_shape_manual(values=c(21,24)) +
  geom_nodepoint(color="darkgray", shape=18, alpha=0.6, size=1.5 , aes(subset=as.numeric(UFboot) >90))+
  new_scale_fill() +
  #geom_tippoint(data=td_filter(node %in% grep("RAB12005|RAB12039|Rab12002|RV2777|RV2778|RV2817|RAB12001|RV2782|RV2776",sequence.isolate, ignore.case=T)), fill="black",col="gray", size=2, shape=24, stroke=0.01, alpha=0.6)+
   geom_tiplab(data=td_filter(node %in% grep("RAB16|PM0",sequence.isolate, ignore.case=T)),aes(label=sequence.isolate),size=1,align=F, linetype="blank")+
  #  geom_fruit(geom=geom_tile, mapping=aes(y=sequence.sequenceID,fill=NGS) ,width=0.001,offset=0.05)+
  #scale_fill_manual(values=c("white", "black","red")) +
  #geom_fruit(geom=geom_tile, mapping=aes(y=sequence.sequenceID,fill=Pemba) ,width=0.001,offset=0.05) +
  # scale_fill_manual(values=c("white", "black","red")) +
  geom_fruit(geom=geom_tile, mapping=aes(y=sequence.sequenceID,fill=sequence.m49_country.display_name), width=0.002,offset=0.05) +
  geom_treescale(fontsize = 6)+
  # new_scale_fill() +
  #geom_fruit(geom=geom_msa,data=dna2, panel = 'msa', font=NULL,border=NA,color = "LETTER")+
  theme(legend.position = "bottom", legend.key.size = unit(0.2, 'cm'),legend.text = element_text(size=10),legend.title=element_blank());g_sub_final

lg1.mrca=getMRCA(g_bs3@phylo,grep("Rab1602|Rab16031|PM023|PM025",g_bs_tree$data$sequence.isolate, ignore.case=T))

lg1<- tree_subset(g_bs3, lg1.mrca,levels_back = 3)
a.g=ggtree(lg1, col="black", size=0.4)%<+% g.annot.mod
#ggtree(lg1)%<+% g.annot.mod+geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate,sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1)+
  geom_fruit(geom=geom_msa,data=dna2.g, panel = 'msa', font=NULL,border=NA,color = "LETTER")

lg2.mrca=getMRCA(g_bs3@phylo,grep("PM015|PM018|PM019|PM024|RAB16033",g_bs_tree$data$sequence.isolate, ignore.case=T))

lg2<- tree_subset(g_bs3, lg2.mrca,levels_back = 0)
b.g=ggtree(lg2, col="black", size=0.4)  %<+% g.annot.mod 


a.g$data$sequence.gb_place_sampled[grep("SD",a.g$data$sequence.isolate)]="Serengeti District"
b.g$data$sequence.gb_place_sampled[grep("SD",b.g$data$sequence.isolate)]="Serengeti District"


a1.g=a.g+
  geom_rootedge(rootedge=0.001)+
  geom_tippoint(size=1)+#geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate, sequence.m49_country.id, sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1,align=T, linesize=0.1)+
  geom_tippoint(data=td_filter(node %in% grep("Rab1602|Rab16031|PM023|PM025",sequence.isolate, ignore.case=T)),col="blue")+
  geom_treescale(fontsize=4,offset=0.1)+
  geom_tiplab(aes(color=(grepl("Rab1602|Rab16031|PM023|PM025",sequence.isolate,ignore.case=T)),label=paste(sequence.isolate,newname, sequence.m49_country.display_name, sequence.gb_place_sampled,sequence.collection_year,sep=" | ")),size=1.5,align=T, linesize=0.1)+
  scale_color_manual(values=c("black","blue"),labels=c("Closest relatives","Pemba Outbreak-cluster 2"))+ guides(color=guide_legend(title=NULL))+
  # geom_fruit(geom=geom_msa,data=dna2.g, panel = 'msa', font=NULL,border=NA,use_dot=T,color = "Chemistry_NT", offset=0.5, ignore_gaps=T)+
  #geom_nodelab(aes(label=UFboot))
  geom_nodelab(aes(label=UFboot, subset=node==c(28:30)),color="dark gray")+
  geom_nodelab(aes(label=UFboot, subset=node==20),color="dark gray")+
 # geom_nodelab(aes(label=label, node= getMRCA(lg2@phylo,tip=lg2@phylo$tip.label)),color="dark gray", alpha=1,size=2.5,hjust=0)+
  #subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color="dark gray", alpha=1,size=2.5,hjust=1)
  theme(legend.position="top")+geom_rootedge(rootedge=0.0001)+
 # geom_rootpoint(pch=0)+
  xlim(-0.001,0.02); a1.g

b1.g=b.g+
  geom_rootedge(rootedge=0.001)+
  geom_tippoint(size=1)+#geom_tiplab(aes(label=paste(sequence.isolate,sequence.isolate, sequence.m49_country.id, sequence.gb_place_sampled,sequence.collection_year,sep="/ ")), size=1,align=T, linesize=0.1)+
  geom_tippoint(data=td_filter(node %in% grep("PM015|PM018|PM019|PM024|RAB16033",sequence.isolate, ignore.case=T)),col="red")+
  geom_treescale(fontsize=4,offset=0.1)+
  geom_tiplab(aes(color=(grepl("PM015|PM018|PM019|PM024|RAB16033",sequence.isolate,ignore.case=T)),label=paste(sequence.isolate,newname, sequence.m49_country.display_name, sequence.gb_place_sampled,sequence.collection_year,sep=" | ")),size=1.5,align=T, linesize=0.1)+
  scale_color_manual(values=c("black","red"),labels=c("Closest relatives","Pemba Outbreak-cluster 1"))+ guides(color=guide_legend(title=NULL))+
  # geom_fruit(geom=geom_msa,data=dna2.g, panel = 'msa', font=NULL,border=NA,use_dot=T,color = "Chemistry_NT", offset=0.5, ignore_gaps=T)+
  #geom_nodelab(aes(label=UFboot))
  geom_nodelab(aes(label=UFboot,subset=node<45|node==62|node==75),color="dark gray")+
 # geom_nodelab(aes(label=UFboot, subset=node==20))+
 # geom_nodelab(aes(label=UFboot, node= getMRCA(lg2@phylo,tip=lg2@phylo$tip.label)),color="dark gray", alpha=1,size=2.5,hjust=0)+
  #geom_nodepoint(aes(label=label, node= getMRCA(lg2@phylo,tip=lg2@phylo$tip.label)),color="dark gray", pch=15)+
  #subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),color="dark gray", alpha=1,size=2.5,hjust=1)
  theme(legend.position="top")+
#  geom_rootpoint(pch=0)+
  xlim(-0.001,0.02);b1.g

ggarrange(a1.g,b1.g, common.legend = T, labels = "AUTO")

ggarrange(b1,e1,a1.g,b1.g, common.legend = T, labels = "AUTO")
#cluster 1
ggarrange(b1,b1.g, common.legend = T, labels = c("A) Nucleoprotein", "B) Glycoprotein"), vjust=1, font.label = list(size = 12, face = "bold"))
#cluster 2
ggarrange(e1,a1.g, common.legend = T, labels = c("A) Nucleoprotein", "B) Glycoprotein"), vjust=1, font.label = list(size = 12, face = "bold"))

#cluster 1
c1=ggarrange(b1,b1.g, common.legend = T, labels = c("A) Nucleoprotein", "B) Glycoprotein"), vjust=1, font.label = list(size = 12, face = "bold"))
#cluster 2
c2=ggarrange(e1,a1.g, common.legend = T, labels = c("C) Nucleoprotein", "D) Glycoprotein"), vjust=1, font.label = list(size = 12, face = "bold"), widths=c(0.9,1))
c3=ggarrange(b1,b1.g, labels="AUTO", font.label = list(size = 12, face = "bold"),common.legend = T)
c4=ggarrange(e1,a1.g, labels=c("C","D"), font.label = list(size = 12, face = "bold"),common.legend = T)

ggarrange(c1,c2,nrow=2)
ggarrange(c3,c4,nrow=2)


