# genetic distances
library(incidence)
library(lubridate)
library(stringr)
library(seqinr)
library(ape)
library(phangorn)
library(dplyr)
library(tidyr)
library(tibble)
library(plotly)
library(RColorBrewer)
library(ade4)
library(adegenet)

#---------
#alignment only- distance matrix based purely on genetic distance
all.fasta=read.dna("pemba_ea.aln.fasta", format = 'fasta')
#subset to pemba
pemba_sub=grep("Pemba",labels(all.fasta), value=T)
pemba_aln=all.fasta[labels(all.fasta) %in% pemba_sub,]
# genetic distance (="TN93" prob best)
alnDist <- dist.dna(all.fasta, model = "TN93", as.matrix = TRUE, pairwise.deletion = T)
pm_alnDist <- dist.dna(pemba_aln, model = "TN93", as.matrix = TRUE, pairwise.deletion = T)
#---------

#----------
# using max lik tree
mltree<- read.tree("pemba_ea_ml.treefile")

#GENETIC DISTANCE using tree information
#patristic distances i.e. takes into account the phylogeny by using branch lengths to represent evolutionary change and genetic difference.Values in the matrix are the sum of the branch lengths separating each pair of sequences .
#some patristic distances from ML trees are greater than the equivalent estimates of genetic distance based on pair-wise comparisons between the raw sequences
pat=cophenetic(mltree)
# to subset to pemba seq only:
sub=grep("Pemba",rownames(pat), value=F)
pat_pemba=pat[sub,sub]

#GENETIC CLUSTERS- ADEGENET
# whole tree
#note: 0.002 cutoff is somewhat arbitrary but output seemed to make sense! A better way might have been to use the 10th percentile of the pat distribution or something (i.e. quantile(pat, probs=0.10)), which is 0.0037
# you can try different cutoffs by changing in the command or run the gengraph code without the cutoff option to do it interactively (turn off show graph to be faster)
clust_all <- gengraph(pat, show.graph=F, col.pal=rainbow,truenames=T,cex=0.3, cutoff=0.002)
#write.csv(clust_all$clust$membership,"tanz_0.002clust_ml.csv")

## plot network graph (slow)
#plot(clust_all$g, edge.label=F, vertex.label.cex=0.3)
# 
# plot mltree coloured by cluster assignment
plot(mltree, tip.color=clust_all$col[clust_all$clust$membership], cex=0.2)
#add.scale.bar(x=0.001380052, y=-29.0968,length=0.01, cex=0.6)

# same for just pemba seq:
clust_pem <- gengraph(pat_pemba, show.graph=T, col.pal=rainbow,truenames=T,cex=0.3,cutoff=0.002) 
#write.csv(genclust_pembaOnly_0.002.csv")
plot(mltree, tip.color=clust_all$col[clust_pem$clust$membership], cex=0.3)

#--------------

# organise a bit for heatmap
patLong <- 
  pat %>% 
  as.data.frame(stringsToFactors = FALSE) %>% 
  rownames_to_column(var = "sample_1") %>% 
  gather(key = "sample_2", value = "distance", -sample_1, na.rm = TRUE) %>% 
  arrange(distance)

#alnDistLong %>% head()

#simple heatmap of pemba seq:
temp <- as.data.frame(as.matrix(pat_pemba))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5)

#fancier interactive heatmap
# organise data a bit first
# is massive with all data - use pat_pemba instead if want a smaller one
patLong <- 
  pat %>% 
  as.data.frame(stringsToFactors = FALSE) %>% 
  rownames_to_column(var = "sample_1") %>% 
  gather(key = "sample_2", value = "distance", -sample_1, na.rm = TRUE) %>% 
  arrange(distance)

# plot (note hold and drag on plot to zoom in)

patLong %>% 
  plot_ly(
    x = ~sample_2,
    y = ~sample_1,
    z = ~distance,
    type = "heatmap", colors = brewer.pal(11, "RdYlBu"), 
    zmin = 0.0, zmax = max(alnDistLong$distance),  xgap = 2, ygap = 1
  ) %>% 
  layout(
    margin = list(l = 100, r = 10, b = 100, t = 10, pad = 4), 
    yaxis = list(tickfont = list(size = 10), showspikes = TRUE),
    xaxis = list(tickfont = list(size = 10), showspikes = TRUE)
  )                     

