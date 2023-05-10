#extract partitions for beast
library(seqinr)
#BiocManager::install("ORFik")
#BiocManager::install("DECIPHER")
library(ORFik)
library(DECIPHER)
library(devtools)
#install_github("thibautjombart/apex")
#install.packages("apex")
library("apex")
library(stringr)
library("Biostrings")
library(sjPlot)
library(ggplot2)
library(sf)
library(reshape)
library("ggspatial")
library(viridis)
library(dplyr)

#input file
file="data/genetic/ALL.aln.fasta"

#output files
dir.create(file.path(paste(dirname(file), "output", sep="/")), showWarnings = FALSE)
# set output file location (same as input) and prefix
newfiles=paste(dirname(file), "output", gsub(".fasta|.fst","",basename(file)), sep="/")

# read in sequences as string set for ORF function
string=readDNAStringSet(file)
#remove any alignment errors (sometimes present if alignment has been subset). Write and replace file
string=RemoveGaps(string, removeGaps = "common")
names(string) <- gsub("New|","",labels(string),fixed=T)
names(string) <- vapply(strsplit(labels(string),".",fixed=T), `[`, 1, FUN.VALUE=character(1)) #extract sampleID from seq name
#writeXStringSet(string, paste0(newfiles, ".fasta"), format="fasta")

# find most complete genome and search for ORFs
chosen=which.max(str_count(string, "A|T|G|C"))
#chosen.seq=RemoveGaps(string[chosen], removeGaps = "common")
genes=as.data.frame(findORFs(string[chosen], startCodon = "ATG", minimumLength =200))

# it can't find the correct M gene start point so have to pull out manually
genes=genes[order(genes$start),]
genes=genes[-3,]
find.m=as.data.frame(findORFs(string[chosen], startCodon = "ATG", minimumLength =200, longestORF = F))
find.m=find.m[order(find.m$start),]
m=find.m[which(find.m$start>=2468 & find.m$width==609),]
#join with other genes
genes=rbind(genes, m)
genes=genes[order(genes$start),]
genes=genes[,-c(1,2)]
genes$gene=NA
genes$gene=c("n","p","m","g","l")

seq=read.fasta(paste0(newfiles, ".fasta"))

# split into coding partitions (5 genes)
# based on ORF positions
n=getFrag(seq, begin=genes$start[1],end=genes$end[1])
p=getFrag(seq, begin=genes$start[2],end=genes$end[2])
m=getFrag(seq, begin=genes$start[3],end=genes$end[3])
g=getFrag(seq, begin=genes$start[4],end=genes$end[4])
l=getFrag(seq, begin=genes$start[5],end=genes$end[5])

# output partitions as fasta files
#write.fasta(n,names=names(seq), paste(newfiles,"n.fasta",sep="_"))
#write.fasta(p,names=names(seq), paste(newfiles,"p.fasta",sep="_"))
#write.fasta(m,names=names(seq), paste(newfiles,"m.fasta",sep="_"))
#write.fasta(g,names=names(seq), paste(newfiles,"g.fasta",sep="_"))
#write.fasta(l,names=names(seq), paste(newfiles,"l.fasta",sep="_"))

summariseAlignment <- function(alignment){
  aligned <- read.alignment(alignment, format = "fasta")
  seq_data <- data.frame(ID = aligned$nam, N = NA, "gap" = NA, bases= NA, ambiguous=NA,
                         Length_before = nchar(aligned$seq[[1]]), Length_after = NA)
  
  for (i in 1:length(aligned$seq)) {
    seq_data$N[i] <- str_count(aligned$seq[[i]], pattern = 'n')
    seq_data$gap[i] <- str_count(aligned$seq[[i]], pattern = '-')
    seq_data$bases[i] <- sum(str_count(aligned$seq[[i]], pattern = c('a','c','t','g')))
    seq_data$ambiguous[i] <- (seq_data$Length_before[i] - seq_data$N[i] - seq_data$gap[i]-seq_data$bases[i])
    seq_data$Length_after[i] <- (seq_data$Length_before[i] - seq_data$N[i] - seq_data$gap[i])
  }
  return(seq_data)
}
pruneAlignment <- function(alignment){
  seq_data <- summariseAlignment(alignment)
  #eliminate bad seqeunces
  badSeq <- seq_data[which(seq_data$Length_after< 100),]
  fasta <- read.fasta(alignment)
  #remove bad seq:
  remove.badSeq=fasta[!names(fasta) %in% badSeq$ID]
  print(paste0("original:", length(fasta), "; removed: ",length(fasta)-length(remove.badSeq), "; remaining: ",length(remove.badSeq)))
  #remove.badEpi=badDate[which(badDate$published=="yes"),]
  return(remove.badSeq)
}

n.2=pruneAlignment(paste(newfiles,"n.fasta",sep="_"))
p.2=pruneAlignment(paste(newfiles,"p.fasta",sep="_"))
m.2=pruneAlignment(paste(newfiles,"m.fasta",sep="_"))
g.2=pruneAlignment(paste(newfiles,"g.fasta",sep="_"))
l.2=pruneAlignment(paste(newfiles,"l.fasta",sep="_"))
# 
# write.fasta(n.2,names=labels(n.2), paste(newfiles,"n_pruned.fasta",sep="_"))
# write.fasta(p.2,names=labels(p.2), paste(newfiles,"p_pruned.fasta",sep="_"))
# write.fasta(m.2,names=labels(m.2), paste(newfiles,"m_pruned.fasta",sep="_"))
# write.fasta(g.2,names=labels(g.2), paste(newfiles,"g_pruned.fasta",sep="_"))
# write.fasta(l.2,names=labels(l.2), paste(newfiles,"l_pruned.fasta",sep="_"))

###
meta=read.csv("data/genetic/AL_Cosmopolitan_AF1b_metadata_manual.csv", na.strings = c("NA","","-"))
set.n=meta[match(labels(n.2),meta$sequence.sequenceID,nomatch = 0),]
missing=labels(n.2)[which(!labels(n.2) %in% meta$sequence.sequenceID)] ;missing
#write.csv(set.n, paste(newfiles, "n_meta.csv", sep="_"), row.names=F)
#write.table(set.n, paste(newfiles, "n_meta.txt", sep="_"), row.names=F, sep="\t", quote = FALSE)

set.p=meta[match(labels(p.2),meta$sequence.sequenceID,nomatch = 0),]
missing=labels(p.2)[which(!labels(p.2) %in% meta$sequence.sequenceID)] ;missing
#write.csv(set.p, paste(newfiles, "p_meta.csv", sep="_"), row.names=F)
#write.table(set.p, paste(newfiles, "p_meta.txt", sep="_"), row.names=F, sep="\t", quote = FALSE)

set.l=meta[match(labels(l.2),meta$sequence.sequenceID,nomatch = 0),]
missing=labels(l.2)[which(!labels(l.2) %in% meta$sequence.sequenceID)] ;missing
#write.csv(set.l, paste(newfiles, "l_meta.csv", sep="_"), row.names=F)
#write.table(set.l, paste(newfiles, "l_meta.txt", sep="_"), row.names=F, sep="\t", quote = FALSE)

set.g=meta[match(labels(g.2),meta$sequence.sequenceID,nomatch = 0),]
missing=labels(g.2)[which(!labels(g.2) %in% meta$sequence.sequenceID)] ;missing
#write.csv(set.g, paste(newfiles, "g_meta.csv", sep="_"), row.names=F)
#write.table(set.g, paste(newfiles, "g_meta.txt", sep="_"), row.names=F, sep="\t", quote = FALSE)

set.m=meta[match(labels(m.2),meta$sequence.sequenceID,nomatch = 0),]
missing=labels(m.2)[which(!labels(m.2) %in% meta$sequence.sequenceID)] ;missing
#write.csv(set.m, paste(newfiles, "m_meta.csv", sep="_"), row.names=F)
#write.table(set.m, paste(newfiles, "m_meta.txt", sep="_"), row.names=F, sep="\t", quote = FALSE)

# summarise data
wgs=meta[meta$sequence.gb_length>=10730,]
a=c(length(n.2),length(p.2),length(m.2),length(g.2),length(l.2), nrow(wgs))
barplot(a)

set.n$length.cat<-cut(set.n$sequence.gb_length, c(0,1216,10730,11923),labels=c("pcr","gene", "wgs"))
set.p$length.cat<-cut(set.p$sequence.gb_length, c(0,804,10730,11923),labels=c("pcr","gene", "wgs"))
set.m$length.cat<-cut(set.m$sequence.gb_length, c(0,548,10730,11923),labels=c("pcr","gene", "wgs"))
set.g$length.cat<-cut(set.g$sequence.gb_length, c(0,1417,10730,11923),labels=c("pcr","gene", "wgs"))
set.l$length.cat<-cut(set.l$sequence.gb_length, c(0,5745,10730,11923),labels=c("pcr","gene", "wgs"))
wgs$length.cat<-cut(wgs$sequence.gb_length, c(0,10730,11923),labels=c("extra", "wgs"))
set.n$gene=NA; set.p$gene=NA;set.m$gene=NA;set.g$gene=NA;set.l$gene=NA; wgs$gene=NA
set.n$gene="n"; set.p$gene="p";set.m$gene="m";set.g$gene="g";set.l$gene="l";wgs$gene="wgs"
all.aln=rbind(set.n, set.p, set.l, set.g, set.m, wgs)


all.aln$gene <- factor(all.aln$gene,levels = c("wgs","n", "p", "m", "g","l"))
all.aln=all.aln[which(all.aln$sequence.m49_country.display_name!="Thailand"),]
ggplot(data=subset(all.aln, gene=="wgs" |length.cat=="gene"|length.cat=="pcr"), aes(x=gene,fill = length.cat)) + geom_bar()+
   ylab("Sequence count")+ labs(fill = "")+scale_fill_viridis(discrete=TRUE)

ggplot(data=all.aln, aes(x=gene,fill =sequence.m49_country.display_name )) + geom_bar()+
  ylab("Sequence count")+ labs(fill = "Country")+
  theme(legend.position="top")

#save_plot(file=paste(newfiles, "GeneSeq_by_extraContext.tif", sep="_"), fig = ggplot2::last_plot(), width = 12, height = 9,dpi = 300, theme = ggplot2::theme_get(), label.color = "black",label.size = 2.4, axis.textsize = 0.8, axis.titlesize = 0.75,legend.textsize = 0.4, legend.titlesize = 0.4, legend.itemsize = 0.6)
ggplot(data=subset(all.aln, gene=="wgs" |length.cat=="gene"|length.cat=="pcr"), aes(x=gene,fill =sequence.m49_country.display_name )) + geom_bar()+
  ylab("Sequence count")+ labs(fill = "Country")+
  theme(legend.position="top")+ geom_vline(xintercept=which(all.aln$gene == 'wgs'))+geom_vline(xintercept =1.5, linetype="dotted")

library(ggplot2)
library(ggmap)
library(scatterpie)
library(rworldmap)
library("rnaturalearth")
library("rnaturalearthdata")
world.sf <- ne_countries(scale = "medium", returnclass = "sf")
plot<-
  ggplot(data = world.sf) +
  geom_sf()
plot
seq.countries=unique(all.aln$sequence.m49_country.display_name)
seq.countries=seq.countries[!is.na(seq.countries)]
africa<- world.sf[world.sf$continent == 'Africa' ,]
seq.countries[seq.countries== "Central African Republic"]="Central African Rep."
seq.countries[seq.countries== "Democratic Republic of the Congo"]="Congo"
#seq.c=world[which(world$name %in% seq.countries),]
# highlight countries with seq data:
#ggplot() + geom_sf(data = africa,fill = ifelse(africa$name %in% seq.countries, "yellow", "gray")) 

plot<-
  ggplot(data = africa) +
  geom_sf(col="white", size=0.5) 
plot   
#ignore some problematic areas when extracting centroids: n/ssudan
#africa <- africa[-c(39,40),]
#africa_points <- cbind(africa[-c(39,40),], st_coordinates(st_centroid(africa[-c(39,40),]$geometry)))
africa_points <- cbind(africa[-c(39,40),]$name, st_coordinates(st_centroid(africa[-c(39,40),]$geometry)))
africa_points=as.data.frame(africa_points)
summary=subset(all.aln, gene=="wgs" |length.cat=="gene"|length.cat=="pcr") %>% count(sequence.m49_country.display_name, length.cat)
summary2<- cast(summary, sequence.m49_country.display_name~length.cat, sum)
head(africa_points)
pies=merge(africa_points, summary2, by.x="V1", by.y="sequence.m49_country.display_name")
pies$X=as.numeric(pies$X)
pies$Y=as.numeric(pies$Y)
pies$radius=NA
pies$cat=NA
pies <- pies%>%
  mutate(radius = select(., pcr:wgs) %>% rowSums(na.rm = TRUE))
pies$radius[which(pies$radius==1)]<-0.5
pies$radius[which(pies$radius %in% 2:10)]<-1
pies$radius[which(pies$radius %in% 11:100)]<-2
pies$radius[which(pies$radius %in% 101:500)]<-3
pies$radius[which(pies$radius %in% 501:2000)]<-4
pies$cat[which(pies$radius==0.5)]<- "1"
pies$cat[which(pies$radius==1)]<- "2-10"
pies$cat[which(pies$radius==2)]<- "11-100"
pies$cat[which(pies$radius==3)]<- "101-500"
pies$cat[which(pies$radius==4)]<- "501-2000"
plot+
  geom_scatterpie(aes(x=X, y=Y,group=V1,r=radius),data=pies, cols=c("pcr", "gene","wgs"), color=NA, alpha=.8)+
  theme(panel.grid.major = element_blank())+ 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        #this bit removes the axis lables     
        axis.title.y=element_blank(),
        axis.title.x=element_blank())+
  labs(fill="Sequence type")+
  coord_sf(xlim = c(0,45), ylim = c(-40,15))+
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.6, "in"), pad_y = unit(0.25, "in"),style = north_arrow_fancy_orienteering) +
  scale_fill_viridis(discrete=TRUE)+
geom_scatterpie_legend(radius=pies$radius, n=5, x=3,y=-1,labeller= function(x) x=unique(pies$cat)[c(2,4,1,3,5)])+
  theme(legend.position = c(0.15, 0.4))

# issue: seperate submissions for sequences from same sample (e.g. if g, n, g-l intergenic submitted)
n_occur <- data.frame(table(meta$sequence.isolate))
n_occur[n_occur$Freq > 1,]
duplicates1=meta[meta$sequence.isolate %in% n_occur$Var1[n_occur$Freq > 1],]


# sort out dataset
refine=all.aln;dim(refine)
#refine=refine[refine$length.cat!="wgs",];dim(refine)
# n_occur <- data.frame(table(refine$sequence.isolate))
# n_occur[n_occur$Freq > 1,]
# duplicates=refine[refine$sequence.isolate %in% n_occur$Var1[n_occur$Freq > 1],];dim(duplicates)
# duplicates=duplicates[duplicates$length.cat!="wgs",];dim(duplicates)
# View(duplicates)
rownames(refine)=NULL
refine$gene=as.character(refine$gene)
refine1=aggregate(refine[which( colnames(refine)=="gene" )], refine[-which(colnames(refine)=="gene" )], FUN = function(X) paste(X, collapse=", "))

# aggregate gene records
refine1=refine %>%
  group_by(across(c(-gene))) %>%
  summarise(gene = paste(sort(gene), collapse = ", "))

refine1=refine1 %>%
  mutate(gene.cat=case_when(
    nchar(gene)>1 & length.cat=="gene"~"multi-gene",
    grepl("wgs",gene)~"wgs",
    nchar(gene)>1 & gene=="g, l" ~"multi-gene",
    nchar(gene)==1 & length.cat=="gene" ~gene,
    nchar(gene)==1 & length.cat=="pcr" ~ gene
    ))
levels(refine1$length.cat) <- c(levels(refine1$length.cat),"fragment")
refine1$length.cat[refine1$length.cat=="pcr"]="fragment"
refine1$gene.cat <- factor(refine1$gene.cat,levels = c("wgs","multi-gene", "n", "p", "m","g","l"))
ggplot(data=refine1, aes(x=gene.cat,fill = length.cat)) + geom_bar()+
  ylab("Sequence count")+ labs(fill = "")+scale_fill_viridis(discrete=TRUE)+theme(axis.text.x = element_text(angle = 45))+xlab("Genome region")

ggplot(data=refine1, aes(x=gene.cat,fill =sequence.m49_country.display_name )) + geom_bar()+
  ylab("Sequence count")+ labs(fill = "Country")+
  theme(legend.position="top")+theme(axis.text.x = element_text(angle = 45))+xlab("Genome region")


## find sequence isolate duplicates
n_occur=refine1 %>%
  group_by(sequence.isolate) %>%
  count()
gene.dup=n_occur[which(n_occur$n>1),]
gene.dup=gene.dup[!is.na(gene.dup$sequence.isolate),]

### record linkage between duplicates and why they differ (i.e.duplicate category)
# store output in new dataframe to keep track of steps
refine2=refine1
#have to define new output columns first first to enable use of ifelse in mutate
refine2$linkage=NA
refine2$duplicate.cat=NA
## loop will compare duplicates, record their linkage and why records differ
for (i in 1:nrow(gene.dup)){
    duplicates=refine2[which(refine2$sequence.isolate==gene.dup$sequence.isolate[i]),]
    compare.duplicates=Filter(function(x) any(x != x[1]), duplicates)
   refine2= refine2 %>% 
      #filter(sequence.isolate == gene.dup$Var1[i])%>%
      mutate(linkage=ifelse(sequence.isolate == gene.dup$sequence.isolate[i],paste(compare.duplicates$sequence.sequenceID,collapse=","),linkage), duplicate.cat=ifelse(sequence.isolate == gene.dup$sequence.isolate[i],paste(names(compare.duplicates)[-1],collapse=","),duplicate.cat))
}
# look at why reasons why there are duplicates:
unique(refine2$duplicate.cat)
#common scenarios:
    # 1) Different genes = gene
    # 2) Different seq length= sequence.gb_length
    # 3) Newer sequence = sequence.gb_create_date,sequence.gb_update_date
    # 4) Same sample, sequenced by different groups =sequence.gb_pubmed_id
    # etc.

