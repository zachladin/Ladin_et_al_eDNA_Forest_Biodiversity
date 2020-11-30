#Ladin et al. eDNA analyses

#clear environment
rm(list=ls())

#may need to install phyloseq from source
# source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",local = TRUE)
# install_phyloseq(branch = "github")

#load packages
library(ggmap)
library(ggplot2)
library(ggsn)
library(reshape2)
library(data.table)
library(taxa)
library(metacoder)
library(readxl)
library(pals)
library(cowplot)
library(maps)
library(reshape2)
library(compositions)
library(plyr)
library(phyloseq)
library(scales)
library(vegan)
library(mvabund)
library(magrittr)
library(dplyr)

#summaryFunction
summaryFunction <- function(DataIn, factor, response){
  summaryOut <- ddply(DataIn, factor, .fun = function(xx){
    c(n = length(xx[,response]),
      mean = mean(xx[,response],na.rm=TRUE),
      var = var(xx[,response],na.rm=TRUE),
      SD = sd(xx[,response],na.rm=TRUE),
      SE = sqrt(var(xx[,response],na.rm=TRUE)/length(xx[,response])),
      CV = sd(xx[,response],na.rm=TRUE)/mean(xx[,response],na.rm=TRUE)*100)
  })
  return(summaryOut)
  dev.off()
}

#set working directory
setwd("YOUR WORKING DIRECTORY")

#####################################################################################
#Figure 1. Study area map

#register Google Maps API key
register_google(key="YOUR GOOGLE KEY")

#read in eDNA sample data
sample.data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Data/ASV_All_Data.csv")

#subset to remove Control 2
#sample.data.sub<-subset(sample.data, c(variable !="ZL_2" & variable !="ZL_11"))

#fix ASVdata (find and replace "DI_H20" with "DI_H2O")
levels(sample.data$Sample_type)<-c("DI_H2O","Rain")
levels(sample.data$Location_type)<-c("Control_DI_H2O", "Control_Rain",  "Forest_DI_H2O",  "Forest_Rain")

location.df<-unique(sample.data[,c("Location","point","Longitude","Latitude")])

#fix Location and point 
location.df$point<-ifelse(is.na(location.df$point),"N_of_EW",location.df$point)

location.df$Location<-ifelse(location.df$point=="N_of_EW","Control",location.df$Location)

#add Treatment
location.df$Treatment<-factor(ifelse(location.df$Location=="Control","Control","Forest"), levels=c("Control","Forest"))

#get map center coords
mapCenter<-c(mean(location.df$Longitude,na.rm=TRUE), mean(location.df$Latitude,na.rm=TRUE)-0.0008)

#create map of study area points
ew.map=get_map(location = mapCenter, maptype="satellite",zoom=16)
#ggmap(ew.map)

myColors<-c("tomato","skyblue3")

#get map extent
mapExtent<-attr(ew.map, "bb")

#make point a factor
location.df$point<-as.factor(location.df$point)

location.df$TreatmentName<-ifelse(location.df$Treatment=="Forest","Throughfall","Rainwater")

#modify for map labels
map.location.df<-location.df
map.location.df$point<-gsub("EW.","",map.location.df$point)
map.location.df$point<-ifelse(map.location.df$point=="N_of_EW","Control 1 & 2",map.location.df$point)

#make plot
study.area.map<-ggmap(ew.map, darken=c(0.6,"gray20"))+
  geom_point(data=location.df, aes(x=Longitude,y=Latitude, color=TreatmentName),size=2,alpha=0.9)+
  labs(x="Longitude",y="Latitude")+
  theme(panel.border=element_rect(color="black",fill="transparent"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position = c(0.85,0.85),
        legend.key.size = unit(0.7,"lines"),
        legend.background=element_rect(fill=alpha("gray90",0.8)),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12),
        legend.key=element_rect(color="transparent",fill = "transparent")  )+
  scale_color_manual(values=myColors, name="Treatment")+
  #scale_fill_manual(values=myColors)+
  geom_label(data=head(map.location.df,-1) , aes(x=Longitude, y=Latitude, label=as.character(point)), nudge_y=-0.00035,label.padding=unit(0.3, "lines"), size=3, alpha=0.7)+
  ggsn::scalebar(x.min=mapExtent$ll.lon, x.max=mapExtent$ur.lon, y.min=mapExtent$ll.lat, y.max=mapExtent$ur.lat,transform=TRUE, dist_unit="km",dist=0.1, st.dist=0.017, height = 0.01,anchor=c(x=-75.739,y=39.6578), st.size = 3, st.color = "white")

#view study area map
study.area.map

#now produce map with inset

#get vector map of USA with State Borders
states <- ggplot2::map_data("state")

USmap <- ggplot(data = states,
            mapping = aes(x = long, y = lat,
                          group = group))

USmap.out<-USmap + geom_polygon(fill = "gray60", color = "white")+
  coord_map(xlim = c(-84, -67), ylim=c(25,47))+
  geom_point(aes(x=-75.7, y=39.5), color="red",shape=0, size=3, stroke=1.5)+
  theme(
    plot.background = element_blank(),
    panel.background=element_rect(fill=alpha("white",0.7),color="transparent"),
    panel.grid = element_blank(),
    panel.border=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.title = element_blank()
  )

USmap.out

#now inset eastern US map
gg_inset_map1 = ggdraw() +
  draw_plot(study.area.map) +
  draw_plot(USmap.out, x = 0.14, y = -0.05, width = 0.25, height = 0.7)

#view study area map with inset
gg_inset_map1

#save plot (Figure 1)
ggsave(gg_inset_map1, file="/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Results/Figures/Figure_1_eDNA_Study_Area_color.png",width=8, height=8, dpi=600)

#############################################################################################################
#Prepare ASV data for analysis
count_table_fname <- "/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Data/asv_seqs.count_table.txt"
tax_table_fname <- "/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Data/sintax_taxonomy_predictions.table.txt"
sample_data_fname <- "/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Data/sample_covariates.txt"

#############################################################################################
#ASV count data
# Read data into data.frame
ct <- read.table(file = count_table_fname, sep = "\t", header = TRUE, comment.char = "")

#name first column 'asv'
names(ct)[1] <- "asv"

# Set rownames and drop asv column
ct <- transform(ct, row.names = asv, asv = NULL)

# Sort ASVs by name to match the tax table
ct <- ct[order(rownames(ct)), ]

# Transpose: will have samples for rows, ASVs for columns.  This is like the vegan "community matrix".
ct <- as.data.frame(t(ct))

#add column with Sample_ID
ct.df<-ct
ct.df$Sample_ID<-row.names(ct.df)

#get count of unique ASVs per sample
ct.df$ASV_count<-rowSums(ct.df[,-1] != 0)

#simplify to include Sample_ID and ASV_count
ct.df.sub<-unique(ct.df[,c("Sample_ID","ASV_count")])
ct.df.sub$Sample_ID<-gsub("ZL","Sample",ct.df.sub$Sample_ID)

#Taxonomic assignments
# Read taxonomy table
tt <- read.table(file = tax_table_fname, sep = "\t", header = TRUE)
#transpose
tt <- transform(tt, row.names = asv, asv = NULL)

# Sort ASV by name to match the count table sorting.
tt <- tt[order(rownames(tt)), ]

#Generate Table S1. with Sample covariate information

#read in sample information
sample.info<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Data/Ladin_eDNA_rain_samples_DNAconc.csv")

#remove Location_Sample_ID column
sample.info$Location_Sample_ID<-NULL

sample.info.merge<-merge(sample.info, ct.df.sub, by="Sample_ID",all.x=TRUE)

sample_dat <- read.table(file = sample_data_fname,sep = "\t",header = TRUE)

sample_dat <- transform(sample_dat, row.names = sample, sample = NULL)

# Sort based on rownames and colnames.
sample_dat <- sample_dat[order(rownames(sample_dat)), order(colnames(sample_dat))]

#make a phyloseq object
phy <- phyloseq(
  otu_table(t(ct), taxa_are_rows = TRUE),
  tax_table(as.matrix(tt)),
  sample_data(sample_dat)
)

print(phy)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4133 taxa and 21 samples ]
# sample_data() Sample Data:       [ 21 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 4133 taxa by 7 taxonomic ranks ]

#view taxonomic assignment table
tax_table(phy)

############################################################################################################
#Table S1. Sample information
#now make ASV info into data.frame for saving .csv 
d <- data.frame(
  Sample = sample_names(phy),
  ReadCount = sample_sums(phy),
  Type = sample_dat[sample_names(phy), "type"]
)

#replace 'ZL' with 'Sample
d$Sample<-gsub("ZL","Sample",as.character(d$Sample))
names(d)[names(d)=="Sample"]<-"Sample_ID"

d.out<-d[,c("Sample_ID","ReadCount")]

#now merge with sample info merge
sample.info.out<-merge(sample.info.merge, d.out, by="Sample_ID",all.x=TRUE)

#save as Table 1.
write.csv(sample.info.out, file="/Users/zach/Dropbox (ZachTeam)/Manuscripts/eDNA_Forests_2020/Scientific_Reports/Tables/Table_S1_Sample_Info_11-30-2020.csv", row.names=FALSE)

######################################################################################################
#Figure S2. Rarefaction curve of phyloseq data

otu.df<-otu_table(ct, taxa_are_rows = FALSE)

S <- specnumber(otu.df) # observed number of species
(raremax <- min(rowSums(otu.df)))
Srare <- rarefy(otu.df, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

rarecurve(otu.df, step = 20, col = "blue", cex = 0.6, xlab="Total reads", ylab="Number of ASVs")
abline(v=2000)

#save figure
png(filename="/Users/zach/Dropbox (ZachTeam)/Manuscripts/eDNA_Forests_2020/Scientific_Reports/Figures/Figure_S3_Rarefaction_ASV_samples.png", width=12, height=7, units="in",res=600)

rarecurve(otu.df, step = 20, col = "royalblue", cex = 0.6, xlab="Total reads", ylab="Number of ASVs")
abline(v=2000)

dev.off()
######################################################################################################
#Figure S3. Plot readCount per sample

# Force ordering by sample read count from low to high.
sample.info.out <- sample.info.out[order(sample.info.out$Sample_type, sample.info.out$ReadCount), ]
sample.info.out$Sample_ID <- factor(sample.info.out$Sample_ID, levels=unique(sample.info.out$Sample_ID))

sample.info.out$Location

#aggregate across within sample replicates, sum ReadCounts
sample.readCount<-aggregate(ReadCount~Location+Plot+Location_type+Sample_type, data=sample.info.out, FUN="sum")

unique(sample.readCount$Location)
#set factor levels of Location
sample.readCount$Location<-factor(sample.readCount$Location, levels=c("Control_1","Control_2","C8_1","C8_2","E5_1","E5_2","G7_1","G7_2"))

#rename factor level Control_1 to Control
levels(sample.readCount$Location)<-c("Control_1","Control_2","C8_1","C8_2","E5_1","E5_2","G7_1","G7_2")

#create Sample_Location column
sample.readCount$Sample_Location<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15")

sample.readCount <- sample.readCount[order(sample.readCount$Sample_type, sample.readCount$ReadCount), ]
sample.readCount$Sample_Location <- factor(sample.readCount$Sample_Location, levels=unique(sample.readCount$Sample_Location))

#modify Sample_type
sample.readCount$Sample_type<-ifelse(sample.readCount$Sample_type=="Rain","Precipitation",sample.readCount$Sample_type)

mycolors=c("gray60","royalblue3")

p <- ggplot(data = subset(sample.readCount,Sample_Location !=2), aes(x = Sample_Location, y = ReadCount)) +
  geom_bar(stat="identity",aes(fill=Sample_type),alpha=0.7)+
  scale_fill_manual(values=mycolors)+
  scale_y_continuous(labels = comma) +
  scale_x_discrete(labels=as.character(sample.readCount$Location))+
  geom_hline(yintercept = 2000,
             linetype = "dashed",
             color = "#666666") +
  labs(x="Sample Location", y="Total Reads")+
  # Rotate x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background=element_rect(fill="white",color="black"),
        panel.border=element_rect(fill="transparent",color="black"),
        legend.position=c(0.2,0.8),
        legend.background = element_rect(fill=alpha("gray80",0.3),color="black"))+
    guides(fill=guide_legend(title="Sample Type"))

print(p)

#save plot
ggsave(p, file="/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Results/Figures/Figure_S3_Total_reads_11-30-2020.png", width=14, height=6, dpi=600)
#############################################################################################################
#remove DI H2O ASVs from data

#add Sample_ID to ct
ct <- read.table(file = count_table_fname, sep = "\t", header = TRUE, comment.char = "")

names(ct)[names(ct)=="X.OTU.ID"]<-"ASV"
names(ct)<-gsub("ZL","Sample",names(ct))

#now melt ct
ct.melt<-melt(ct, id.vars = "ASV")

# Sort ASV by name to match the count table sorting.
tt$ASV <- row.names(tt)

#remove genus and species
tt.sub<-tt

#merge with ct.melt
ct.merge<-merge(ct.melt, tt.sub, by="ASV", all.x=TRUE)


#change 'variable' column header to Sample_ID
names(ct.merge)[names(ct.merge)=="variable"]<-"Sample_ID"

#change 'value' column header to 'Abundnace
names(ct.merge)[names(ct.merge)=='value']<-"Count"

#now merge with sample.info
ASV.all<-merge(ct.merge, sample.info.out, by="Sample_ID",all.x=TRUE)

#now aggregate pooling replicates
ASV.all.agg<-aggregate(Count~ASV+Location+Plot+Location_type+Sample_type, data=ASV.all, FUN="sum")

#make ASV_Location_Sample_type column
ASV.all.agg$Sample_type<-ifelse(ASV.all.agg$Sample_type=="Rain","Precipitation",ASV.all.agg$Sample_type)
ASV.all.agg$Location_Sample<-paste(ASV.all.agg$Location, ASV.all.agg$Sample_type, sep="_")

#merge with 
sample.readCount$Location_Sample<-paste(sample.readCount$Location, sample.readCount$Sample_type, sep="_")

ASV.all.merge<-merge(ASV.all.agg, sample.readCount, by=c("Location_Sample","Location","Plot","Location_type","Sample_type"),all.x=TRUE)

####################################################################
#Table S2. ASVs removed from further anlayses: subset ddH2O ASVs to inlcude in Supplementary Materials

#now subset only unique ASVs from within DI_H2O
ASV_DI_H2O_full<-subset(ASV.all.merge, c(Sample_type=="DI_H2O" & Count>0))

#get ASV_DI_H2O_list
ASV_DI_H2O_full_list<-sort(unique(ASV_DI_H2O_full$ASV))

#get list of all ASVs in non-DI_H2O samples
ASV_rain_full<-subset(ASV.all.merge, c(Sample_type!="DI_H2O" & Count>0))

#get ASV_rain_list
ASV_rain_full_list<-unique(ASV_rain_full$ASV)

#get unique list of ASVs in DI_H2O
ASVremoveDIH2OList_full<-sort(setdiff(ASV_DI_H2O_full_list, ASV_rain_full_list))

length(ASVremoveDIH2OList_full)

#get ASV info for table within Supplementary Materials
ct.merge.ASVremove<-subset(ct.merge, c(ASV %in% ASVremoveDIH2OList_full & Count > 0))
ct.merge.ASVremove$Count<-NULL
ct.merge.ASVremove$Sample_ID<-NULL
#get unique ASVs
ct.merge.ASVremove.out<-unique(ct.merge.ASVremove)

#remove any that lack coarsest taxonomic assignment
ct.merge.ASVremove.out<-subset(ct.merge.ASVremove.out, ! is.na(Domain))

#save as .csv
write.csv(ct.merge.ASVremove.out, file="/Users/zach/Dropbox (ZachTeam)/Manuscripts/eDNA_Forests_2020/Scientific_Reports/Tables/Table_S2_ASV_ddH2O_removed.csv",row.names=FALSE)

###################################################################
#subset  only those with > 2000 readCount
ASV.all.merge.sub<-subset(ASV.all.merge, ReadCount>2000)

#unique samples
unique(ASV.all$Sample_ID)

#now subset only unique ASVs from within DI_H2O
ASV_DI_H2O<-subset(ASV.all.merge.sub, c(Sample_type=="DI_H2O" & Count>0))

#get ASV_DI_H2O_list
ASV_DI_H2O_list<-sort(unique(ASV_DI_H2O$ASV))

#get list of all ASVs in non-DI_H2O samples
ASV_rain<-subset(ASV.all.merge.sub, c(Sample_type!="DI_H2O" & Count>0))

#get ASV_rain_list
ASV_rain_list<-unique(ASV_rain$ASV)

#get unique list of ASVs in DI_H2O
ASVremoveDIH2OList<-sort(setdiff(ASV_DI_H2O_list, ASV_rain_list))

#now filter these from ASV.all
ASV.all.keep<-subset(ASV.all.merge.sub, ! ASV  %in% ASVremoveDIH2OList)

#now remove DI_H2O from ASV.all.keep and discard any ASVs with counts <5
ASV.all.keep.sub<-subset(ASV.all.keep, c(Sample_type !="DI_H2O" & Count >4))

#merge with taxa info
ASV.all.keep.sub<-merge(ASV.all.keep.sub, tt.sub, by="ASV",all.x=TRUE)

#remove any contaminants from Control 2 if needed?
Control_2_data<-subset(ASV.all.keep.sub, Location=="Control_2")
Control_2_ASVlist<-unique(Control_2_data$ASV)

Control_2_unique<-setdiff(Control_2_ASVlist, unique(ASV.all.keep.sub$ASV))
#character(0) = no unique ASVs
#############################################################################################
#summarize data for manuscript

#average number of reads
#first aggregate to get totals per pooled sub samples
reads.agg<-aggregate(ReadCount~Location+Plot+Location_type, data=sample.readCount, FUN="sum")
reads.agg.sub<-subset(reads.agg, ReadCount>2000)

min.reads<-min(reads.agg.sub$ReadCount) #4922
max.reads<-max(reads.agg.sub$ReadCount) #196015
mean.reads<-mean(reads.agg.sub$ReadCount) #93200.14

#remove any samples with < 2000 reads
reads.agg.sub<-subset(sample.readCount, ReadCount>2000)

reads.summary<-summaryFunction(DataIn = reads.agg.sub, factor="Location_type",response="ReadCount")

mod.reads<-lm(ReadCount~Location_type, data=reads.agg.sub)
summary(mod.reads) #not significant

#########################################################################3
#average ASV abundance

#remove zeros
ASV.all.keep.sub.noZeros<-subset(ASV.all.keep.sub, Count>5)

#first aggregate to get totals per pooled sub samples
asv.agg<-aggregate(ASV~Location+Plot+Location_type, data=ASV.all.keep.sub.noZeros, FUN="length")
#asv.agg.sub<-subset(asv.agg, ASV_count>500)
length(unique(ASV.all.keep.sub.noZeros$ASV))
min.asv<-min(asv.agg$ASV) #140
max.asv<-max(asv.agg$ASV) #2615
mean.asv<-mean(asv.agg$ASV) #1471

asv.summary<-summaryFunction(DataIn = asv.agg, factor="Location_type",response="ASV")

# Location_type n   mean      var       SD       SE        CV
# 1       Control 2  589.5 404100.5 635.6890 449.5000 107.83528
# 2        Forest 5 1824.0 524713.0 724.3708 323.9485  39.71331

#use non-parametric test to see if number of ASVs differ between Rain and Throughfall
mod.asv<-kruskal.test(ASV~Location_type, data=asv.agg)
mod.asv

#Not significant
# Kruskal-Wallis rank sum test
# 
# data:  ASV by Location_type
# Kruskal-Wallis chi-squared = 2.4, df = 1, p-value = 0.1213
#############################################################################################
#get numbers of unique ASVs

#how many unique ASVs?
ASV.all.noZeros<-subset(ASV.all.keep.sub, Count>5)
length(unique(ASV.all.noZeros$Family))

#group ASVs by Forest, Forest & Control, and Control factors
#get unique list of ASVs in control only
control.data.asv<-subset(ASV.all.noZeros, Location_type=="Control")
controlListAllasv<-unique(control.data.asv$ASV)

forest.data.asv<-subset(ASV.all.noZeros, Location_type=="Forest")
forestListAllasv<-unique(forest.data.asv$ASV)

forestOnlyListASV<-setdiff(forestListAllasv, controlListAllasv)
controlOnlyListASV<-setdiff(controlListAllasv, forestListAllasv)

#create factors
ASV.all.noZeros$Sample_Group<-ifelse(ASV.all.noZeros$ASV %in% forestOnlyListASV, "Throughfall",
                                     ifelse(ASV.all.noZeros$ASV %in% controlOnlyListASV, "Rainwater", "Throughfall & Rainwater"))

#set up factor levels
ASV.all.noZeros$Sample_Group<-factor(ASV.all.noZeros$Sample_Group, levels=rev(c("Rainwater","Throughfall & Rainwater", "Throughfall")))

#step 1, sum across replicates
ASV.all.sum.replicates<-aggregate(Count~ASV+Plot+Location+Location_type+Sample_Group, data=ASV.all.noZeros, FUN="sum")

#pool data by Sample
asv.count.data<-aggregate(ASV~Location+Location_type+Location_Sample+Sample_Group, data=ASV.all.noZeros, FUN="length")
names(asv.count.data)[names(asv.count.data)=="ASV"]<-"ASV_Count"

#get summaries 
asv.count.summary<-summaryFunction(DataIn=asv.count.data, response="ASV_Count", factor="Sample_Group")

# Sample_Group n      mean       var        SD       SE        CV
# 1             Throughfall 5 1170.4000 354195.80 595.14351 266.1563  50.84958
# 2 Throughfall & Rainwater 7  612.1429  71251.81 266.93035 100.8902  43.60589
# 3               Rainwater 2   81.0000   6728.00  82.02439  58.0000 101.26467

#use non-parametric Kruskal-Wallis test due to unbalanced and small sample sizes
mod.1<-kruskal.test(ASV_Count ~ Sample_Group, data = asv.count.data)
mod.1

# Kruskal-Wallis rank sum test
# 
# data:  ASV_Count by Sample_Group
# Kruskal-Wallis chi-squared = 6.8155, df = 2, p-value = 0.03312

#pair-wise comparison

### Dunn test
library(FSA)

DunnPostHoc = dunnTest(ASV_Count ~ Sample_Group, data = asv.count.data, method="bh") 
# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Benjamini-Hochberg method.
# 
# Comparison         Z    P.unadj      P.adj
# 1               Rainwater - Throughfall -2.514286 0.01192738 0.03578213
# 2   Rainwater - Throughfall & Rainwater -1.405528 0.15986420 0.15986420
# 3 Throughfall - Throughfall & Rainwater  1.667986 0.09531853 0.14297780
#############################################################################################################
#Figure 2. Boxplots comparing Throughfall and Rainwater

#define colors for plot
mycolors<-c("tomato","skyblue3")

#add columns showing Throughfall and Rainwater
asv.count.data$Location_type_factor<-ifelse(asv.count.data$Location_type=="Forest","Throughfall","Rainwater")

asvCountBoxPlot<-ggplot(data=asv.count.data, aes(x=Sample_Group, y=ASV_Count))+
  geom_boxplot(aes(fill=Location_type_factor), alpha=0.95, outlier.shape=NA)+
  geom_point(aes(fill=Location_type_factor),color="black",shape=21, alpha = 0.8,
             position = position_jitter(w = 0.05, h = 0))+
  scale_fill_manual(values=mycolors)+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    legend.position="none"
    # legend.background = element_rect(fill="white",color="black"),
    # legend.key=element_blank()
  )+
  labs(x="Sample groups", y="Count of ASVs")+
  guides(fill=guide_legend(title="Treatment"), color="none")+
  annotate("text", x = 1, y=2000, label = "a")+
  annotate("text",x=2, y=1150, label="ab")+
  annotate("text",x=3, y=300, label="b")
asvCountBoxPlot

ggsave(asvCountBoxPlot, file="/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Results/Figures/Figure_2_Number_of_ASVs_per_Location_Type_11-30-2020.png", width=5.5, height=5, dpi=600)

#########################################################################################
#Figure 3. Aitchison distance analysis


#how many unique ASVs?
ASV.all.noZeros<-subset(ASV.all.keep.sub, Count>5)
length(unique(ASV.all.noZeros$Family))

#create factor Sample_Group
ASV.all.noZeros$Sample_Group<-ifelse(ASV.all.noZeros$ASV %in% forestOnlyListASV, "Throughfall",
                                     ifelse(ASV.all.noZeros$ASV %in% controlOnlyListASV, "Rainwater", "Throughfall & Rainwater"))

#set up factor levels
ASV.all.noZeros$Sample_Group<-factor(ASV.all.noZeros$Sample_Group, levels=rev(c("Rainwater","Throughfall & Rainwater", "Throughfall")))

#Compute Aitchison distance = euclidean distance of clr(transformed) Among Sampling locations

#set up community matrix (using ASVs)
ASV.all.sum.replicates.sub<-ASV.all.noZeros[,c("ASV", "Count","Plot","Location" ,"Location_type","Sample_Group","Sample_Location")]

#cast data into wide format
ASV_Location_cast<-data.table::dcast(Location+Location_type+Sample_Location~ASV, value.var="Count", data=ASV.all.sum.replicates.sub, fun.aggregate = sum, fill=0)

unique(ASV_Location_cast$Location)

#reorder factors
ASV_Location_cast$Location<-factor(ASV_Location_cast$Location, levels=c("C8_1","E5_1","E5_2","G7_1","G7_2","Control_1","Control_2"))

#reorder data
ASV_Location_cast<-ASV_Location_cast[order(ASV_Location_cast$Location),]

levels(ASV_Location_cast$Location)

#Try and impute zeros using other methods
library(zCompositions)

ASV_Location_cast_values<-ASV_Location_cast[,-1:-3]

ASV_Location_ReplZerosGB<-cmultRepl(ASV_Location_cast_values,method="SQ", output="p-counts")

#prepare for compositional analysis add 1 for imputation of zeros
ASV.acomp<-acomp(ASV_Location_ReplZerosGB)

#centered log-ratio transform
ASV.clr<-clr(ASV.acomp)

#get matrix of Aitchison distances
ASV.distance.location.matrix<-as.data.frame(as.matrix(dist(ASV.clr, method="euclidean")))

names(ASV.distance.location.matrix)<-unique(as.character(ASV_Location_cast$Location))

#get lower triangle
ASV.distance.location.matrix[lower.tri(ASV.distance.location.matrix,diag=FALSE)] <- NA

rownames(ASV.distance.location.matrix)<-c(3,4,5,6,7,1,2)

#melt matrix (now lower triangle)
ASV.distance.location.melt<-melt(as.matrix(ASV.distance.location.matrix),na.rm=TRUE)

levels(ASV.distance.location.melt$Var2)

#set up factor levels
ASV.distance.location.melt$Var2<-factor(ASV.distance.location.melt$Var2, levels=rev(levels(ASV.distance.location.melt$Var2)))
levels(ASV.distance.location.melt$Var2)<-c("Control_2","Control_1","G7_2","G7_1","E5_2","E5_1","C8_1")

unique(ASV.distance.location.melt$Var1)
ASV.distance.location.melt$Var1<-factor(ASV.distance.location.melt$Var1,levels=c(3,4,5,6,7,1,2))

levels(ASV.distance.location.melt$Var1)<-c("C8_1","E5_1","E5_2","G7_1","G7_2","Control_1","Control_2")

range(ASV.distance.location.melt$value)

#make a matrix location
distMatrixLocation<-ggplot(data=ASV.distance.location.melt, aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(low="gray90",high="seagreen4",limits=c(0, 272), breaks=c(0,68,136,204,272))+
  geom_text(aes(label=round(value,1)))+
  labs(x="Sample",y="Sample")+
  guides(fill=guide_colorbar(title="Aitchison\ndistance"))+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    axis.title = element_blank(),
    legend.position=c(0.85,0.7),
    legend.key.size=unit(2, "lines"),
    legend.text = element_text(size=12),
    legend.title = element_text(size=14),
    #legend.background = element_rect(fill="transparent",color="gray10")
    #plot.background = element_rect(fill="transparent",color="black")
  )+
  ggtitle("A) Sampling Locations")
distMatrixLocation

#get summaries
ASV.distance.location.melt.noZeros<-subset(ASV.distance.location.melt, value>0)
range(ASV.distance.location.melt.noZeros$value)#169.3563 271.4427

#Compte Aitchison distance = euclidean distance of clr(transformed) among plots

#set up community matrix (using ASVs)
ASV.all.sum.replicates.sub<-ASV.all.sum.replicates[,c("ASV", "Count","Plot","Location" ,"Location_type","Sample_Group")]

#cast data into wide format
ASV_Plot_cast<-dcast(Plot+Location_type~ASV, value.var="Count", data=ASV.all.sum.replicates.sub, fun.aggregate = sum, fill=0)
#reorder factors
ASV_Plot_cast$Plot<-factor(ASV_Plot_cast$Plot, levels=c("C8","E5","G7","N_of_EW"))
#reorder data
ASV_Plot_cast<-ASV_Plot_cast[order(ASV_Plot_cast$Plot),]

ASV_Plot_cast_values<-ASV_Plot_cast[,-1:-2]

ASV_Plot_ReplZerosGB<-cmultRepl(ASV_Plot_cast_values,method="SQ", output="p-counts")

#prepare for compositional analysis add 1 for imputation of zeros
ASV.acomp<-acomp(ASV_Plot_ReplZerosGB)

#prepare for compositional analysis add 1 for imputation of zeros
#ASV.acomp<-acomp(ASV_Plot_cast[,-1:-2]+1)

#centered log-ratio transform
ASV.clr<-clr(ASV.acomp)

levels(ASV_Plot_cast$Plot)

#get matrix of Aitchison distances
ASV.distance.plot.matrix<-as.data.frame(as.matrix(dist(ASV.clr, method="euclidean")))

names(ASV.distance.plot.matrix)<-unique(as.character(ASV_Plot_cast$Plot))

#get lower triangle
ASV.distance.plot.matrix[lower.tri(ASV.distance.plot.matrix,diag=FALSE)] <- NA

#melt matrix (now lower triangle)
ASV.distance.plot.melt<-melt(as.matrix(ASV.distance.plot.matrix),na.rm=TRUE)

#set up factor levels
ASV.distance.plot.melt$Var2<-factor(ASV.distance.plot.melt$Var2, levels=rev(levels(ASV.distance.plot.melt$Var2)))
levels(ASV.distance.plot.melt$Var2)<-c("Control","G7","E5","C8")

ASV.distance.plot.melt$Var1<-factor(ASV.distance.plot.melt$Var1,levels=rev(c(4,3,2,1)))
levels(ASV.distance.plot.melt$Var1)<-c("C8","E5","G7","Control")

#make a matrix plot
distMatrixPlot<-ggplot(data=ASV.distance.plot.melt, aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(low="gray90",high="seagreen4",limits=c(0, 272))+
  geom_text(aes(label=round(value,1)))+
  labs(x="Plot",y="Plot")+
  guides(fill=guide_colorbar(title="Aitchison\ndistance"))+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    axis.title =element_blank(),
    legend.position="none"
  )+
  ggtitle("B) Plots")
distMatrixPlot

#get summaries
ASV.distance.plot.melt.noZeros<-subset(ASV.distance.plot.melt, value>0)
range(ASV.distance.plot.melt.noZeros$value)#228.0017 252.0726

#now compare Forest and Control Treatments Aitchison Distances

#set up community matrix (using ASVs)
ASV.all.sum.replicates.sub<-ASV.all.noZeros[,c("ASV", "Count","Plot","Location" ,"Location_type","Sample_Group")]

#cast data into wide format
ASV_Treatment_cast<-dcast(Location_type~ASV, value.var="Count", data=ASV.all.sum.replicates.sub, fun.aggregate = sum, fill=0)

#reorder factors
ASV_Treatment_cast$Location_type<-factor(ASV_Treatment_cast$Location_type, levels=c("Control","Forest"))

#reorder data
ASV_Treatment_cast<-ASV_Treatment_cast[order(ASV_Treatment_cast$Location_type),]

ASV_Treatment_cast_values<-ASV_Treatment_cast[,-1]

ASV_Treatment_ReplZerosGB<-cmultRepl(ASV_Treatment_cast_values,method="SQ", output="p-counts")

#prepare for compositional analysis add 1 for imputation of zeros
ASV.acomp<-acomp(ASV_Treatment_ReplZerosGB)

#centered log-ratio transform
ASV.clr<-clr(ASV.acomp)

levels(ASV_Treatment_cast$Location_type)

#get matrix of Aitchison distances
ASV.distance.treatment.matrix<-as.data.frame(as.matrix(dist(ASV.clr, method="euclidean")))

names(ASV.distance.treatment.matrix)<-unique(as.character(ASV_Treatment_cast$Location_type))

#get lower triangle
ASV.distance.treatment.matrix[lower.tri(ASV.distance.treatment.matrix,diag=FALSE)] <- NA

#melt matrix (now lower triangle)
ASV.distance.treatment.melt<-melt(as.matrix(ASV.distance.treatment.matrix),na.rm=TRUE)

#set up factor levels
ASV.distance.treatment.melt$Var2<-factor(ASV.distance.treatment.melt$Var2, levels=rev(levels(ASV.distance.treatment.melt$Var2)))
levels(ASV.distance.treatment.melt$Var2)<-c("Control","Forest")

ASV.distance.treatment.melt$Var1<-factor(ASV.distance.treatment.melt$Var1,levels=rev(c(2,1)))
levels(ASV.distance.treatment.melt$Var1)<-c("Forest","Control")

#switch to Throughfall and Rainwater
levels(ASV.distance.treatment.melt$Var1)<-c("Throughfall","Rainwater")
levels(ASV.distance.treatment.melt$Var2)<-c("Rainwater","Throughfall")

#make a matrix treatment
distMatrixTreatment<-ggplot(data=ASV.distance.treatment.melt, aes(x=Var1, y=Var2, fill=value))+
  geom_tile()+
  scale_fill_gradient(low="gray90",high="seagreen4",limits=c(0, 272))+
  geom_text(aes(label=round(value,1)))+
  labs(x="Treatment",y="Treatment")+
  guides(fill=guide_colorbar(title="Aitchison\ndistance"))+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    axis.title=element_blank(),
    legend.position = "none"
  )+
  ggtitle("C) Treatments")
distMatrixTreatment

#box plot of Aitchison distances (Treatment)

#set up data.frame
ASV.distance.plot.melt$Location_type_diff<-ifelse(
  c(substr(ASV.distance.plot.melt$Var1,1,7)=="Control" & substr(ASV.distance.plot.melt$Var2,1,7)=="Control"),"Control_Control",
  ifelse(c(substr(ASV.distance.plot.melt$Var1,1,7)=="Control" & substr(ASV.distance.plot.melt$Var2,1,7)!="Control"),"Forest_Control",
         ifelse(c(substr(ASV.distance.plot.melt$Var1,1,7)!="Control" & substr(ASV.distance.plot.melt$Var2,1,7)=="Control"),"Forest_Control","Forest_Forest")))

#now remove duplicates
ASV.distance.plot.melt.sub<-unique(ASV.distance.plot.melt[,c("value","Location_type_diff")])
names(ASV.distance.plot.melt.sub)[names(ASV.distance.plot.melt.sub)=="value"]<-"Aitchison_Distance"

#remove zeros (self-similar comparisons)
ASV.distance.plot.melt.sub<-subset(ASV.distance.plot.melt.sub, Aitchison_Distance != 0)

unique(ASV.distance.plot.melt.sub$Location_type_diff)
ASV.distance.plot.melt.sub$Location_type_diff<-factor(ASV.distance.plot.melt.sub$Location_type_diff, levels=c("Forest_Forest","Forest_Control","Control_Control"))

mod.2<-kruskal.test(Aitchison_Distance ~ Location_type_diff, data = ASV.distance.plot.melt.sub)
mod.2

# Kruskal-Wallis rank sum test
# 
# data:  Aitchison_Distance by Location_type_diff
# Kruskal-Wallis chi-squared = 0.047619, df = 1, p-value = 0.8273

#try with location
#set up data.frame
ASV.distance.location.melt$Location_type_diff<-ifelse(
  c(substr(ASV.distance.location.melt$Var1,1,7)=="Control" & substr(ASV.distance.location.melt$Var2,1,7)=="Control"),"Control_Control",
  ifelse(c(substr(ASV.distance.location.melt$Var1,1,7)=="Control" & substr(ASV.distance.location.melt$Var2,1,7)!="Control"),"Forest_Control",
         ifelse(c(substr(ASV.distance.location.melt$Var1,1,7)!="Control" & substr(ASV.distance.location.melt$Var2,1,7)=="Control"),"Forest_Control",
                "Forest_Forest")))

#now remove duplicates
ASV.distance.location.melt.sub<-unique(ASV.distance.location.melt[,c("value","Location_type_diff")])
names(ASV.distance.location.melt.sub)[names(ASV.distance.location.melt.sub)=="value"]<-"Aitchison_Distance"

#remove zeros (self-similar comparisons)
ASV.distance.location.melt.sub<-subset(ASV.distance.location.melt.sub, Aitchison_Distance != 0)

unique(ASV.distance.location.melt.sub$Location_type_diff)
ASV.distance.location.melt.sub$Location_type_diff<-factor(ASV.distance.location.melt.sub$Location_type_diff, levels=c("Forest_Forest","Forest_Control","Control_Control"))

mod.3<-kruskal.test(Aitchison_Distance ~ Location_type_diff, data = ASV.distance.location.melt.sub)
mod.3
# Kruskal-Wallis rank sum test
# 
# data:  Aitchison_Distance by Location_type_diff
# Kruskal-Wallis chi-squared = 3.6208, df = 2, p-value = 0.1636

#fill box plot as function of meadian value
box.data = ASV.distance.location.melt.sub %>% mutate(Location_type_diff = Location_type_diff) %>%
  group_by(Location_type_diff) %>% 
  mutate(medianValue = median(Aitchison_Distance))

#merge with data.frame
median.df<-data.frame(Location_type_diff=box.data$Location_type_diff, Median=box.data$medianValue)
ASV.distance.location.melt.sub.merge<-merge(ASV.distance.location.melt.sub, unique(median.df), by="Location_type_diff",all.x=TRUE)

levels(ASV.distance.location.melt.sub.merge$Location_type_diff)<-c("Throughfall-Throughfall","Throughfall-Rainwater","Rainwater-Rainwater")


#plot boxplot
AitchisonDistBoxplotPlot<-ggplot(data=subset(ASV.distance.location.melt.sub.merge, Location_type_diff!="Rainwater-Rainwater"), aes(x=Location_type_diff, y=Aitchison_Distance))+
  geom_boxplot(aes(fill=Median),alpha=0.8,outlier.shape=NA)+
  geom_point(aes(fill=Median), color="black",shape=21,alpha=0.8,
             position = position_jitter(w = 0.05, h = 0))+
  scale_fill_gradient(low="gray90",high="seagreen4",limits=c(0, 272))+
  theme(
    panel.background = element_rect(fill="white", color="black"),
    panel.border = element_rect(fill="transparent",color="black"),
    legend.position = "none"
  )+
  #ylim(c(0,60))+
  labs(x="Treatment Comparison", y="Aitchison Distance")+
  ggtitle("D) Community dissimilarity")

AitchisonDistBoxplotPlot

#save combined plot
library(ggpubr)

matrixPlotsCombine<-ggarrange(distMatrixLocation, distMatrixPlot, distMatrixTreatment, AitchisonDistBoxplotPlot, common.legend = FALSE)
matrixPlotsCombine

#save distance matrix
ggsave(matrixPlotsCombine, file="/Users/zach/Dropbox (ZachTeam)/Manuscripts/eDNA_Forests_2020/Scientific_Reports/Figures/Figure_3_Aitchison_distance_matrices_plot_and_Boxplot_11-30-2020.png", width=11, height=11, dpi=600)

#########################################################################################
#Examine alpha and beta diversity patterns with DivNet package

# remotes::install_github("adw96/breakaway")
# remotes::install_github("adw96/DivNet")
library(DivNet) #version 0.3.6

#Estimate Alpha and Beta Diveristy with DivNet

#how many unique ASVs?
ASV.all.noZeros<-subset(ASV.all.keep.sub, Count>5)

ASV.all.sum.replicates.sub<-ASV.all.noZeros[,c("ASV", "Count","Plot","Location" ,"Location_type")]

ASV_Location_cast<-dcast(Location+Location_type~ASV, value.var="Count", data=ASV.all.sum.replicates.sub, fun.aggregate = sum, fill=0)

#reorder factors
ASV_Location_cast$Location<-factor(ASV_Location_cast$Location, levels=c("C8_1","E5_1","E5_2","G7_1","G7_2","Control_1","Control_2"))

#reorder data
ASV_Location_cast<-ASV_Location_cast[order(ASV_Location_cast$Location),]

levels(ASV_Location_cast$Location)

ASV_Location_cast_values<-ASV_Location_cast[,-1:-3]

#Try and impute zeros using other methods
library(zCompositions)

#Impute zeros
ASV_Location_ReplZerosGB<-cmultRepl(ASV_Location_cast_values,method="SQ", output="p-counts")
ASV_Location_ReplZerosGB$Location<-ASV_Location_cast$Location

#now prepare for phyloseq object
otu.data<-ASV_Location_ReplZerosGB
row.names(otu.data)<-otu.data$Location
otu.data$Location<-NULL

#taxa table
tt.data<-unique(data.frame(ASV=ASV.all.noZeros$ASV, ASV.all.noZeros[,10:16]))
row.names(tt.data)<-tt.data$ASV
tt.data$ASV<-NULL
#sample covariates
sample.data<-unique(ASV.all.noZeros[,c("Location","Plot","Location_type")])
row.names(sample.data)<-sample.data$Location

#make columns factors
sample.data$Location<-as.factor(sample.data$Location)
sample.data$Plot<-as.factor(sample.data$Plot)
sample.data$Location_type<-as.factor(sample.data$Location_type)


#make a phyloseq object
phy.new <- phyloseq(
  otu_table(t(otu.data), taxa_are_rows = TRUE),
  tax_table(as.matrix(tt.data)),
  sample_data(sample.data)
)

#treatment Throughfall vs. Rain
divnet_treatment <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                         X = "Location_type",
                         tuning = "careful",
                         ncores = 3)


alphaDiv_treatment<-testDiversity(divnet_treatment, h0="shannon")
alphaDiv_treatment
# Estimates Standard Errors p-values
# (Intercept)      3.1402762      0.03596808        0
# predictorsForest 0.5423623      0.03845448        0

estimates_treatment <- divnet_treatment$shannon %>% summary %$% estimate
ses <- sqrt(divnet_treatment$`shannon-variance`)
X <- breakaway::make_design_matrix(phy.new, "Location_type")
betta(estimates_treatment, ses, X)$table
betta(estimates_treatment, ses, X)$global

#gather estimates
alphaDiv_treatment_summary<-data.frame(divnet_treatment$shannon %>% summary)
#simplify
alphaDiv_treatment_summary_sub<-subset(alphaDiv_treatment_summary, c(sample_names =="C8_1" | sample_names=="Control_1"))
#add Treatment
alphaDiv_treatment_summary_sub$Treatment<-factor(ifelse(alphaDiv_treatment_summary_sub$sample_names=="C8_1","Throughfall","Rainwater"),levels=c("Throughfall","Rainwater"))

# estimate      error    lower    upper sample_names   name     model   Treatment
# 1 3.682638 0.08598683 3.510665 3.854612         C8_1 DivNet Aitchison Throughfall
# 6 3.140276 0.14379424 2.852688 3.427865    Control_1 DivNet Aitchison   Rainwater

#now plot
alphaDiv_treatmentFig<-ggplot(data=alphaDiv_treatment_summary_sub, aes(x=Treatment, y=estimate))+
  geom_errorbar(aes(ymin=estimate-error, ymax=estimate+error,color=Treatment),width=0, size=1.3)+
  geom_point(aes(color=Treatment),size=3)+
  scale_color_manual(values=rev(myColors))+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    legend.position=c(0.8,0.8),
    legend.background = element_rect(fill="white",color="black"),
    legend.key=element_blank()
  )+
  labs(x="Treatment", y="Shannon's Diversity Estimate")+
  guides(color=guide_legend(title="Treatment"))+
  annotate("text", x = 1, y=3.83, label = "a")+
  annotate("text",x=2, y=3.34, label="b")+
  ggtitle("A) Alpha Diversity Between Treatments")
alphaDiv_treatmentFig

###############################################################
#Alpha diversity among Plots

#Test differences between C8
phy.new@sam_data$Plot<-factor(phy.new@sam_data$Plot, levels=c("C8","E5","G7","N_of_EW"))
divnet_plot_C8 <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                            tuning="careful", 
                            X = "Plot",
                            ncores = 3)

alphaDiv_plot_C8<-testDiversity(divnet_plot_C8, h0="shannon")
alphaDiv_plot_C8
# Hypothesis testing:
#   p-value for global test: 0 
# > alphaDiv_plot
# Estimates Standard Errors p-values
# (Intercept)        3.9630091      0.06027797    0.000
# predictorsE5      -0.1198057      0.08083718    0.138
# predictorsG7      -0.2985582      0.12865543    0.020
# predictorsN_of_EW -0.8245994      0.15127298    0.000

#get summary of estimates
estimates_C8 <- divnet_plot_C8$shannon %>% summary
#     estimate error lower upper sample_names name   model    
# 1     3.96 0.164  3.64  4.29 C8         DivNet Aitchison
# 2     3.84 0.147  3.54  4.13 E5         DivNet Aitchison
# 4     3.65 0.249  3.16  4.15 G7         DivNet Aitchison
# 7     3.14 0.302  2.54  3.74 Control    DivNet Aitchison

estimates_C8 <- divnet_plot_C8$shannon %>% summary %$% estimate
ses <- sqrt(divnet_plot_C8$`shannon-variance`)
X <- breakaway::make_design_matrix(phy.new, "Plot")
betta(estimates_C8, ses, X)$table
betta(estimates_C8, ses, X)$global
# 3785.884    0.000

#Test differences between E5
phy.new@sam_data$Plot<-factor(phy.new@sam_data$Plot, levels=c("E5","C8","G7","N_of_EW"))
divnet_plot_E5 <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                          tuning="careful", 
                          X = "Plot",
                          ncores = 3)

alphaDiv_plot_E5<-testDiversity(divnet_plot_E5, h0="shannon")
alphaDiv_plot_E5
# Hypothesis testing:
#   p-value for global test: 0 
# > alphaDiv_plot
# Estimates Standard Errors p-values
# (Intercept)        3.8304914      0.06158367    0.000
# predictorsC8       0.1325177      0.17006597    0.436
# predictorsG7      -0.1718159      0.11207120    0.125
# predictorsN_of_EW -0.6920638      0.20135236    0.001

estimates_E5 <- divnet_plot_E5$shannon %>% summary %$% estimate
ses <- sqrt(divnet_plot_E5$`shannon-variance`)
X <- breakaway::make_design_matrix(phy.new, "Plot")
betta(estimates_E5, ses, X)$table
betta(estimates_E5, ses, X)$global
#3683.114    0.000

#Test differences between G7
phy.new@sam_data$Plot<-factor(phy.new@sam_data$Plot, levels=c("G7","C8","E5","N_of_EW"))
divnet_plot_G7 <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                          tuning="careful",
                          X = "Plot",
                          ncores = 3)

alphaDiv_plot_G7<-testDiversity(divnet_plot_G7, h0="shannon")
alphaDiv_plot_G7
# Hypothesis testing:
#   p-value for global test: 0 
# > alphaDiv_plot
# Estimates Standard Errors p-values
# Estimates Standard Errors p-values
# (Intercept)        3.6639598      0.05801467    0.000
# predictorsC8       0.2990493      0.15705915    0.057
# predictorsE5       0.1728746      0.10265382    0.092
# predictorsN_of_EW -0.5427089      0.21620184    0.012

estimates_G7 <- divnet_plot_G7$shannon %>% summary %$% estimate
ses <- sqrt(divnet_plot_G7$`shannon-variance`)
X <- breakaway::make_design_matrix(phy.new, "Plot")
betta(estimates_G7, ses, X)$table
betta(estimates_G7, ses, X)$global
#4125.384    0.000


#Test differences between Control
phy.new@sam_data$Plot<-factor(phy.new@sam_data$Plot, levels=c("N_of_EW","C8","E5","G7"))
divnet_plot_Control <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                          tuning = "careful",
                          X = "Plot",
                          ncores = 3)

alphaDiv_plot_Control<-testDiversity(divnet_plot_Control, h0="shannon")
alphaDiv_plot_Control
# Hypothesis testing:
#   p-value for global test: 0 
# > alphaDiv_plot
# Estimates Standard Errors p-values
# (Intercept)  3.1367930      0.05543939        0
# predictorsC8 0.8262161      0.18357651        0
# predictorsE5 0.6936129      0.10865222        0
# predictorsG7 0.5220343      0.11101540        0


estimates_Control <- divnet_plot_Control$shannon %>% summary %$% estimate
ses <- sqrt(divnet_plot_Control$`shannon-variance`)
X <- breakaway::make_design_matrix(phy.new, "Plot")
betta(estimates_Control, ses, X)$table
betta(estimates_Control, ses, X)$global
#4072.625    0.000


#gather estimates
alphaDiv_plot_summary<-data.frame(divnet_plot_C8$shannon %>% summary)

#add Plot name
alphaDiv_plot_summary$Plot<-factor(c("C8","E5","E5","G7","G7","Control","Control"),levels=c("C8","E5","G7","Control"))

#simplify
alphaDiv_plot_summary$sample_names<-NULL
alphaDiv_plot_summary_sub<-unique(alphaDiv_plot_summary)

#add Treatment
alphaDiv_plot_summary_sub$Treatment<-factor(c("Throughfall","Throughfall","Throughfall","Rainwater"),levels=c("Throughfall","Rainwater"))

#now plot
alphaDiv_plotFig<-ggplot(data=alphaDiv_plot_summary_sub, aes(x=Plot, y=estimate))+
  geom_errorbar(aes(ymin=estimate-error, ymax=estimate+error,color=Treatment),width=0, size=1.3)+
  geom_point(aes(color=Treatment),size=3)+
  scale_color_manual(values=rev(myColors))+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    legend.position="none"
    # legend.background = element_rect(fill="white",color="black"),
    # legend.key=element_blank()
  )+
  labs(x="Plot", y="Shannon's Diversity Estimate")+
  guides(fill=guide_legend(title="Plot"), color="none")+
  annotate("text", x = 1, y=4.26, label = "a")+
  annotate("text",x=2, y=4.02, label="ab")+
  annotate("text",x=3, y=3.95, label="bc")+
  annotate("text",x=4, y=3.42, label="d")+
  ggtitle("B) Alpha Diversity Among Plots")
    
alphaDiv_plotFig

#alpha Diversity among sampling locations

#C8_1
phy.new@sam_data$Location<-factor(phy.new@sam_data$Location, levels=c("C8_1","E5_1","E5_2","G7_1","G7_2","Control_1", "Control_2"))
divnet_location_C8_1 <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                       X = "Location",
                       tuning="careful",
                       ncores = 3)

alphaDiv_location_1_C8<-testDiversity(divnet_location_C8_1, h0="shannon")
alphaDiv_location_1_C8
# Hypothesis testing:
#   p-value for global test: 0 
#> alphaDiv_location_1_C8
# Estimates Standard Errors p-values
# (Intercept)          3.963009064       0.1220285    0.000
# predictorsE5_1      -0.196036408       0.3226284    0.543
# predictorsE5_2       0.006820705       0.3226714    0.983
# predictorsG7_1      -0.330446136       0.3231464    0.307
# predictorsG7_2       0.002886549       0.3228144    0.993
# predictorsControl_1 -0.880239694       0.3234101    0.006
# predictorsControl_2 -0.054445440       0.3226709    0.866

estimates_C8_1 <- divnet_location_C8_1$shannon %>% summary %$% estimate
ses <- sqrt(divnet_location_C8_1$`shannon-variance`)
X <- breakaway::make_design_matrix(phy.new, "Location")
betta(estimates_C8_1, ses, X)$global
#953.4297   0.0000

#E5_1
phy.new@sam_data$Location<-factor(phy.new@sam_data$Location, levels=c("E5_1","E5_2","G7_1","G7_2","Control_1", "Control_2"))
divnet_location_E5_1 <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                                X = "Location",
                                tuning="careful",
                                ncores = 3)

alphaDiv_location_1_E5<-testDiversity(divnet_location_E5_1, h0="shannon")
alphaDiv_location_1_E5
# Hypothesis testing:
#   p-value for global test: 0 
#> alphaDiv_location_1_E5
# Estimates Standard Errors p-values
# (Intercept)          3.963009064       0.1220285    0.000
# predictorsE5_1      -0.196036408       0.3226284    0.543
# predictorsE5_2       0.006820705       0.3226714    0.983
# predictorsG7_1      -0.330446136       0.3231464    0.307
# predictorsG7_2       0.002886549       0.3228144    0.993
# predictorsControl_1 -0.880239694       0.3234101    0.006
# predictorsControl_2 -0.054445440       0.3226709    0.866

estimates_E5_1 <- divnet_location_E5_1$shannon %>% summary %$% estimate
ses <- sqrt(divnet_location_E5_1$`shannon-variance`)
X <- breakaway::make_design_matrix(phy.new, "Location")
betta(estimates_E5_1, ses, X)$global
#953.4297   0.0000


#########################################################################################################
#Alpha diversity plot among locations
phy.new@sam_data$Location<-factor(phy.new@sam_data$Location, levels=c("C8_1","E5_1","E5_2","G7_1","G7_2","Control_1", "Control_2"))
divnet_location <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                                X = "Location",
                                tuning="careful",
                                ncores = 3)

#gather estimates
alphaDiv_location_summary<-data.frame(divnet_location$shannon %>% summary)

alphaDiv_location_summary$Location<-factor(alphaDiv_location_summary$sample_names, levels=c("C8_1","E5_1","E5_2","G7_1","G7_2","Control_1", "Control_2"))

#add Treatment
alphaDiv_location_summary$Treatment<-factor(c("Throughfall","Throughfall","Throughfall","Throughfall","Throughfall","Rainwater","Rainwater"),levels=c("Throughfall","Rainwater"))

#now plot
alphaDiv_locationFig<-ggplot(data=alphaDiv_location_summary, aes(x=Location, y=estimate))+
  geom_errorbar(aes(ymin=estimate-error, ymax=estimate+error,color=Treatment),width=0, size=1.3)+
  geom_point(aes(color=Treatment),size=2)+
  scale_color_manual(values=rev(myColors))+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    legend.position="none"
    # legend.background = element_rect(fill="white",color="black"),
    # legend.key=element_blank()
  )+
  labs(x="Plot", y="Shannon's Diversity Estimate")+
  guides(fill=guide_legend(title="Sampling Location"), color="none")+
  ylim(0,4)
  # annotate("text", x = 1, y=4.25, label = "a")+
  # annotate("text",x=2, y=3.95, label="ab")+
  # annotate("text",x=3, y=3.88, label="bc")+
  # annotate("text",x=4, y=3.45, label="c")
alphaDiv_locationFig

###############################################################
#Beta diversity with DivNet

#need to change column headers to compare Throughfall-Throughfall and Throughfall-Control beta diversity
row.names(otu.data)

betaDiv_treatment<-simplifyBeta(divnet_treatment, phy.new, measure="bray-curtis", x="Location_type")
betaDiv_treatment
# Covar1 Covar2  beta_est  beta_var     lower     upper
# 1 Control Forest 0.3870081 0.0100319 0.1866894 0.5873268

######################################################################
#Plot-level
phy.new@sam_data$Plot<-factor(phy.new@sam_data$Plot, levels=c("C8","E5","G7","N_of_EW"))
divnet_plot <-  divnet(tax_glom(phy.new, taxrank="Genus"),
                       X = "Plot",
                       tuning="careful",
                       ncores = 3)

betaDiv_plot<-simplifyBeta(divnet_plot, phy.new, measure="bray-curtis", x="Plot")
betaDiv_plot
# Covar1 Covar2  beta_est    beta_var     lower     upper
# 1      E5     C8 0.3333768 0.030404352 0.0000000 0.6821137
# 2      G7     C8 0.5040055 0.001749873 0.4203425 0.5876684
# 3 N_of_EW     C8 0.4374381 0.014371035 0.1976796 0.6771966
# 4      G7     E5 0.4084018 0.018354882 0.1374414 0.6793622
# 5 N_of_EW     E5 0.4463761 0.012898584 0.2192322 0.6735200
# 6 N_of_EW     G7 0.4823235 0.007198741 0.3126327 0.6520143

betaDiv_plot$Covs<-c("C8-E5","C8-G7","C8-Control","E5-G7","E5-Control","G7-Control","Control-Control")
betaDiv_plot$Covs<-factor(betaDiv_plot$Covs, levels=c("C8-G7","E5-G7","C8-E5","C8-Control","E5-Control","G7-Control","Control-Control"))

betaDiv_plot$Location_Type_Compare<-as.factor(c("Throughfall_Throughfall", "Throughfall_Throughfall","Throughfall_Rainwater", "Throughfall_Throughfall","Throughfall_Rainwater","Throughfall_Rainwater","Rainwater_Rainwater"))

#Hypothesis testing

# H0: BC(E5-C8) = BC(G7-C8))
betta(betaDiv_plot %>% pull(beta_est) %>% .[c(1,2)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(1,2)], 
      cbind("E5-C8" = c(1,1), "G7-C8" = c(0,1)))$table
# Estimates Standard Errors p-values
# E5-C8 0.3333070      0.09610084    0.001
# G7-C8 0.1692416      0.12185441    0.165

# H0: BC(E5-C8) = BC(G7-E5))
betta(betaDiv_plot %>% pull(beta_est) %>% .[c(1,3)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(1,3)], 
      cbind("E5-C8" = c(1,1), "Control-C8" = c(0,1)))$table
# Estimates Standard Errors p-values
# E5-C8      0.3333070       0.0901954    0.000
# Control-C8 0.1011799       0.1321601    0.444

# H0: BC(E5-C8) = BC(G7-E5))
betta(betaDiv_plot %>% pull(beta_est) %>% .[c(1,4)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(1,4)], 
      cbind("E5-C8" = c(1,1), "G7-E5" = c(0,1)))$table
# Estimates Standard Errors p-values
# E5-C8 0.33330697      0.08013867    0.000
# G7-E5 0.07216009      0.11392507    0.526

# H0: BC(E5-C8) = BC(Control-E5))
betta(betaDiv_plot %>% pull(beta_est) %>% .[c(1,5)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(1,5)], 
      cbind("E5-C8" = c(1,1), "Control-E5" = c(0,1)))$table
# Estimates Standard Errors p-values
# E5-C8      0.3333768       0.1124793    0.003
# Control-E5 0.1129992       0.1388632    0.416

# H0: BC(E5-C8) = BC(Control-G7))
betta(betaDiv_plot %>% pull(beta_est) %>% .[c(1,6)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(1,6)], 
      cbind("E5-C8" = c(1,1), "Control-G7" = c(0,1)))$table
#           Estimates Standard Errors p-values
# E5-C8      0.3333768       0.1126736    0.003
# Control-G7 0.1489467       0.1352453    0.271

# H0: BC(G7-C8) = BC(Control-G7))
betta(betaDiv_plot %>% pull(beta_est) %>% .[c(2,3)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(2,3)], 
      cbind("G7-C8" = c(1,1), "Control-C8" = c(0,1)))$table
# Estimates Standard Errors p-values
# G7-C8       0.50400547      0.05657163    0.000
# Control-C8 -0.06656734      0.12878913    0.605

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(2,4)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(2,4)], 
      cbind("G7-C8" = c(1,1), "G7-E5" = c(0,1)))$table
# Estimates Standard Errors p-values
# G7-C8  0.50400547      0.07038573    0.000
# G7-E5 -0.09560367      0.15140975    0.528

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(2,5)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(2,5)], 
      cbind("G7-C8" = c(1,1), "Control-E5" = c(0,1)))$table
# Estimates Standard Errors p-values
# G7-C8       0.50400547      0.05256594    0.000
# Control-E5 -0.05762939      0.12066133    0.633

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(2,6)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(2,6)], 
      cbind("G7-C8" = c(1,1), "Control-G7" = c(0,1)))$table
# Estimates Standard Errors p-values
# G7-C8       0.50400547      0.03958053    0.000
# Control-G7 -0.02168198      0.08621946    0.801

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(3,4)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(3,4)], 
      cbind("Control-C8" = c(1,1), "G7-E5" = c(0,1)))$table
# Estimates Standard Errors p-values
# Control-C8  0.43743814      0.09096186    0.000
# G7-E5      -0.02903634      0.13702714    0.832

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(3,5)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(3,5)], 
      cbind("Control-C8" = c(1,1), "Control-E5" = c(0,1)))$table
# Estimates Standard Errors p-values
# Control-C8 0.437438137       0.0825685    0.000
# Control-E5 0.008937949       0.1137476    0.937

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(3,6)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(3,6)], 
      cbind("Control-C8" = c(1,1), "Control-G7" = c(0,1)))$table
# Estimates Standard Errors p-values
# Control-C8 0.43743814      0.07314933     0.00
# Control-G7 0.04488535      0.09058746     0.62

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(4,5)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(4,5)], 
      cbind("G7-E5" = c(1,1), "Control-E5" = c(0,1)))$table
# Estimates Standard Errors p-values
# G7-E5      0.40840180      0.08914168    0.000
# Control-E5 0.03797429      0.11670307    0.745

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(4,6)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(4,6)], 
      cbind("G7-E5" = c(1,1), "Control-G7" = c(0,1)))$table
# Estimates Standard Errors p-values
# G7-E5      0.40840180      0.08216683    0.000
# Control-G7 0.07392169      0.09965415    0.458

betta(betaDiv_plot %>% pull(beta_est) %>% .[c(5,6)], 
      betaDiv_plot %>% pull(beta_var) %>% sqrt %>% .[c(5,6)], 
      cbind("Control-E5" = c(1,1), "Control-G7" = c(0,1)))$table
# Estimates Standard Errors p-values
# Control-E5 0.4463761      0.07048166    0.000
# Control-G7 0.0359474      0.08857115    0.685

#beta diverstiy plot comparison figure

#add betaDiv between treatments for comparison
betaDiv_treatment$SE<-sqrt(betaDiv_treatment$beta_var)/sqrt(7)
  
#add SE column
betaDiv_plot$SE<-sqrt(betaDiv_plot$beta_var)/sqrt(2)

#now plot
betaDiv_plotFig<-ggplot(data=head(betaDiv_plot,6), aes(x=Covs, y=beta_est))+
  geom_hline(yintercept = betaDiv_treatment$beta_est + betaDiv_treatment$SE, color="gray80",alpha=0.7)+
  geom_hline(yintercept = betaDiv_treatment$beta_est, color="gray80",alpha=0.7, linetype=2)+
  geom_hline(yintercept = betaDiv_treatment$beta_est - betaDiv_treatment$SE, color="gray80",alpha=0.7)+
  geom_errorbar(aes(ymin=beta_est-SE, ymax=beta_est+SE,color=Location_Type_Compare),width=0, size=1.3)+
  geom_point(aes(color=Location_Type_Compare),size=3)+
  scale_color_manual(values=myColors)+
  theme(
    panel.background=element_rect(fill="white",color="black"),
    panel.border=element_rect(fill="transparent",color="black"),
    legend.position="none"
    # legend.background = element_rect(fill="white",color="black"),
    # legend.key=element_blank()
  )+
  labs(x="Plot", y="Bray-Curtis Dissimilarity Estimate")+
  guides(fill=guide_legend(title="Pair-wise Plot Comparison"), color="none")+
  ggtitle("C) Beta Diversity of Paired Plots")
  # annotate("text", x = 1, y=4.26, label = "a")+
  # annotate("text",x=2, y=4.02, label="ab")+
  # annotate("text",x=3, y=3.95, label="bc")+
  # annotate("text",x=4, y=3.42, label="d")
betaDiv_plotFig

####################################################################################################
#combine figs 4A, 4B and 4C

library(ggpubr)

divPlotsCombine<-ggarrange(alphaDiv_treatmentFig, alphaDiv_plotFig, betaDiv_plotFig,common.legend = FALSE,widths = c(0.25, 0.5),ncol=1)
divPlotsCombine

#save distance matrix
ggsave(divPlotsCombine, file="/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Results/Figures/Figure_4_Diversity_plots_7-9-2020.png", width=6, height=10, dpi=600)

####################################################################################################################################################
#Figure 5. Mean and CV of ASV counts by genera

#how many unique ASVs?
ASV.all.noZeros<-subset(ASV.all.keep.sub, Count>30)

ASV.all.noZeros.2<-subset(ASV.all.noZeros, c(Genus !="NA" & Genus !="Na"))

#test for differences among Genera
mod.4<-kruskal.test(Count~Genus, data=ASV.all.noZeros.2)
mod.4

#	Kruskal-Wallis rank sum test (if Count > 30)

# data:  Count by Genus
# Kruskal-Wallis chi-squared = 270.58, df = 149, p-value = 4.37e-09

# Kruskal-Wallis rank sum test (if Count > 5)
# 
# data:  Count by Genus
# Kruskal-Wallis chi-squared = 791.5, df = 271, p-value < 2.2e-16

#make Genus_Location_type colum
ASV.all.noZeros.2$Genus_Location_type<-paste(ASV.all.noZeros.2$Genus, ASV.all.noZeros.2$Location_type,sep="+")

#step 1, sum across replicates
ASV.all.sum.replicates<-aggregate(Count~Domain+Phylum+Class+Order+Family+Genus+Plot+Location+Location_type+Sample_type+Genus_Location_type, data=ASV.all.noZeros.2, FUN="sum")

#get mean and summary stats for ASV_count per family per Location_type
counts.per.genus.data<-summaryFunction(DataIn = ASV.all.sum.replicates, response="Count", factor="Genus_Location_type")

#make a plot of counts.per genus for suplementary materials

#merge with variables
ASV.genus.covs<-unique(ASV.all.noZeros[,c("Genus_Location_type","Genus","Location_type")])

counts.per.genus.merge<-merge(counts.per.genus.data, ASV.genus.covs, by="Genus_Location_type",all.x=TRUE)

#order the plot
counts.per.genus.order<-counts.per.genus.merge[order(counts.per.genus.merge$mean, decreasing=FALSE),]
counts.per.genus.order$Genus<-factor(counts.per.genus.order$Genus, levels=unique(as.character(counts.per.genus.order$Genus)))

length(unique(counts.per.genus.order$Genus))

#break groups into 3 groups: Just Forest, Both Forest and Control, and just Control

#get unique list of families in control only
control.data<-subset(counts.per.genus.order, Location_type=="Control")
controlListAll<-unique(control.data$Genus)

forest.data<-subset(counts.per.genus.order, Location_type=="Forest")
forestListAll<-unique(forest.data$Genus)

forestOnlyList<-setdiff(forestListAll, controlListAll)
controlOnlyList<-setdiff(controlListAll, forestListAll)

#create factors
counts.per.genus.order$Sample_Group<-ifelse(counts.per.genus.order$Genus %in% forestOnlyList, "Forest",
                                            ifelse(counts.per.genus.order$Genus %in% controlOnlyList, "Control", "Forest & Control"))

#set up factor levels
counts.per.genus.order$Sample_Group<-factor(counts.per.genus.order$Sample_Group, levels=rev(c("Control","Forest & Control", "Forest")))

counts.per.genus.order$Location_type<-factor(counts.per.genus.order$Location_type, levels=c("Control","Forest"))

levels(counts.per.genus.order$Location_type)

#subset data with only Forest & Control
forest.and.control.data<-subset(counts.per.genus.order, Sample_Group=="Forest & Control")

#get sum of mean by family
SumMean.genus<-aggregate(mean~Genus, FUN="sum",data=counts.per.genus.order)
colnames(SumMean.genus)<-c("Genus","SumMean")

counts.per.genus.order.merge<-merge(counts.per.genus.order, SumMean.genus, by="Genus",all.x=TRUE)

#order by SumMean
counts.per.genus.order.merge<-counts.per.genus.order.merge[order(counts.per.genus.order.merge$SumMean,decreasing=FALSE),]
counts.per.genus.order.merge$Genus<-factor(counts.per.genus.order.merge$Genus, levels=unique(counts.per.genus.order.merge$Genus))

#how many unique genera (total= 300)
length(unique(counts.per.genus.order.merge$Genus))

#forest and control (N = 246)
length(unique(forest.and.control.data$Genus))

#forest only (N = 49)
length(forestOnlyList)

#control only (N = 5)
length(controlOnlyList)

mycolors=c("tomato","skyblue3")

#set data up to plot all together (stack mean and CV)
#get data for mean
plot.genus.data.mean<-counts.per.genus.order.merge[,c("Genus","mean","Location_type","Sample_Group")]
plot.genus.data.mean$Variable<-"Mean"
names(plot.genus.data.mean)[names(plot.genus.data.mean)=="mean"]<-"Value"

#get data for CV
plot.genus.data.cv<-counts.per.genus.order.merge[,c("Genus","CV","Location_type","Sample_Group")]
plot.genus.data.cv$Variable<-"CV"
names(plot.genus.data.cv)[names(plot.genus.data.cv)=="CV"]<-"Value"

#combine data
plot.genus.data.combine<-rbind(plot.genus.data.mean, plot.genus.data.cv)

#reorder factors for Variable
plot.genus.data.combine$Variable<-factor(plot.genus.data.combine$Variable,levels=c("Mean","CV"))

#change treatment to Rainwater and Throughfall
levels(plot.genus.data.combine$Location_type)<-c("Rainwater","Throughfall")

#plot
countsPerGenusPlot<-ggplot(data=plot.genus.data.combine,aes(x=Genus, y=Value))+
  geom_bar(stat="identity", aes(fill=Location_type),position=position_stack(reverse=TRUE), alpha=0.95)+
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors)+
  coord_flip()+
  theme(
    panel.background = element_rect(fill="white",color="black"),
    panel.border= element_rect(fill="transparent",color="black"),
    legend.position=c(0.8,0.25),
    legend.key.size=unit(1.75,"line"),
    legend.background=element_rect(fill="transparent",color="black"),
    axis.title.x=element_blank()
  )+
  labs(x="Genus",y="CV count of reads")+
  guides(fill=guide_legend(title="Treatment"))
countsPerGenusPlot

#now facet plot
countsPerGenusPlotFacet<-countsPerGenusPlot+facet_grid(~Variable, scales="free",drop=TRUE)+theme(
  strip.text.x=element_text(size=14),
  strip.background.x=element_rect(fill="white"),
  strip.text.y=element_text(size=12)
)
countsPerGenusPlotFacet

#save plot
ggsave(countsPerGenusPlotFacet, file="/Users/zach/Dropbox (ZachTeam)/Manuscripts/eDNA_Forests_2020/Scientific_Reports/Figures/Figure_5_count_of_genera_mean_CV_11-30-2020.png", width=10, height=14, dpi=600)

################################################################################################
#Figure 6. Phylogenetic tree figure from metacoder package

#Set up data frame with the following columns: ASV, Locations: Control_1, . . ., lineageGenus (capitalized)
#now remove DI_H2O from ASV.all.keep
#ASV.all.keep.sub<-subset(ASV.all.keep, Sample_type !="DI_H2O")

#get numbers of unique ASVs

#how many unique ASVs?
ASV.all.noZeros.tree<-subset(ASV.all.keep.sub, Count>11)
length(unique(ASV.all.noZeros.tree$Genus))

#remove NA and Na
ASV.all.keep.tree<-subset(ASV.all.noZeros.tree, c(Genus !="NA" & Genus !="Na"))

length(unique(ASV.all.keep.tree$Genus))#150

#make Genus_Location_type colum
ASV.all.keep.tree$Genus_Location_type<-paste(ASV.all.keep.tree$Genus, ASV.all.keep.tree$Location_type,sep="+")

#step 1, sum across replicates
ASV.all.sum.replicates.tree<-aggregate(Count~Plot+Location+Location_type+ASV, data=ASV.all.keep.tree, FUN=sum)

#cast data
ASV.genera.1<-reshape2::dcast(data=ASV.all.sum.replicates.tree, formula=ASV~Location, value.var="Count", fill=0, fun.aggregate = sum)

names(ASV.all.sum.replicates.tree)
#subset to remove rows with all zeros
ASV.genera.1$rowSums<-rowSums(ASV.genera.1[,-1])
ASV.genera.2<-subset(ASV.genera.1, rowSums > 0)
ASV.genera.2$rowSums<-NULL

#convert any reads <= 5 to 0
#ASV.genera.2[ASV.genera.2 <=5] <- 0

#now bring in taxonomic info
taxa_predictions_1<-read.table(file="/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/zach_16s_pipeline-master/pipeline_output/sintax_taxonomy_predictions.txt",sep="\t",fill=TRUE)

#simplify and rename columns
taxa_predictions_2<-taxa_predictions_1[,c("V1","V4")]
colnames(taxa_predictions_2)<-c("ASV","lineage")

#merge with ASV.genera.2
ASV.genera.3<-merge(ASV.genera.2, taxa_predictions_2, by="ASV",all.x=TRUE)

#remove any rows where taxonomic assignment is missing
ASV.genera.4<-subset(ASV.genera.3, lineage!="")

#remove all species info in lineage strings for Family
ASV.genera.4$lineageFamily<-gsub(",g:.*","",gsub(",s:.*","",ASV.genera.4$lineage))

#remove species to create column with lineage up to Genus
ASV.genera.4$lineageGenus<-gsub(",s:.*","",ASV.genera.4$lineage)

#function to capitilize first character in a string
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  c.text<-read.table(text=y, sep=":")
  c.text.out<-paste(c.text$V1, paste(toupper(substring(c.text$V2, 1,1)), substring(c.text$V2,2), sep=""), sep=":")
  return(c.text.out)
}

#Capitalize taxonomic ranks
taxaNames<-ASV.genera.4$lineageGenus
taxaCapitalize<-list()
for(i in 1:length(taxaNames)){
  print(i)
  new.taxa<-taxaNames[i]
  
  taxa.names<-read.table(text=as.character(new.taxa), sep=",",fill=TRUE)
  
  taxa.caps<-apply(taxa.names, FUN=CapStr, MARGIN=2)
  
  taxa.caps.df<-data.frame(as.character(unlist(paste(taxa.caps, collapse=","))))
  
  colnames(taxa.caps.df)<-"Genus_Caps"
  
  taxaCapitalize<-rbind(taxaCapitalize, taxa.caps.df)
}

#remove any NAs

#now make new column with Capitalized Names
ASV.genera.4$Genus_Caps<-taxaCapitalize$Genus_Caps


#create object and parse out taxonomic info
obj <- parse_tax_data(ASV.genera.4,
                      class_cols = "Genus_Caps", # the column that contains taxonomic information
                      class_sep = ",", # The character used to separate taxa in the classification
                      class_regex = "^(.+):(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name"))

print(obj)

my.sample.info<-subset(sample.info.out, c(Location!="C8_2" & Sample_type !="DI_H2O"))
#set up factor for Location
my.sample.info$Location<-factor(my.sample.info$Location, levels=unique(my.sample.info$Location))
#set up factor for Location_type
my.sample.info$Location_type<-factor(my.sample.info$Location_type, levels=unique(my.sample.info$Location_type))

#set levels to Rainwater and Throughfall
levels(my.sample.info$Location_type)<-c("Rainwater","Throughfall")

#check that all ASVs have reads
no_reads <- rowSums(obj$data$tax_data[, as.character(my.sample.info$Location)]) == 0

sum(no_reads) #0

#check for uneven sampling
obj$data$tax_data <- calc_obs_props(obj, "tax_data")

print(obj)

#get per taxon abundance
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",cols = levels(my.sample.info$Location))

#now get Forest vs. Rain number of taxa
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = levels(my.sample.info$Location_type))

#comparing two groups
obj$data$diff_table <- compare_groups(obj,dataset = "tax_abund",
                                      cols = my.sample.info$Location, # What columns of sample data to use
                                      groups = my.sample.info$Location_type) # What category each sample is assigned to
print(obj$data$diff_table)

range(obj$data$diff_table$log2_median_ratio, na.rm=TRUE, finite=TRUE)

set.seed(999)
heatTreePlot<-heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_interval = c(-1.5,7), # The range of `log2_median_ratio` to display
          node_color_range = c("skyblue3", "gray", "tomato"), # The color palette used
          node_size_axis_label = "ASV count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          margin_size=c(0.1,0.1,0,0))
#save with ggsave
ggsave(heatTreePlot, file="./Results/Figures/Figure_6_Compare_Forest_and_Rain_Genera_heat_tree_6-15-2020.png", width=9, height=8, dpi=600)

##################################################################################################################3
#End script