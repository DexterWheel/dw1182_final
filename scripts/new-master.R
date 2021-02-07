 
#######################################################################################
#                                                                                     #
#                               IMPORTING THE DATA                                    #  
#                                                                                     #
#######################################################################################

library(ggplot2)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(patchwork)
library(presetplots)
library(Rmisc)
library(car)

#Old dataset from workshop
load(url("https://www-users.york.ac.uk/~dj757/BIO00047I/data/yeast_data.28-02-2020.Rda"))

#data from database
angeli <- read.delim("data-raw/AnGeLiDatabase.txt",h=T)

ls()
class(gene)
nrow(gene)
ncol(gene)
head(gene)

class(angeli)
nrow(angeli)
ncol(angeli)
head(angeli)

#selecting the localisation data based on the data in row 3. 
#Group refers to first column which allows datasets to be merged

cols <- c("Group","Protein Localisations (ORFeome)")

angeli2 <- angeli %>%
  select(grep("FYPO|PF|GO.",colnames(angeli),invert=T)) %>%
  t() %>%
  data.frame() %>%
  filter(X3 %in% cols ) %>%
  t() %>%
  data.frame()

#removing info columns  
angeli2 <- angeli2[8:7012,]

#changing the first column to the same name as in the gene dataset
names(angeli2)[1] = "gene"

#add the other columns to this dataframe

#use the colnames function to identify the column number of the desired columns
colnames(gene)
gene <- gene[,c(1, 4, 7, 8, 9, 11, 22, 25)]

#combining the two datasets
gene_f <- merge(gene,angeli2, by = "gene", all = T)
save(gene_f, file = "data-processed/gene_f.Rda")

class(gene_f)
#should be data.frame
nrow(gene_f)#how many rows?
#should be 7009
ncol(gene_f)#columns?
#should be 47
head(gene_f)

#I can use this file for future use without having to import the whole database again
load(file = "data-processed/gene_f.Rda")


#######################################################################################
#                                                                                     #
#                                 DATA ANALYSIS                                       #  
#                                                                                     #
#######################################################################################


#####Section 1: Are Essential Genes Expressed Differently to Non-Essential Genes?####

load(file = "data-processed/gene_f.Rda")

nrow(gene_f)
class(gene_f$essential)
class(gene_f$gene.expression.RPKM)

#only protein coding genes, removing missing values
prot <- gene_f %>%
  subset(protein_coding ==1) %>%
  subset(!is.na(essential)) %>%
  subset(!is.na(gene.expression.RPKM))

nrow(prot)
#around 2000 rows removed

#vector for choosing the plot output, see plots.R in presetplots package provided
plot_choice <- c("box", "jit", "vio")

plots(plot_choice,
      prot, 
      prot$essential, 
      log10(prot$gene.expression.RPKM), 
      prot$essential,
      "Boxplot",
      "Essentiality",
      "Gene Expression log10(RPKM)",
      c("non-Essential","Esssential"))

# subset for wilcoxon
ess <- subset(prot, essential==1)
non.ess <- subset(prot, essential==0)

#wilcoxon signed rank test to determine significance
wilcox.test(ess$gene.expression.RPKM, non.ess$gene.expression.RPKM)
#the difference observed is very unlikely to occur by chance 

##(Null hypothesis rejected)##

###essential genes are more highly expressed than non-essential genes in general###

#####Section 2: The relationship between mRNA Stability and RPKM in Essential and Non-Essential mRNAs####

###scatterplot: Plotting the relationship between mRNA stability and RPKM
plot_sctr(prot,
          log10(prot$mRNA.stabilities),
          log10(prot$gene.expression.RPKM),
          "Boxplot",
          "gene.expression.RPKM",
          "mRNA.stabilities")

#logs of both datasets were taken for visualisation

#another wilcoxon signed rank test to determine significance
wilcox.test(prot$mRNA.stabilities, prot$gene.expression.RPKM)
#the difference observed is very unlikely to occur by chance 

##(Null hypothesis rejected)##

###More stable mRNAs tend to be more highly expressed###

#making a subset removing missing essentiality and mrna stability data
prot2 <- gene_f %>%
  subset(protein_coding ==1) %>%
  subset(!is.na(essential)) %>%
  subset(!is.na(mRNA.stabilities))

plots(plot_choice,
      prot2,
      prot2$essential,
      log10(prot2$mRNA.stabilities),
      prot2$essential,
      "Boxplot",
      "Essentiality",
      "mRNA.stabilities",
      c("non-Essential","Esssential") )

#Third wilcoxon signed rank test to determine significance,
wilcox.test(ess$mRNA.stabilities, non.ess$mRNA.stabilities)
#the difference observed is very unlikely to occur by chance

##Null hypothesis rejected##

##mRNA of essential genes is more stable on average##

#####Section 3: The Relationship Between mRNA Stability and Localisation in Essential and Non-Essential mRNAs####

nams<- names(gene_f[,1:8])

#making localisation data a single variable
locale <- gene_f %>%
  subset(!is.na(mRNA.stabilities)) %>% 
  subset(!is.na(essential)) %>% 
  pivot_longer(names_to = "localisation",
               values_to = "present",
               cols = -nams  ) %>% 
  subset(present == 1)

#I create summary data for the plot
summary <- summarySE(locale, measurevar = "mRNA.stabilities", groupvars = c("localisation", "essential"))
summary

mod <- aov(data = locale, mRNA.stabilities ~ localisation * essential)
summary(mod)

#Significant main effect of localisation on mRNA stability

#Significant main effect of essentiality on mRNA stability

#There was an interaction

TukeyHSD(mod)

#there are more than 5000 values so a shapiro test cannot be done with all of them
res <- mod$residuals

#Testing the normal distribution assumption
shapiro.test(res[1:5000])
#Null hypothesis rejected

##this assumption is not met as the data is not normally distributed##

#Testing equal variances assumption
#for the leveneTest function which can test for heterogeneity of variance
leveneTest(mRNA.stabilities ~ localisation * essential, data = locale_mrna)
#Null hypothesis rejected

##equal variances assumption not met##


#names for barplot axis#####
nams2<- names(gene_f[,9:47])

nams2

nams2<- names(gene_f[,9:47])
nams2[1]="Amb"
nams2[2]="cdc"
nams2[3]="cdcts"
nams2[4]="cdp"
nams2[5]="cds"
nams2[6]="cdot"
nams2[7]="cmt"
nams2[8]="cGTn"
nams2[9]="cyt"
nams2[10]="cyQQn"
nams2[11]="ER"
nams2[12]="Fil"
nams2[13]="gol"
nams2[14]="mitu"
nams2[15]="mito"
nams2[16]="n.a.s"
nams2[17]="n.d"
nams2[18]="ndot"
nams2[19]="nen"
nams2[20]="no
GT
c"
nams2[21]="no
GT
GT
n"
nams2[22]="no
GT
n
GT
c"
nams2[23]="no
GT
n
GT
GT
c"
nams2[24]="n
GT
n"
nams2[25]="no
GT
n
QQ
c"
nams2[26]="nol"
nams2[27]="n
GT
c"
nams2[28]="n
GTQQ
c"
nams2[29]="nuc"
nams2[30]="o.n.d"
nams2[31]="pcts"
nams2[32]="pct"
nams2[33]="ps"
nams2[34]="per"
nams2[35]="sep"
nams2[36]="SPB"
nams2[37]="spmi"
nams2[38]="vac"
nams2[39]="vac.m"

#####
#barplot of mRNA stability with respect to localisation an essentiality

plot_bar(summary,
         summary$localisation,
         summary$mRNA.stabilities,
         summary$essential, 
         "barplot",
         "Localisation",
         "mRNA Stability",
         nams2)
