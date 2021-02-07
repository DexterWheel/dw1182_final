 
#######################################################################################
#                                                                                     #
#                               IMPORTING THE DATA                                    #  
#                                                                                     #
#######################################################################################

setwd("~/Biology/R/Biology/Big Data Biology/Data")

#First I load in the two data sets required for my analysis

load(url("https://www-users.york.ac.uk/~dj757/BIO00047I/data/yeast_data.28-02-2020.Rda"))
angeli <- read.delim("data-raw/AnGeLiDatabase.txt",h=T)

#I then check to see if the data has been imported correctly
ls()
class(gene)
nrow(gene)
ncol(gene)
head(gene)

class(angeli)
nrow(angeli)
ncol(angeli)
head(angeli)

#Here I insert cytosol an ER localisation data from the angeli database into the gene dataframe

#First I remove gene ontology columns
#I get a the vector of the column names
nam <- names(angeli)

#Using grep I can locate gene ontology or "GO." columns
#I invert the grep, so it does not locate these columns
not.go.columns <- grep("GO.",nam,invert=T)

#using this I create a subset of angeli without GO columns 

angeli2 <- angeli[,not.go.columns]

#I also need to remove 'FYPO', (fission yeast phenotype ontology) or PFAM groups
#so I remove them with grep also
not.fypo.or.pfam.columns <- grep("FYPO|PF",names(angeli2),invert=T)
angeli3 <- angeli2[,not.fypo.or.pfam.columns]

#I create a subset with just the information rows so I can locate the data I need
info <- angeli3[1:7,]
View(info)
#I also create a subset with only the rows with data
gene.data <- angeli3[8:7012,]

#renaming first column
names(gene.data)[1]="gene"

#I create subset of the data containing just the two columns I need
my.data.col <- which(names(gene.data) == "Cytosol")
my.data.col2 <- which(names(gene.data) == "ER")
angeli4 <- gene.data[,c(1,my.data.col, my.data.col2)]
View(angeli4)

#I combine angeli4 with the gene dataframes, and save it as an Rda file to save 
#time when reopening the script
gene_final<-merge(gene,angeli4,by="gene",all=T)
save(gene_final,file="gene_finaldata.Rda")

#Finally I check that the dataframe has been generated correctly
class(gene_final) #data.frame
nrow(gene_final)#how many rows?
#7009
ncol(gene_final)#columns?
#25
head(gene_final)#first few rows

#I can use this file for future use without having to import the whole database again
load(file="gene_finaldata.Rda")

#######################################################################################
#                                                                                     #
#                                 DATA ANALYSIS                                       #  
#                                                                                     #
#######################################################################################

######
######Section 1: Are Essential Genes Expressed Differently to Non-Essential Genes?
######

#first I make a new data frame with just protein coding data

prot_pre <- subset(gene_final, protein_coding ==1)
#Checking rows have been removed
nrow(prot_pre)

#I make the essential column into a factor (it's currently numeric) This will change the 
#data from continuous to discrete categorical which is what it should be.
prot_pre$essential <- as.factor(prot_pre$essential)

#I use this to make sure any rows that are missing data are removed
prot <- subset(prot_pre, !is.na(essential))
nrow(prot) #fortunately none were removed at this stage

#I now split the protein coding data into essential and non essential proteins
ess <- subset(prot, essential==1)
non.ess <- subset(prot, essential==0)

#Plotting the difference in expression of essential and non-essential genes
boxplot(
  log10(ess$gene.expression.RPKM),
  log10(non.ess$gene.expression.RPKM),
  cex.axis=1.05,
  cex.lab=1.5,
  ylab="Gene Expression log10(RPKM)",
  xlab="Essentiality",
  names=c("Esssential","non-Essential"),
  outline=F)
#log 10 in order to visualise data more clearly
#Outliers removed for a clearer plot

#I use a wilcoxon signed rank test to determine significance
wilcox.test(ess$gene.expression.RPKM, non.ess$gene.expression.RPKM)
#the difference observed is very unlikely to occur by chance (Null hypothesis rejected) (W = 3039370, p-value < 2.2e-16)
#essential genes are more highly expressed than non-essential genes in general

######
######Section 2: The relationship between mRNA Stability and RPKM in Essential and Non-Essential mRNAs
######

#Plotting the relationship between mRNA stability and RPKM

plot(log10(prot$mRNA.stabilities), log10(prot$gene.expression.RPKM), 
     xlab = "log10(mRNA Stability)",
     ylab= "Gene Expression Log10(RPKM)",
     cex.axis=1.05,
     cex.lab=1.5)
#logs of both datasets were taken to visualise the trend easier 

#I use another wilcoxon signed rank test to determine significance
wilcox.test(prot$mRNA.stabilities, prot$gene.expression.RPKM)
#Null hypothesis rejected (W = 9156457, p-value < 2.2e-16)
#More stable mRNAs tend to be more highly expressed

#Plotting the relationship between esentiality and mRNA stability
boxplot(ess$mRNA.stabilities, non.ess$mRNA.stabilities,
        outline=F,
        ylab="mRNA stability", 
        cex.axis=1.05,
        xlab="Essentiality",
        cex.lab=1.6,
        names=c("Esssential","non-Essential"))
#Outliers removed for a clearer plot

#Third wilcoxon signed rank test to determine significance
wilcox.test(ess$mRNA.stabilities, non.ess$mRNA.stabilities)
#Null hypothesis rejected (W = 2454908, p-value = 3.743e-08)
#mRNA of essential genes is more stable on average 

######
######Section 3: The Relationship Between mRNA Stability and Localisation in Essential and Non-Essential mRNAs
######

#First I need to make localisation a group variable

#I create subsets of each localisation 

Prot_dotssubset <- subset(prot, Nuclear_dots ==1)
Prot_envsubset <- subset(prot, Nuclear_envelope ==1) 
Prot_nolussubset <- subset(prot, Nucleolus ==1)
Prot_nucsubset <- subset(prot, Nucleus ==1)
Prot_cysubset <- subset(prot, Cytosol ==1)
Prot_ersubset <- subset(prot, ER ==1)

#And change the data from binary (0 or 1) to a category e.g nucleolus ("nolus") 

Prot_dotssubset$localisation="dots" #nuclear dots
Prot_envsubset$localisation="n.en" #nucler envelope
Prot_nolussubset$localisation="nolus" #nucleolus
Prot_nucsubset$localisation="nuc" #nucleus
Prot_cysubset$localisation="cyt" #cytoplasm
Prot_ersubset$localisation="ER" #ER


#I re-merge the subsets, but now all localisations are in one column
library(dplyr)#for the bind_rows function


locale <- bind_rows(Prot_dotssubset,Prot_envsubset,
                    Prot_nolussubset,Prot_nucsubset,
                    Prot_cysubset,Prot_ersubset,)

##################################################################################

#I remove N/a values so that the summary function works
locale_mrna<- subset(locale, !is.na(mRNA.stabilities))

library(Rmisc)#for the summarySE function

#I create summary data for the plot
summary <- summarySE(locale_mrna, measurevar = "mRNA.stabilities", groupvars = c("localisation", "essential"))
summary

#Two-way ANOVA is carried out with mRNA stability being an explanatory variable and 
#localisation and essentiality as response variables
mod <- aov(data = locale_mrna, mRNA.stabilities ~ localisation * essential)
summary(mod)
#Significant main effect of localisation on mRNA stability (ANOVA: F = 43.510; d.f. = 5; p< 2e-16)
#No significant main effect of essentiality on mRNA stability (ANOVA: F = 1.357; d.f. = 1; p= 0.24422)
#There was an interaction  (ANOVA: F = 3.709; d.f. = 5; p= 0.00241)
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#TukeyHSD provides detail on the differences of each localisations
TukeyHSD(mod)
#I chose not to refer to them in the text as the intention was just to determine 
#whether there was a main effect and an interaction
#Interesting that the Nucleoulus pattern is completely different from the rest however

#Testing the normal distribution assumption
shapiro.test(mod$residuals)
#Null hypothesis rejected (W = 0.81428, p-value < 2.2e-16)
#this assumption is not met as the data is not normally distributed

#Testing equal variances assumption
library(car)#for the leveneTest function which can test for heterogeneity of variance
leveneTest(mRNA.stabilities ~ localisation * essential, data = locale_mrna)
#Null hypothesis rejected (p-value < 2.2e-16)
#equal variances assumption not met

library(ggplot2)#For ggplot function

#barplot of mRNA stability with respect to localisation an essentiality
ggplot(data = summary, aes(x=localisation, y=mRNA.stabilities, fill=essential,)) +
  geom_bar(stat= "identity",position=position_dodge())+
  geom_errorbar(aes( ymin = mRNA.stabilities - se, 
                     ymax = mRNA.stabilities + se), 
                width = 0.2, 
                position = position_dodge(0.8), 
                colour = c("Black"), size = 0.5) +
   ylab("mRNA Stability") +
  scale_x_discrete(breaks=c("dots", "n.en", "nolus", "nuc", "cyt", "ER"),
                   labels=c("Nuclear Dots", "Nuclear Envelope", "Nucleolus", "Nucleus", "Cytosol", "ER"))+
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 20, r = 40, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 1, r = 20, b = 0, l = 0)),
        legend.title = element_text(size= 20),
        legend.text = element_text(size = 13), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        legend.key = element_blank())+
  xlab("Localisation")