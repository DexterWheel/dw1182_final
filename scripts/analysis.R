######
#Section 1: Are Essential Genes Expressed Differently to Non-Essential Genes?
######

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

load(file = "data-processed/gene_f.Rda")

nrow(gene_f)
class(gene_f$essential)
class(gene_f$gene.expression.RPKM)


#only protein coding genes
prot <- gene_f %>%
  subset(protein_coding ==1) %>%
  subset(!is.na(essential)) %>%
  subset(!is.na(gene.expression.RPKM))

#around 80000 rows removed
nrow(prot) 

#need to remove background on all of these
#need to change fonts and positions
#need to change 1 and 0 to essential and non-essential
#need to make functions with presets for all three boxplot types
#combine all the plots into one, could potentially make the function apply all three at once. High level... you could even code a choice for which of the three boxplot types you want to use

#I use a wilcoxon signed rank test to determine significance

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

ess <- subset(prot, essential==1)
non.ess <- subset(prot, essential==0)

wilcox.test(ess$gene.expression.RPKM, non.ess$gene.expression.RPKM)
#the difference observed is very unlikely to occur by chance (Null hypothesis rejected) (W = 3039370, p-value < 2.2e-16)
#essential genes are more highly expressed than non-essential genes in general
######
#Section 2: The relationship between mRNA Stability and RPKM in Essential and Non-Essential mRNAs
######

#Plotting the relationship between mRNA stability and RPKM
plot_sctr(prot,
         log10(prot$mRNA.stabilities),
         log10(prot$gene.expression.RPKM),
         "Boxplot",
         "Essentiality",
         "mRNA.stabilities")

#logs of both datasets were taken to visualise the trend easier 

#making a subset removing missing essentialityy and mrna stability data
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

#I use another wilcoxon signed rank test to determine significance
wilcox.test(prot2$mRNA.stabilities, prot2$gene.expression.RPKM)
#Null hypothesis rejected (W = 7856176, p-value < 2.2e-16)
#More stable mRNAs tend to be more highly expressed


#Third wilcoxon signed rank test to determine significance
wilcox.test(ess$mRNA.stabilities, non.ess$mRNA.stabilities)
#Null hypothesis rejected (W = 2372773, p-value = 3.728e-08)
#mRNA of essential genes is more stable on average

######
#Section 3: The Relationship Between mRNA Stability and Localisation in Essential and Non-Essential mRNAs
######

nams<- names(gene_f[,1:8])

locale <- gene_f %>% subset(!is.na(mRNA.stabilities)) %>% subset(!is.na(essential)) %>% pivot_longer(names_to = "localisation",
                                                                       values_to = "present",
                                                                       cols = -nams  ) %>% subset(present == 1)

#I create summary data for the plot
summary <- summarySE(locale, measurevar = "mRNA.stabilities", groupvars = c("localisation", "essential"))
summary

mod <- aov(data = locale, mRNA.stabilities ~ localisation * essential)
summary(mod)

#Significant main effect of localisation on mRNA stability (ANOVA: F = 16.860; d.f. = 38; p< 2e-16)

#Significant main effect of essentiality on mRNA stability (ANOVA: F = 9.218; d.f. = 1; p= 0.00241)

#There was an interaction  (ANOVA: F = 2.691; d.f. = 35; p= 2.97e-07)

#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#####
TukeyHSD(mod)

#there are more than 5000 values so a shapiro test cannot be done with all of them
res <- mod$residuals

#Testing the normal distribution assumption
shapiro.test(res[1:5000])
#Null hypothesis rejected (W = 0.81428, p-value < 2.2e-16)
#this assumption is not met as the data is not normally distributed

#Testing equal variances assumption
#for the leveneTest function which can test for heterogeneity of variance
leveneTest(mRNA.stabilities ~ localisation * essential, data = locale_mrna)
#Null hypothesis rejected (p-value < 2.2e-16)
#equal variances assumption not met

#####
#names for barplot axis
#####
nams2<- names(gene_f[,9:47])

nams2

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
nams2[20]="noGTc"
nams2[21]="noGTGTn"
nams2[22]="noGTnGTc"
nams2[23]="noGTnGTGTc"
nams2[24]="nGTn"
nams2[25]="noGTnQQc"
nams2[26]="nol"
nams2[27]="nGTc"
nams2[28]="nGTQQc"
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


