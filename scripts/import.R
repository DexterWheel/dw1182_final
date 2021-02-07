library(dplyr)

#Old dataset from workshop
load(url("https://www-users.york.ac.uk/~dj757/BIO00047I/data/yeast_data.28-02-2020.Rda"))

#data from database
angeli <- read.delim("data-old/AnGeLiDatabase.txt",h=T)

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

#next I need add the other columns to this dataframe

#I use the colnames function to identify the column number of the desired columns
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


