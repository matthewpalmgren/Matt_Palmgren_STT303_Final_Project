library(tidyverse)
library(ggplot2)
install.packages('readxl')
library(readxl)
install.packages("zip")
library('zip')
library(dplyr)
library(data.table)
library(tibble)
install.packages('factoextra')
library(factoextra)
install.packages('pheatmap')
library(pheatmap)
install.packages('caret')
library(caret)
install.packages("randomForest")
library(randomForest)
install.packages('nnet')
library(nnet)
# 1 CLEANING AND TIDYING DATA

# Reading Data Sets into R

genefile<-"C:/Users/matth/OneDrive/Documents/GTEx_Analysis_gene_median_tpm.gct.txt"
gene_data<-read.delim(genefile, skip=2, header = T, stringsAsFactors = F)

hippo_file<-"C:/Users/matth/OneDrive/Documents/Brain hippocampus gene expression by subject.txt"
hippo_data<-read.delim(hippo_file,skip = 2,header = T, stringsAsFactors = F)

hippo_meta<-fread("https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")



subjIDS<-colnames(hippo_data)[-c(1:2)]
CLsubjIDS<-sub("\\..*","",subjIDS)
Vhippo_meta<-hippo_meta[hippo_meta$SUBJID %in% CLsubjIDS, ]

trpsd_hippo_data<-as.data.frame(t(hippo_data))
#colnames(trpsd_hippo_meta)<-Vhippo_meta$SUBJID
trpsd_hippo_data<-trpsd_hippo_data[-1,]
correct_hippo_data <- trpsd_hippo_data|>
  rownames_to_column(var = "rowname")
colnames(correct_hippo_data) <- as.character(correct_hippo_data[1,])
correct_hippo_data <- correct_hippo_data[-1,]

colnames(correct_hippo_data)[1] <- "SUBJID"
colnames(correct_hippo_data) <- make.names(colnames(correct_hippo_data), unique = TRUE)
correct_hippo_data$SUBJID <- substr(correct_hippo_data$SUBJID,1,10)
correct_hippo_data$SUBJID <- substr(correct_hippo_data$SUBJID,6,10)
hippo_meta$SUBJID <- substr(hippo_meta$SUBJID,6,10)

joined_data <- inner_join(hippo_meta,correct_hippo_data,by = "SUBJID")

Median_data<-t(gene_data[,-c(1,2)])
colnames(Median_data)<-gene_data$Description
Median_data1<-apply(Median_data, 2, function(x) as.numeric(x))
Median_data1<-apply(Median_data1, c(1,2), function(x) ifelse(x==0, NA, x))
rownames(Median_data1)<-rownames(Median_data)
#   02   EXPLORATORY ANLAYSIS

#Preliminary Data Structure Analysis

#PCA - Median Gene Expression Data 

CLmedGnData<-Median_data[-c(1,2),]
numCLmedData<-apply(CLmedGnData,2,as.numeric)
blankGenes<-apply(numCLmedData,2,function(x) var(x, na.rm = T)==0)
blankTab<-which(blankGenes)
refinedPCAmedData<-numCLmedData[, !blankGenes]
PCAmedData<-prcomp(refinedPCAmedData, scale.=F)
summary(PCAmedData)


print(PCAmedData$rotation)

ExplainedVar<-PCAmedData$sdev^2
ProportionVar<-ExplainedVar/sum(ExplainedVar)
ScreeDta<-data.frame(PC=seq_along(ProportionVar), Variance=ProportionVar)

ggplot(ScreeDta, aes(x=PC, y= Variance))+
  geom_point(color='seagreen', size=3)+
  geom_line(color='turquoise')+
  theme_minimal()+
  labs(title = 'Scree Plot', x='Principle Component', y='Proportion of Variance')


#Most of variance can be explained by first four PCs covering the 40% and the rest
#'comprising of 42%
#'PC1 could be like specific tissues with similar expression profiles or specific genes
#'contributing heavily to variance



rownames(refinedPCAmedData) <- paste0("Sample_", seq_len(nrow(refinedPCAmedData)))
biplot(PCAmedData)
#'Black group - Tissues - pretty clustered together and no real detectable variation
#'Red group - Genes - Few long arrows pointing up to the right from the main group
#'are the biggest thing to notice (Drivers of PC1)
#'
#'Overall very clustered to toward the center at 0.0 meaning not a great amount of
#'genes are influenced by PC1 or PC2


# EXploring the upper right facing genes
PC1Ldings<-PCAmedData$rotation[,1]
TopPC1Gns<-sort(PC1Ldings,decreasing = T)[1:10]



#' Which Tissues Express the Most Genes

MaxDtaExprsmed<-Median_data[-c(1,2),]
expressCNT<-apply(MaxDtaExprsmed,1,function(row) sum(row>1.0))
TissueMaxExpression<-names(which.max(expressCNT))
cat('Tissue with most expressed genes (>1.0 TPM)', TissueMaxExpression, '\n')

# Most expressed genes is the in the Testis
barplot(expressCNT,
        main = 'Expressed genes per tissue',
        xlab = 'Tissue Type',
        ylab = 'Number of Expressed genes (>1.0 TPM)',
        las=2,
        col='darkgreen')

# Which Tissue has the Least amount of Genes Expressed
TissueMinExpression<-names(which.min(expressCNT))
cat('Tissue with least expressed genes (>1.0 TPM)', TissueMinExpression, '\n')


# Tissue with the greatest range in gene expression values
ExprRanges<-apply(Median_data1, 1, function (row) max(row, na.rm=T) - min(row, na.rm=T))
GreatestRangeTissue<-names(which.max(ExprRanges))
cat('Tissue with the greatest range of gene expression values (non-zero values)', GreatestRangeTissue, '\n')

# Tissue with the smallest range in gene expression values
SmallestRangeTissue<-names(which.min(ExprRanges))
cat("Tissue with the smallest range of gene expression values (non-zero values)", SmallestRangeTissue, '\n')


# Visualizing Ranges
ExprRangeDF<-data.frame(
  Tissue=names(ExprRanges),
  Range= ExprRanges
)
ggplot(ExprRangeDF,
       aes(x=Tissue, y= Range))+
  geom_bar(stat = 'identity', fill='forestgreen')+theme_minimal()+
  labs(title='Range of gene Expression Per Tissue (TPM)',
       x='Tissues',
       y='Range')

# Which gene is expressed in the most tissues
XprsdGeneCount<-apply(Median_data1, 2, function(col) sum(!is.na(col)))
MostXprsd<-names(which.max(XprsdGeneCount))
cat('Most expressed gene across tissue samples:', MostXprsd, '\n')

# Which gene is least expressed across tissue types
LeastXprsd<-names(which.min(XprsdGeneCount))
cat('Least expressed gene across tissue types:', LeastXprsd, '\n')

#Visualizing gene expression across tissues
XprsXTissueDF<-data.frame(
  Gene=names(XprsdGeneCount),
  Tissue_Types=XprsdGeneCount
)
ggplot(XprsXTissueDF, aes(x=Gene, y=Tissue_Types))+
  geom_bar(stat='identity', fill='lightblue')+
  theme_minimal()+
  labs(title = 'Number of tissues in which the gene is expressed',
       x='Genes',
       y='Number of Tissues')
# Above is nearly impossible to read, way too many genes and the scale is completley off

# Comparing gene expression profiles of sex organs
SxTissDta<-Median_data1[c('Testis','Ovary'),]
SxTissDF<-data.frame(
  Gene=colnames(SxTissDta),
  Testis=SxTissDta[1,],
  Ovary=SxTissDta[2,]
)

# Visualization
ggplot(SxTissDF, aes(x=Testis, y=Ovary))+geom_point(alpha=0.6, color='blue')+
  theme_minimal()+
  labs(title = 'Gene Expression Comparison of Males (Testis) and Female(Ovary) Sex Organs',
       x='Expression in Testis (TPM)',
       y='Expression in Ovary (TPM)')+
  geom_abline(slope = 1, intercept = 0, linetype='dashed', color='red')
#' points close to the line indicate similar levels of expression
#' The few points that are very close to the x-axis and far out from center line
#' idicate these are genes that are highly expressed in Ovaries and hardly expressed in Testis
#' Dense cluster of points going up the Y-axis means there are a lot of genes that are more
#' expressed in the testis than there are in the ovaries
#' Very FEW genes are similar in levels of expression meaning its more likley that they
#' dissimilar


# Calculating Correlation Coefficient

SxTissCor<-cor(SxTissDta[1,], SxTissDta[2,], use = 'complete.obs')
cat('Correlation between Testis gene expression and Ovary gene expression:', SxTissCor, '\n')
#' closer to one means more correlated however if you look at the plot
#' Most of the data points are clustered very close in the left corner of the graph
#' meaning low expression values for both tissues
#' 


# Determining which tissues have the most similar and least similar gene expression
# profiles

# Whole Data Set Correlation Matrix
#*Since Data set is very large I will just analyze top 100 most variable genes


VarGenes<-apply(Median_data1, 2, var, na.rm=T)
Var100<-order(VarGenes,decreasing = T)[1:100]

Var100Med<-Median_data1[,Var100]
Tvar100Med<-t(Var100Med)
Var100Mtrx<-cor(Tvar100Med, use='pairwise.complete.obs')

# Visualizing Correlation in Gene expression profiles
pheatmap(Var100Mtrx,
         color = colorRampPalette(c('blue', 'white', 'red'))(100),
         cluster_rows = T,
         cluster_cols = T,
         main = 'Tissue to Tissue Correlation Heatmap'
         
         )

# It looks like a vast majority of the genes are highly correlated except for the last
#' which show very little correlation to each other
#' Possibly two distinct groups then???
#' 

# Which two tissues are most correlated and least correlated

MostCorTiss<-which(Var100Mtrx==max(Var100Mtrx[lower.tri(Var100Mtrx)]), arr.ind = T)
cat('Most similar tissues: ', rownames(Var100Mtrx)[MostCorTiss[,1]],'\n')

LeastCorTiss<-which(Var100Mtrx==min(Var100Mtrx[lower.tri(Var100Mtrx)]), arr.ind = T)
cat('Least similar tissues: ', rownames(Var100Mtrx)[LeastCorTiss[,1]],'\n')








# Determining if Gene Expression levels can predict age, sex etc.

# Data Prep

SSAD<-joined_data[,c('SUBJID', 'SEX', 'AGE', 'DTHHRDY')]
JeenXprsn<-joined_data[,-c(1:4)]

CLHippoDta<-na.omit(joined_data)


#Normalizing gene expression values??

JeenXprsnNum<-apply(JeenXprsn, 2, function(x) as.numeric(x))
JeenXprsn<-scale(JeenXprsnNum)


#Visualizing Data

table(SSAD$SEX)
table(SSAD$AGE)
table(SSAD$DTHHRDY)

boxplot(JeenXprsn[,'WASH7P'] ~ SSAD$SEX,
        main='WASH7P by Sex',
        xlab = 'SEX',
        ylab = 'Expression Level'
        )
boxplot(JeenXprsn[,'WASH7P'] ~ SSAD$AGE,
        main='WASH7P by Age',
        xlab = 'Age',
        ylab = 'Expression Level')


# Statistical Tests for WASH7P

WASHttestSex<-t.test(JeenXprsn[,'WASH7P']~SSAD$SEX)



WASHAgeReg<-lm(JeenXprsn[,'WASH7P']~SSAD$AGE, data = SSAD)
summary(WASHAgeReg)
#' Essentially, moderate significance found in 60-69, and 70-79 age groups, 
#' howeverAdjusted and multiple R squared values indicate the model is a poor fit and
#' large p-value shows model is not statistically significant


# Finding Heat Shock Proteins
HSPindex<-which(startsWith(colnames(joined_data), 'HSP'))

boxplot(JeenXprsn[,'HSPE1P8'] ~ SSAD$SEX,
        main='HSPE1P8 by Sex',
        xlab = 'SEX',
        ylab = 'Expression Level'
)
#Linear Regression for HSPE1P8

HSPE1P8RegTest<-lm(JeenXprsn[,'HSPE1P8']~SSAD$SEX, data = SSAD)
summary(HSPE1P8RegTest)
# Doesn't appear to be a great model



# Predicting based off meta data

# Data Prep
Vari<-apply(JeenXprsn,2,var, na.rm=T)
Top1000VarGn<-order(Vari, decreasing = T)[1:1000]
TopJeenExprsn<-JeenXprsn[,Top1000VarGn]

# PCA with Sex
HippoPCA<-prcomp(t(TopJeenExprsn), scale. = T)
plot(HippoPCA$x, col=SSAD$SEX, pch=19, main='PCA of Hippocampus Gene Expression by Sex')
VarExplained<-HippoPCA$sdev^2/sum(HippoPCA$sdev^2)
#PC1=12% total variance
#'PC2=5.5% total Variance
#'PC3=2.6% total Variance
#'PC4=2% total Variance
#'
#'Essentially, data is mostly explained by first 2 PCAs and many samples have similar
#'profiles along the PC1 axis 
#'EXpression data has similar shared variance across across all samples
#' Sex does not have strong influence on gene expression across genes in hippocampus
#' tissue 
#' 
# PCA with AGE
SSAD$AGE<-as.factor(SSAD$AGE)
plot(HippoPCA$x, col=as.numeric(SSAD$AGE), pch=19, main = 'PCA of Hippocampus Gene EXpression by Age')
#' Basically the same story as with Sex in that Age has no real influence on 
#' Gene expression in hippocampus tissue
#' 
#' 

# Kmeans Clustering
TopJeenExprsnScale<-scale(TopJeenExprsn)

#Kmeans w 3 clusters
set.seed(123)
kmeans_TopJeens<-kmeans(TopJeenExprsnScale, centers = 3, nstart=25)
print(kmeans_TopJeens$cluster)

#Visualizing Kmeans
Kpca<-prcomp(TopJeenExprsnScale, scale. = T)
KpcaDta<-data.frame(Kpca$x, Cluster=as.factor(kmeans_TopJeens$cluster))

ggplot(KpcaDta, aes(x=PC1, y=PC2, color=Cluster))+
  geom_point(size=3)+
  labs(title = 'K-means Clustering (PCA visual)',
       x='PC1',
       y='PC2')+
  theme_minimal()

# Hierarchical Clustering


DistMtrx<-dist(TopJeenExprsnScale)
HigherTopJeensExprsn<-hclust(DistMtrx, method = 'complete')
plot(HigherTopJeensExprsn, main = 'Higherarchical Clustering of Top 1000 Most Variable Genes Expressed in Hippocampus Tissue',
     xlab = "",sub = "")
ClusteredTopJeens<-cutree(HigherTopJeensExprsn, k=3)
table(ClusteredTopJeens)

# Comparison of Hierarchical clustering with Meta Data
SSAD$TopJeensClustered<-as.factor(ClusteredTopJeens)

table(SSAD$TopJeensClustered, SSAD$SEX)
table(SSAD$TopJeensClustered, SSAD$AGE)

# Clustering Visualization
NewHippoPCA<-prcomp(TopJeenExprsn, scale. = T)


hPCADta<-data.frame(NewHippoPCA$x[,1:2], Cluster=as.factor(ClusteredTopJeens))
colnames(hPCA)
ggplot(hPCADta, aes(x=PC1,y=PC2, color=Cluster))+
  geom_point(size=3)+
  labs(title = 'Hierarchical Clustering Visual',
       x='PC1',
       y='PC2')+
  theme_minimal()




# Determining Gene's that are driving clusters

HippoLoadings<-HippoPCA$rotation
topGenesPc1<-head(order(abs(HippoLoadings[,1]), decreasing = T),10)
topGenesPc2<-head(order(abs(HippoLoadings[,2]), decreasing = T),10)

TopGeneNames1<-colnames(TopJeenExprsnScale)[topGenesPc1]
TopGeneNames2<-colnames(TopJeenExprsnScale)[topGenesPc2]


# Classification Modeling


# Data Prep


set.seed(123)
Subj100<-sample_n(joined_data1, 100)
GNCol<-setdiff(colnames(Subj100), c("SUBJID", "SEX", "AGE", "DTHHRDY"))

set.seed(123)
SELgns<-sample(GNCol, 50)
SmallHippoDta<-Subj100|>
  select(SUBJID, SEX, AGE, DTHHRDY, all_of(SELgns))

# Possibly Random Forest Modeling
SmallHippoDta$AGE<-as.factor(SmallHippoDta$AGE)

set.seed(123)
trainingIDX<-createDataPartition(SmallHippoDta$AGE, p=0.8, list=F)
trainingDta<-SmallHippoDta[trainingIDX,]
testingDta<-SmallHippoDta[-trainingIDX,]

RFmodel<-randomForest(AGE ~., data = trainingDta, importance=T, ntree=500)
RFPrdct<-predict(RFmodel, newdata = testingDta)
confusionMatrix(RFPrdct, testingDta$AGE)
varImpPlot(RFmodel)

#FullJeenXSC<-cbind(TopJeenExprsnScale,SSAD)
#OversampleJeenAge<-upSample(FullJeenXSC, FullJeenXSC$AGE, yname = 'Age')
#UndersampleJeenAge<-downSample(FullJeenXSC,FullJeenXSC$AGE, yname = 'Age')


# Training(70%) vs Testing Data(30%)
#set.seed(123)
#trainingIDX<-createDataPartition(SSAD$AGE, p=0.7, list=F)
#trainingDta<-[trainingIDX,]
#testingDta<-OversampleJeenAge[-trainingIDX,]
#trainingLbl<-SSAD$AGE[trainingIDX]
#testingLbl<-SSAD$AGE[-trainingIDX]


# Random Forest
#RFmodel<-randomForest(as.factor(AGE)~.,data=trainingDta, ntree=100)
