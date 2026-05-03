#############################################################################################
#################################### MLTEST Enas HASSIS #####################################
############ Multi class analysis (Normal, Small Adenoma, Large Adenoma, Cancer)#############
#############################################################################################


library("dplyr")
library("caret")
filePathOutput = "C://Users//Desktop//RESULTSmulticlass_MLTEST_Enas"
setwd("C:/Users/hassi/Desktop/RESULTSmulticlass_MLTEST_Enas")

################ Read the microbiome species abundance data


AllSamples <- as.data.frame(t(all_samples.sat))
rownames(AllSamples) <- AllSamples$V1
AllSamples$V1 <- NULL



################ create column for the samples ID
AllSamples$ID <- rownames(AllSamples)
AllSamples <- AllSamples[,c(ncol(AllSamples),1:ncol(AllSamples)-1)]
rownames(AllSamples) <- NULL
colnames(AllSamples)[1] <- "ID"

################ order the sample ID
AllSamples <- AllSamples[order(AllSamples$ID, decreasing = FALSE),]
names(AllSamples) <- lapply(AllSamples[1, ], as.character)
AllSamples <- AllSamples[-1,]
colnames(AllSamples)[1] <- "ID"


AllSamples[,2:length(AllSamples)] <-
  as.data.frame(as.matrix((sapply(AllSamples[,2:length(AllSamples)],as.numeric))))

############## Read the patients meta data
library(data.table)
Meta <- meta.tsv

############## biosample_accession and  AllSamples$ID must be identical
Meta <- setDT(Meta)[biosample_accession %chin% AllSamples$ID]
AllSamples$ID == Meta$biosample_accession #### just to verify 
AllSamples <- AllSamples[,-1]

DiagnosisCase <- data.frame(Meta$diagnosis)

####################### Removing the species that contain zeros for all patients ################
####################### for the purpose of P values calculation
AllSamples <- AllSamples[, colSums(AllSamples != 0) > 0]


################################################################################################
########## four cases analysis: (Normal, Small Adenoma, Large Adenoma, Cancer) #################
################################################################################################


########### 1:Normal, 2:Small Adenoma, 3:Large Adenoma, 4:Cancer
DiagnosisCase <-  DiagnosisCase %>%
  mutate(diagnosis = case_when(Meta.diagnosis %in% c("Normal") ~ '1', 
                               
                               Meta.diagnosis %in% c('Small adenoma') ~ '2',
                               
                               Meta.diagnosis %in% c("Large adenoma") ~ '3',
                               
                               Meta.diagnosis %in% c("Cancer") ~ '4'
                               
                               
  ))




####### Now connect species table with daignostic table 
My_data_diagnostic <- bind_cols(AllSamples,DiagnosisCase$diagnosis)
names(My_data_diagnostic)[names(My_data_diagnostic) == '...1704'] <- 'Cases'
My_data_diagnostic$Cases <- as.numeric(My_data_diagnostic$Cases)
My_data_diagnostic <- My_data_diagnostic[order(My_data_diagnostic$Cases, decreasing = FALSE),]  

####### count the number of patients in each group
Normal_COUNT<- nrow(DiagnosisCase[DiagnosisCase$diagnosis== 1,])
SmallADENOMA_COUNT <- nrow(DiagnosisCase[DiagnosisCase$diagnosis== 2,])
LargeADENOMA_COUNT <- nrow(DiagnosisCase[DiagnosisCase$diagnosis== 3,])
CancerCOUNT_COUNT <- nrow(DiagnosisCase[DiagnosisCase$diagnosis== 4,])

##############################################################################################################
############################# Analysis of Variance for ordered factor with 4 levels ###########################
###############################################################################################################

ANOVA <- list()
Pvalues <- NULL
samples <- NULL


for (i in colnames(My_data_diagnostic[,1:(ncol(My_data_diagnostic)-1)]))
{
  set.seed(123)
  
  ANOVA[[i]] <- aov(My_data_diagnostic[, i] ~ ordered(Cases), data = My_data_diagnostic)
  
  Pvalues[[i]] <- summary(ANOVA[[i]])[[1]]$"Pr(>F)"[1] ###### extract p value for each marker

}
Pvalues <- as.data.frame(Pvalues)
Pvalues <- reshape2::melt(Pvalues)
Pvalues$variable <- colnames((AllSamples))
PvaluesNonSig <- Pvalues[Pvalues$value>0.054,]
PvaluesSig <- Pvalues[Pvalues$value<0.055,]


#################################################################################################################
################## Representation of p values; for the significant species and non significant ###################
##################################################################################################################
###################### significant species ordered regarding the level of p values ##############################

pdf("P values of sig Species by ordered ANOVA.pdf",width = 40,height = 10)
ggplot(PvaluesSig, aes(variable,-log10(value))) + 
  geom_point(size = 2, alpha=1) +
  annotate("label",x=5, y = -log10(.001), label = "***",color="black",size=4 , fontface="bold") +
  annotate("label", x = 5, y = -log10(.01), label = "**",color="black",size=4 , fontface="bold") +
  annotate("label", x = 5, y = -log10(.0544), label = "*",color="black",size=4 , fontface="bold")+
  theme(axis.text.x = element_text(size = 10,margin = margin(t = 40),angle = 45, hjust = 1,color = "black"))+
  theme(legend.title=element_text(size=7))+ theme(legend.text = element_text(size = 6))  +
  
  labs(x = "Marker", y = "-log10(p.values)")+
  labs(title=paste0("Significants Species:",nrow(PvaluesSig),"\n","P.values:ANOVA for ordered factor with 4 levels\nNormal: 62, Small Adenoma:23, Large Adenoma:15, Cancer:87\nSignif.codes:'***' 0.001 '**' 0.01 '*' 0.05") )

dev.off()

######### plot nonsignificant species

pdf("P values of NONsig markers by ordered ANOVA.pdf",width = 100,height = 10)
ggplot(PvaluesNonSig, aes(variable,-log10(value))) + 
  geom_point(size = 2, alpha=1) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color = "black"))+
  theme(legend.title=element_text(size=7))+ theme(legend.text = element_text(size = 6))  +
  labs(x = "Marker", y = "-log10(p.values)")+
  labs(title=paste0("NonSignificant Species:",nrow(PvaluesNonSig),"\n","P.values:ANOVA for ordered factor with 4 levels\nNormal: 62, Small Adenoma:23, Large Adenoma:15, Cancer:87\nSignif.codes:'***' 0.001 '**' 0.01 '*' 0.05") )

dev.off()

###########################################################################################
########## visualizing the p values across the 4 groups by violin plot ####################
########## the visualization for the higher significant species ###########################
########## when P value <.001, but you can change the thresold of p values ################


Pvalue_Extreme<- Pvalues[Pvalues$value<.001,]


Pvalue_Extreme$variable <-  sub("%", "", (Pvalue_Extreme$variable)) 
cols <- (Pvalue_Extreme$variable)

library(vioplot)
for (i in c(cols))
{
  png(paste0(i,"violin plot.png"), width = 900,height = 800)
  print(vioplot(My_data_diagnostic[,i] ~ My_data_diagnostic$Cases,
                ##col = 1:length(levels(My_data_diagnostic$Cases)),
                xlab = "cases", ylab = i))
  dev.off()
}


