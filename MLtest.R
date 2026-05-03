#############################################################################################
#################################### MLTEST Enas HASSIS #####################################
################ Build ML classification algorithm that can predict patient status ##########
###############from their relative species abundances, mainly Cancer or not Cancer ##########
#############################################################################################

library("dplyr")
library("caret")
filePathOutput = "C://Users//Desktop//RESULTSmltest"
setwd("C:/Users/hassi/Desktop/RESULTSmltest")

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

############## Read the patients meta data
library(data.table)
Meta <- meta.tsv

############## biosample_accession and  AllSamples$ID must be identical
Meta <- setDT(Meta)[biosample_accession %chin% AllSamples$ID]
AllSamples$ID == Meta$biosample_accession #### just to verify 
AllSamples <- AllSamples[,-1]



###########################################################################################
######################### Cancer/notcancer classification step ############################
###########################################################################################

Case <- data.frame(Meta$diagnosis)

############  Here I treated all these cases as notcancer ; Normal, Small Adenoma, Large Adenoma 
############ give Cancer: 1  notcancer: 0
Case$Diagnosis <- ifelse(Case$Meta.diagnosis =="Cancer",1,0)
LogicalCol <- data.frame(Case$Diagnosis)
names(LogicalCol)[names(LogicalCol) == "Case.Diagnosis"] <- "Diagnosis"
Case$Diagnosis==LogicalCol    #### to verify

################ connect Allsamples data to the diagnostic column and put them in another table
My_data_diagnostic <- bind_cols(AllSamples,LogicalCol)


#############################################################################################
################# seperate the data into 2 tables Cancer/Notcancer ##############################
#############################################################################################

InfectedTable <- My_data_diagnostic[My_data_diagnostic$Diagnosis== 1,]
NoninfectedTable <- My_data_diagnostic[My_data_diagnostic$Diagnosis== 0,]

InfectedTable2 <- subset (InfectedTable, select = -Diagnosis)
NoninfectedTable2 <- subset(NoninfectedTable,select=-Diagnosis)

##########################################################################################
######################### define titles for all plots ####################################
##########################################################################################

DataName <- "Microbiome Species Abundance Data"  
NumMarker <- paste0(length(AllSamples)," Microbiome species")
NumPatients <- paste0(nrow(AllSamples)," Patients")
NumNash <- paste0(nrow(InfectedTable)," Cancer")
NumNonnash <- paste0(nrow(NoninfectedTable)," NotCancer")



##########################################################################################
############################# Build  machine learning models #############################
##########################################################################################

##########################################################################################
########## Random forest using train function 5 folds cross validation ##################
##########################################################################################

library("caret")
set.seed(1) 
modelLookup("rf")

########### generate the model parameters
trainforest <- caret::trainControl(method="cv",  #### resampling method CV
                                   summaryFunction=twoClassSummary,# Use AUC to pick the best model
                                   classProbs=TRUE,number = 5,  # 5 folds CV
                                   savePredictions = TRUE)


########## convert all columns to numeric 
My_data_diagnostic <-  lapply( My_data_diagnostic, function(x) as.numeric(as.character(x)))
My_data_diagnostic <- as.data.frame(My_data_diagnostic)

########### change the level of diagnosis
My_data_diagnostic$Diagnosis <- as.factor(My_data_diagnostic$Diagnosis)
levels((My_data_diagnostic$Diagnosis))
levels(My_data_diagnostic$Diagnosis) <- c("No","Yes")


########### fit predictive models 
fitforest<- caret::train(Diagnosis~.,data=My_data_diagnostic,
                         method = "rf",   
                         metric="ROC",
                         trControl=trainforest )


########## plot the model
pdf("Random Forest By applying cross validation.pdf")
plot(fitforest,main=paste0("ROC: ",DataName,",",NumMarker,";",NumPatients,
                           "\n Model:Random forest model by applying cross validation 5 folds"))

dev.off()

########## Save the markers (microbiome species) important
write.csv2(as.matrix(varImp(fitforest)$importance,scale=FALSE), row.names = TRUE,
          "Random Forest; All Variables Important Table .csv")


########## plot the most important variables, here i selected the first 40 var
library(ggplot2)
pdf(" Random Forest, Variables Importance; 40 Variables.pdf",height = 12,width = 10)
plot(varImp(fitforest), top = 40,pch=16,main=paste0("Variables Important: ",DataName,"\n",NumMarker,";",NumPatients,
                                                   "\n Model:Random forest model by applying cross validation 5 folds")) 
dev.off()

################## AUC for RANDOM FOREST model using evalm function 
library(MLeval)
evalRandom<- evalm(fitforest,percent = 95, showplots = TRUE, positive = NULL, plots ="r")
ci <- evalRandom$optres$`Group 1`$CI[13]

pdf("Random Forest by applying cross validation, AUC.pdf")
evallog<- (evalm(fitforest,percent = 95, showplots = TRUE, positive = NULL, plots ="r",
                 title = paste0("AUC: ",DataName,"\n"
                                ,NumMarker,";",NumPatients,"\n Random Forest model by applying cross validation 5 folds \n"
                 ),labs(caption=paste0("CI: ",ci))))
dev.off()






##########################################################################################
########## Logistic Regression using train function 5 folds cross validation #############
##########################################################################################

set.seed(1)
trainlogistic <- caret::trainControl(method="cv", summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                                     classProbs=TRUE,number = 5,savePredictions = TRUE)

fitlogistic<- caret::train(Diagnosis~.,data=My_data_diagnostic,
                           method = "glm", 
                           family="binomial",
                           metric="ROC",
                           trControl=trainlogistic,preProcess="scale")

################## AUC for logistic model by evalm function 
library(MLeval)

evallog<- evalm(fitlogistic,percent = 95, showplots = TRUE, positive = NULL, plots ="r")
ci <- evallog$optres$`Group 1`$CI[13]


pdf("Logistic Regression by applying cross validation, AUC.pdf")
evallog<- (evalm(fitlogistic,percent = 95, showplots = TRUE, positive = NULL, plots ="r",
                 title = paste0("AUC: ",DataName,"\n",NumMarker,";",NumPatients,
                                "\n Miltivariete LR model by applying cross validation 5 folds \n"),
                 labs(caption=paste0("CI: ",ci))))
dev.off()



###########################################################################################
############## logistic regression with lasso penalty: CV.glmnet function ##################
############################################################################################

################ Find the best lambda using cross-validation

set.seed(1) 
x <- model.matrix(Diagnosis~., My_data_diagnostic)[,-1]
y <-My_data_diagnostic$Diagnosis

library(glmnet)
cv.lasso <- cv.glmnet((x), (y), alpha = 1, family = "binomial",
                      type.measure  = "auc", nfolds = 5)

cv.lasso$lambda

pdf("Lasso Models,cv.glmnet.pdf")
plot(cv.lasso)
title( paste0("Lasso Model using cv.glmnet function(CV;5folds)\n",DataName,",",NumMarker,";",NumPatients),
       line=2,cex=.3)
dev.off()

##################### cvm here contains all auc for all models and maximum is the best one 
(cv.lasso)
AUC <- max(cv.lasso$cvm)
cv.lasso$lambda.1se
cv.lasso$lambda.min
#The sign of coefficients tells you if the independent variable
#is positively or negatively related to the outcome 

######## coef with LAMBDA 1SE
coef.lassomodel <- as.matrix(coef(cv.lasso)) 


write.csv2(coef.lassomodel, row.names = TRUE,
           "Lasso Model,1se;coefficient.csv")

######## coef with lambda min

coefficientmin<- coef(cv.lasso,cv.lasso$lambda.min) 
write.csv2(as.matrix(coefficientmin),
           "Lasso Model,lambdamin;coefficient.csv"
                  ,row.names = TRUE)



##########################################################################################
#####################  lasso model: train function 5 folds CV ############################
##################### apply alpha=1 and lambda generated from cv.glmnet ##################
##########################################################################################

modelLookup("glmnet")
set.seed(1) 
library(caretStack)
trainlasso <- caret::trainControl(method="cv",  
                                  summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                                  classProbs=TRUE,number = 5,savePredictions = TRUE)

fitlasso <- caret::train(Diagnosis~.,data=My_data_diagnostic,
                         method = "glmnet", 
                         metric="ROC",preProcess="scale",
                         trControl=trainlasso,
                         tuneGrid = expand.grid(alpha = 1,lambda = cv.lasso$lambda))

fitlasso$results$lambda <- sprintf("%.4f", fitlasso$results$lambda) 
coefficientscaret<-coef(fitlasso$finalModel,fitlasso$bestTune$lambda)
write.csv2(as.matrix(coefficientscaret),
           "Lasso Model,CARET;coefficient.csv"
                  ,row.names = TRUE)


########### all lambda or level of penalization from 2 methods
cv.lasso$lambda.1se
cv.lasso$lambda.min
fitlasso$bestTune$lambda
############### save these lambdas
cat(paste0("cv.lasso$lambda.1se\n",cv.lasso$lambda.1se,"\ncv.lasso$lambda.min",
           cv.lasso$lambda.min,"\nfitlasso$bestTune$lambda",
           fitlasso$bestTune$lambda), file = "Lasso,lambda from 2 packeges.txt")

################## plot the model
pdf("lasso model, caret package.pdf")
ggplot(fitlasso)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(text = element_text(size=8))+
  ggtitle(paste0("lasso Model using caret package:",DataName,";",NumMarker,",",NumPatients))
dev.off()

################## AUC for LAsso model by evalm function 

evallasso <- evalm(fitlasso,percent = 95, showplots = TRUE, positive = NULL, plots ="r")
ci <- evallasso$optres$`Group 1`$CI[13]


pdf("Lasso Model using cross validation(caret package), AUC.pdf")
evallog <- (evalm(fitlasso,percent = 95, showplots = TRUE, positive = NULL, plots ="r",
                  title = paste0("AUC: ",DataName,"\n" ,NumMarker,";",NumPatients,"\n Lasso model using caret package(CV;5folds) \n"
                  ),labs(caption=paste0("CI: ",ci) )))
dev.off()

########## plot the most important variables, here i selected the first 40 var
pdf("Lasso, Variables Importance; 40 Variables.pdf",height = 12,width = 10)
plot(varImp(fitlasso), top = 40,pch=16) 
dev.off()







##########################################################################################################
##################  secondary objective: representation of the high important species ####################
#################  involved in the signature of CRC in the form of a network  ###########################
##########################################################################################################
################### simple network based on the most significant species in RF model######################
library(igraph)

########## Determine the scores of the species(markers)
important_markers <- varImp(fitforest, scale = FALSE)$importance
important_markers <- as.data.frame(important_markers)
important_markers$Name <- rownames(important_markers)


########## Connect markers with higher importance scores in RF model; select the higher imp species
threshold_importance <- 0.5  ######## or other threshold
high_importance_markers <- important_markers$Name[important_markers$Overall > threshold_importance]

######### Create an empty graph under selected markers length
marker_graph <- graph.empty(n = length(high_importance_markers), directed = FALSE)

########## Add nodes; vertices (markers,species)
V(marker_graph)$name <- high_importance_markers


########### fill the graph with nodes and edges
for (i in 1:(length(high_importance_markers) - 1)) {
  for (j in (i + 1):length(high_importance_markers)) {
    marker_graph <- add_edges(marker_graph, c(high_importance_markers[i], high_importance_markers[j]))
  }
}


########### Plot the network with nodes colored based on importance scores in RF model
pdf("Network of Important Species, simple.pdf")

plot(marker_graph, layout = layout_with_fr,main = "Network of Important Species (Importance Scores From RF Model)", 
     vertex.color = heat.colors(nrow(important_markers))[rank(-important_markers$Overall)],
     vertex.label.cex = .7,vertex.label.family = "sans",vertex.label.font = 1)
dev.off()


############################################################################################
##############Network by calculating the distance and generating adjacencey matrix #########
############################################################################################
############# More complicated Network based on the distance betweeen species###############



########## Extract importance scores from the RF model
importance_scores <- varImp(fitforest, scale = FALSE)$importance
importance_scores$Name <- rownames(importance_scores)

########### Sort species by importance scores
sorted_species <- importance_scores[order(importance_scores$Overall, decreasing = TRUE), ]

######## Select the top 40 species (higher important markers)
top_species <- row.names(sorted_species)[1:40]

####### Extract importance scores for top species
top_species_importance <- sorted_species$Overall[1:40]


###### Calculate distances based on importance scores of top species
top_species_distances <- as.matrix(dist(top_species_importance))

########## Set a threshold for distances, here i selected 0.5
distance_threshold <- 0.5

######### Create the graph
species_network <- graph.empty(n = length(top_species), directed = FALSE)

########### Add vertices (species)
V(species_network)$name <- top_species
vertex_mapping <- 1:length(top_species)

######### Create adjacency matrix based on distance threshold
adjacency_matrix <- (top_species_distances <= distance_threshold)
adjacency_matrix <- as.matrix(adjacency_matrix)

######### Convert adjacency matrix to edges
edges <- which(adjacency_matrix, arr.ind = TRUE)

######### Add edges to the graph
species_network <- add_edges(species_network, edges)

########## Plot the network
pdf("Network of Important Species,based on disatnce.pdf")
plot(species_network, layout = layout_with_fr,
     vertex.size = 3, edge.arrow.size = 0.2,
     vertex.label.dist = 1, vertex.label.cex = 0.8,
     edge.color = "gray", edge.lty = "solid",main = "Network of Higher 40 Important Species Based on RF model")

dev.off()






