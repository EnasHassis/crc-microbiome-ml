# This project demonstrates the use of machine learning models to extract microbiome-based signatures for colorectal cancer prediction.
# Microbiome-based classification of colorectal cancer

## Project overview
This project applies machine learning to microbiome species abundance data to classify colorectal cancer (CRC) status (Cancer vs Non-Cancer).

## Methods
- Logistic Regression
- Random Forest
- Support Vector Machine (SVM)

## Workflow
1. Data preprocessing (normalization, filtering)
2. Feature selection
3. Model training and evaluation (cross-validation)
4. Model comparison

## Results
- Best model: Random Forest
- ROC-AUC: .87

## Key findings
- Identified microbial species associated with CRC
- Machine learning can effectively classify disease status

## Project structure
- notebooks/: analysis notebooks
- src/: scripts
- results/: figures and outputs

## Note
Dataset not included due to usage restrictions.

## Results visualization

![ROC Curve](results/roc.png)
![Feature Importance](results/importance.png)
