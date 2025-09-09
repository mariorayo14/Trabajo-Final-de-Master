rm(list=ls())


####### 1. Establecimiento de la ruta de directorio de trabajo########
getwd()
path <- "C:/Users/Asus/Desktop/TFM/Scripts"
setwd (path)

.libPaths("C:/Users/Asus/Desktop/TFM/Scripts/R_libreria")

save.image("C:/Users/Asus/Desktop/TFM/Scripts/Script_TFM.RData")

load("C:/Users/Asus/Desktop/TFM/Scripts/Script_TFM.RData")






##########       2. Instalación y carga de librerias  ##############################

#Instalación de librerias

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager")}
BiocManager::install("TCGAbiolinks")
BiocManager::install("org.Hs.eg.db") #Paquete para cambiar los ENSEMBL id al SYMBOL GENE
BiocManager::install("DESeq2")      #Paquete para la normalizacion y realización de los analisis de DEGs
BiocManager::install("apeglm")
if (!require("clusterProfiler")) {
  BiocManager::install("clusterProfiler")}
BiocManager::install("EnhancedVolcano")
install.packages("janitor")
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
packageVersion("DESeq2")


install.packages("dplyr")         # Paquete para manipulación de datos
install.packages("SummarizedExperiment")  # Paquete para almacenar datos omics
install.packages("GEOquery")      # Paquete para descargar datos de GEO (Gene Expression Omnibus)
install.packages("tidyverse")
install.packages("RColorBrewer")
install.packages("ggrepel")
install.packages("patchwork")
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)
version
packageVersion("TCGAbiolinks")

install.packages("patchwork", dependencies = TRUE)
install.packages ("ggplot2")
install.packages("smotefamily")
install.packages ("glmnet")
install.packages ("caret")
install.packages ("PRROC")
install.packages ("rattle")
install.packages("rlang")
#Carga de librerias
library(ggplot2)
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library (AnnotationDbi)
library(DESeq2)
library (clusterProfiler)
library (janitor)
library(tidyverse)
library(EnhancedVolcano)
library(apeglm)
library(RColorBrewer) 
library(ggrepel)
library(patchwork)# for nice annotations
library(enrichplot)
library(VennDiagram)
library(biomaRt)
library(glmnet)
library(caret)
library(smotefamily)
library(rpart) # DT
library(rpart.plot) # DT plot
library(rattle) # DT plot
library(pROC) # ROC
library(PRROC) # PR-Curve
library(rlang)
######### 3. Descarga de los datos clinicos de pacientes TNBC y tejidos normales #####
# Datos clinicos de los pacientes con cancer de mama
clin_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
View(getResults(clin_query))
# Descargar los datos clinicos
#GDCdownload(clin_query)

# Preparar los datos clinicos
clin_data <- GDCprepare(clin_query)

#Seleccionamos la tabla con la informacion clinica de los pacientes
patient_data <- clin_data$clinical_patient_brca


# Filtrar pacientes con cáncer de mama triple negativo
# Los criterios son: ER negativo, PR negativo, HER2 negativo
triple_negative_patients <- patient_data [
  patient_data$er_status_by_ihc== "Negative" &
    patient_data$pr_status_by_ihc == "Negative" &
    patient_data$her2_status_by_ihc == "Negative",
]

#Guardamos los identificadores de los pacientes tnbc para la expresion genica
tnbc_barcodes <- triple_negative_patients$bcr_patient_barcode

##### Datos clinicos de los tejidos normales ######
#No lo he conseguido
clin_query_normal <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab",
  sample.type = "Solid Tissue Normal"
)


###### 4. Descarga de los datos transcriptomicos ######

########################    DATOS DE EXPRESIÓN TNBC #####################
Expresion_data_tnbc <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type =  "Primary Tumor",
  barcode = tnbc_barcodes
)

# Download the data
#GDCdownload(Expresion_data_tnbc)

# Prepare the data for analysis
tnbc_expression_data <- GDCprepare(Expresion_data_tnbc)

#Guardamos el RDS en la carpeta RDS
#saveRDS(tnbc_expression_data, "C:/Users/Asus/Desktop/TFM/Scripts/RDS/tnbc_expresion_data_SummarizedExperiment.rds")

#Leemos el documento RDS
TNBC_RDS<- readRDS("C:/Users/Asus/Desktop/TFM/Scripts/RDS/tnbc_expresion_data_SummarizedExperiment.rds")

# Nos guardamos solo las raw counts matrix
TNBC_mat<- assay(TNBC_RDS)

######### DATOS DE EXPRESION DE TEJIDOS NORMALES   ######
Expresion_data_sanos <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal",
)

#Descargamos la expresion en formato GDC
#GDCdownload(Expresion_data_sanos)

#Guardamos la expresion
expresion_sanos <- GDCprepare(Expresion_data_sanos)

#Guardamos el RDS en la carpeta RDS

#saveRDS(expresion_sanos, "C:/Users/Asus/Desktop/TFM/Scripts/RDS/expresion_sanos_SummarizedExperiment.rds")

Sanos_RDS<- readRDS("C:/Users/Asus/Desktop/TFM/Scripts/RDS/expresion_sanos_SummarizedExperiment.rds")

# Nos guardamos solo las raw counts matrix
Sanos_mat<- assay(Sanos_RDS)


######### 5.  Dataset conjunto con la información Coldata de los datos de expresion del tumor
#TNBC
coldata_tnbc = colData(tnbc_expression_data)
coldata_tnbc = as.data.frame(coldata_tnbc)

#NORMAL
coldata_sano = colData(expresion_sanos)
coldata_sano = as.data.frame(coldata_sano)

########################     DATAFRAME  CON LOS COUNTS BRUTOS DE LOS DATOS DE EXPRESION ########################
#TNBC
rawcounts_tnbc_matrix = assay(tnbc_expression_data) 
rawcounts_tnbc_df = as.data.frame(rawcounts_tnbc_matrix)
# SANO
rawcounts_sanos_matrix = assay(expresion_sanos)
rawcounts_sanos_df = as.data.frame(rawcounts_sanos_matrix)



####### 6. Cambiamos los simbolos ENSEMBL a los nombres de los GENES ########
#TNBC

# Eliminamos el numero de la version final del ENSEMBL ID
rawcounts_tnbc_matrix_genes<- rownames(rawcounts_tnbc_matrix) %>%
  tibble::enframe() %>%
  mutate(ENSEMBL =stringr::str_replace(value, "\\.[0-9]+", ""))

head(rawcounts_tnbc_matrix_genes)

## Hay gene symbols duplicados para diferentes ENSEMBL ID
clusterProfiler::bitr(rawcounts_tnbc_matrix_genes$ENSEMBL, 
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db) %>%
  janitor::get_dupes(SYMBOL) %>%
  head()

# Nos quedamos con uno solo
TNBC_gene_map<- clusterProfiler::bitr(rawcounts_tnbc_matrix_genes$ENSEMBL,
                                           fromType = "ENSEMBL",
                                           toType = "SYMBOL",
                                           OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE) 

#Unimos la matriz de genes con la de los gene simbols
TNBC_gene_map<- TNBC_gene_map %>%
  left_join(rawcounts_tnbc_matrix_genes)

head(TNBC_gene_map)

#Transformamos la matriz de expresion genica reemplazando los rownames por los gene symbols
TNBC_mat<- rawcounts_tnbc_matrix[TNBC_gene_map$value, ]
row.names(TNBC_mat)<- TNBC_gene_map$SYMBOL


#Normal
#Al ser los mismos genes solo tenemos que hacer el paso final
Normal_mat<- rawcounts_sanos_matrix[TNBC_gene_map$value, ]
row.names(Normal_mat)<- TNBC_gene_map$SYMBOL

#Check que los genes son los mismos
all.equal(rownames(Normal_mat), rownames(TNBC_mat))

########## 7. Preprocesamiento para normalización #########
#Convertimos la matriz de expresion de TNBC a dataframe para poder unir la columna de Estadio
TNBC_mat_df <- as.data.frame(t(TNBC_mat))

# Añadir la columna de Estadio como primera posición
TNBC_df <- cbind(Estadio = triple_negative_patients$ajcc_pathologic_tumor_stage, TNBC_mat_df)
# Hacer una copia de seguridad de tus datos
TNBC_backup <- TNBC_df

# Pasar los stage a I. IIA,IIB III y IV
TNBC_df <- TNBC_df %>%
  mutate(Estadio = case_when(
    Estadio == "Stage I" | Estadio == "Stage IA" | Estadio == "Stage IB" ~ "I",
    Estadio == "Stage II" | Estadio == "Stage IIA" | Estadio == "Stage IIB" ~ "II",
    Estadio == "Stage III" | Estadio == "Stage IIIA" | Estadio == "Stage IIIB" | Estadio == "Stage IIIC" ~ "III",
    Estadio == "Stage IV" ~ "IV",
    TRUE ~ "No disponible"  #Los valores NA pasan a llamarse No disponible
  ))
#Añadimos la columna Muestra = TNBC)
TNBC_df <- cbind(Muestra = "TNBC", TNBC_df)

#Añadimos la columna Fase para temprana o tardía (Estadio I y II es para temprana, Estadio III Y IV para fase tardía)
TNBC_df <- TNBC_df %>%
  mutate(Fase = case_when(
    Estadio %in% c("I", "II") ~ "Temprana",
    Estadio %in% c("III", "IV") ~ "Tardía",
    TRUE ~ NA_character_  # para el resto 
  ))
#La quiero en la tercera posicion
TNBC_df <- TNBC_df %>%
  relocate(Fase, .after = Estadio)

#Eliminamos las filas que no tienen disponible el Estadio
TNBC_df_purificado <- TNBC_df[TNBC_df$Estadio != "No disponible", ]

#Ahora lo mismo con los normal
Normal_mat_df <- as.data.frame(t(Normal_mat))
# Añadir las columnas de Muestra, Estadio y Fase como primeras posiciones (todas "Normal")
Normal_df <- cbind (Muestra = "Normal", Estadio = "Normal",Fase ="Normal", Normal_mat_df)
#Unimos los dos dataframes
TNBC_normal_df <- rbind (TNBC_df_purificado, Normal_df)

######## 8. Analisis de Expresion Diferencial Génica por DESeq2 de las Fases Normales, Tardías y Tempranas ########
# Para normalizar necesitamos una matriz de rawcounts y un dataframe con los datos en columnas 

# Cogemos las primeras dos columnas con los metadatos
colData_df <- TNBC_df_purificado[, c("Muestra", "Estadio", "Fase")]
# Supongamos que tienes el objeto con datos clínicos o de expresión
rownames(colData_df)  # Esto te da los IDs de las muestras, tipo: "TCGA-XX-XXXX-01A"

# Si tienes el objeto dds (DESeq2)
colnames(dds)
# Obtener los nombres de las columnas
nombres_columnas <- colnames(dds)

# Opción 1: Usar el paquete flextable
library(flextable)
library(officer)

# Crear una tabla con los nombres
tabla_nombres <- data.frame(Columnas = nombres_columnas)
ft <- flextable(tabla_nombres)

# Guardar en Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "nombres_columnas.docx")

#Para el analisis con DESeq2, es importante pasar a factor la metadata y poner Normal como referencia (el primero)
colData_df$Muestra <- factor(colData_df$Muestra, levels = c("TNBC"))
colData_df$Estadio <- factor(colData_df$Estadio, levels = c("I", "II", "III", "IV"))
colData_df$Fase <- factor(colData_df$Fase, levels = c("Temprana", "Tardía"))

rawcounts_df <- TNBC_df_purificado[, -(1:3)]  # eliminar muestras, estadio y fase
rawcounts_df <- t(rawcounts_df)  # transponer: DESeq2 espera genes en filas y muestras en columnas

#Comprobación que el nombre de las filas de colData es el mismo que las columnas de rawcounts
all.equal(rownames(colData_df),colnames(rawcounts_df))

####En este caso haremos el Análisis de Expresión Diferencial únicamente para la fase temprana y tardía

#DESeq2 para estadio del cáncer
dds <- DESeqDataSetFromMatrix(countData = rawcounts_df, 
                              colData = colData_df,
                              design = ~ Fase)

# # Filtrado opcional: eliminar genes con muy pocos counts
dds_keep <- dds[rowSums(counts(dds)) > 10, ]  # por ejemplo, al menos 10 reads en total por gen


# Ejecución de DESeq2
dds <- DESeq(dds_keep)

res_temprana_tardia <- results(dds, contrast = c("Fase", "Tardía", "Temprana"))
# Volcano plot para cada comparación
plot_temprana_tardia <- EnhancedVolcano(res_temprana_tardia,
                                     lab = rep("", nrow(res_temprana_tardia)),
                                     x = 'log2FoldChange',
                                     y = 'padj',
                                     title = 'Estadio Tardío vs Temprano',
                                     subtitle = "",
                                     pCutoff = 0.05,
                                     FCcutoff = 1.0,
                                     pointSize = 2.5,
                                     labSize = 1.0)

plot_temprana_tardia
# Guardar en PNG
ggsave("volcano_plot_temprana_tardia.png",
       plot = plot_temprana_tardia,
       width = 18, height = 10, dpi = 300)



# Listas de genes significativos (p-adj < 0.05 y |log2FC| > 1)
DEGs_significativos <- rownames(subset(res_temprana_tardia, padj < 0.05 & abs(log2FoldChange) > 1))
length(DEGs_significativos)

#Exportación en formato .csv
write.csv(DEGs_significativos, "DEGs_significativos.csv", row.names = FALSE)

#Guardamos los nombres de los genes para el analisis funcional
writeLines(DEGs_significativos, "DEGs_significativos.txt")


#Top10 upregulated y downregulated
# Convertir el objeto DESeq2 'results' a data frame
res_df <- as.data.frame(res_temprana_tardia)
res_df$gene <- rownames(res_df)

# Filtrar genes significativos (opcional pero recomendado)
res_sig <- subset(res_df, padj < 0.05)

# Ordenar por log2FoldChange descendente (mayor expresión en Temprana)
top10_up <- res_sig[order(-res_sig$log2FoldChange), ][1:10, ]

# Ordenar por log2FoldChange ascendente (mayor expresión en Tardía)
top10_down <- res_sig[order(res_sig$log2FoldChange), ][1:10, ]


# Creación de una tabla con los 3 grupos con los 20 DEGs y su FC
tabla_final <- data.frame(
  DEG_up = rownames(top10_up),
  FC_up   = round(top10_up$log2FoldChange, 2),
  DEG_down  = rownames(top10_down),
  FC_down   = round(top10_down$log2FoldChange, 2)
)


print(tabla_final)

install.packages("officer")
install.packages("flextable")
library(flextable)
library(officer)

# Crear tabla
tabla_word <- flextable(tabla_final)

# Guardar en Word
doc <- read_docx()
doc <- body_add_flextable(doc, value = tabla_word)
print(doc, target = "tabla_editable.docx")


##### Filtro de genes no codificantes ######
#Conexión a Ensembl y el dataset
# Conexión a Ensembl
library(biomaRt)

# Usar un servidor espejo de Ensembl (por ejemplo el de EE. UU.)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")  # También puedes probar "asia" o "uswest"

# Lista de DEGs
gene_list <- DEGs_significativos  # formato SYMBOL

# Obtener información de tipo de gen
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = gene_list,
               mart = ensembl)

# Filtrar solo genes codificantes para proteínas
genes_codificantes <- annot[annot$gene_biotype == "protein_coding", "hgnc_symbol"]

#Guardamos los nombres de los genes codificantes para el analisis funcional
writeLines(genes_codificantes, "DEGs_codificantes.txt")

library(org.Hs.eg.db)
library(dplyr)

# Tus genes DEGs (en SYMBOL)
genes_symbol <- DEGs_significativos 

# Mapeo a ENTREZID y tipo de gen
gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = genes_symbol,
                                   columns = c("SYMBOL", "GENETYPE", "ENTREZID"),
                                   keytype = "SYMBOL")

# Filtrar por los que codifican para proteínas
genes_codificantes2 <- gene_info %>%
  filter(GENETYPE == "protein-coding" & !is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

genes_codificantes2 <- genes_codificantes2$SYMBOL

writeLines(genes_codificantes2, "DEGs_codificantes2.txt")

Genes_up <- rownames(subset(res_sig[res_sig$log2FoldChange > 1, ]))# Sobreexpresados en TARDÍA
Genes_down <- rownames(subset(res_sig[res_sig$log2FoldChange < -1, ]))  # Sobreexpresados en TEMPRANA

#Exportación en formato .csv
write.csv(Genes_up, "Genes_up.csv", row.names = FALSE)
write.csv(Genes_down, "Genes_down.csv", row.names = FALSE)


#Guardamos los nombres de los genes para el analisis funcional
writeLines(Genes_up, "Genes_up.txt")
writeLines(Genes_down, "Genes_down.txt")

###### Machine Learning ######
#### Preprocesamiento de datos #####

#Nuestro dataset inicial es dds_keep
#Normalizamos el dataset
vsd <- vst(dds_keep, blind = FALSE)

#Transformamos los datos de genes x muestras en una matriz
expr_matrix <- assay(vsd)  

# Filtrar solo los genes de interés 
expr_DEGs <- expr_matrix[DEGs_significativos, ]

# Transponemos la matriz
expr_DEGs_t <- t(expr_DEGs)

#Escalamos la matriz
scaled_expr_DEGs <- scale(expr_DEGs_t)

#Convertimos las etiquetas en binario (Tardía = 1; Temprana = 0)
colData_df$Fase_bin <- ifelse(colData_df$Fase == "Tardía", 1, 0)
colData_df$Fase_bin 

#Comprobamos que son las mismas filas 
all.equal(rownames(colData_df),rownames(scaled_expr_DEGs))

#Definimos las variables para LASSO
x <- scaled_expr_DEGs
y <- as.numeric(colData_df$Fase == "Tardía")

# Aseguramos que x sea numérica
x <- apply(x, 2, as.numeric)

# Seleccionamos características usando Lasso
#Para la reproductibilidad
set.seed(123)

# LASSO con validación cruzada, binomial al haber solo dos tipos 
cv_lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")

# En lugar de usar lambda.min (el mínimo error), prueba con lambda.1se
coef_lasso <- coef(cv_lasso, s = "lambda.1se")
selected_genes <- rownames(coef_lasso)[which(coef_lasso != 0)][-1]

selected_genes
ml_data_lasso <- as.data.frame(x[, selected_genes])
ml_data_lasso$Fase <- y

##### División del set ####
# Dividir el conjunto de datos en conjuntos de entrenamiento y prueba
set.seed (123) #Determinamos la semilla de aleatoriedad

#Aseguramos que la variable de salida sea factor
ml_data_lasso$Fase <- as.factor(ml_data_lasso$Fase)

#Comprobación de desbalance
table (ml_data_lasso$Fase)

#División en entrenamiento y set (70/30)
trainIndex <- createDataPartition(ml_data_lasso$Fase, p = 0.7, list = FALSE)
train_data <- ml_data_lasso[trainIndex, ]
test_data <- ml_data_lasso[-trainIndex, ]

#Para SMOTE, se necesita una data unicamente numérica (X)
X <- train_data[, sapply(train_data, is.numeric)]

#Aplicación de SMOTE al conjunto de entrenamiento (por desbalance de muestras)
train_data_smote <- SMOTE(X = X, 
                          target = train_data$Fase,
                          K = 5, dup_size = 2)

train_data_smote <- train_data_smote$data
train_data_smote$Fase <- as.factor(train_data_smote$class)
train_data_smote$class <- NULL

# Crear la formula sumando cada gen
length(train_data_smote) #Para saber cuantas variables se encuentran
names <- colnames(train_data_smote[1:28]) #La variable Fase no se incluye
formula <- as.formula(paste("Fase ~", paste(names, collapse = "+"))) #Se obtiene la formula del diagnosis junto a la suma de todos los parametros
formula

str(train_data_smote)
table (train_data_smote$Fase)

colnames(train_data_smote) <- make.names(colnames(train_data_smote))
colnames(test_data) <- make.names(colnames(test_data))

train_data_smote$Fase <- factor(train_data_smote$Fase, levels = c("0", "1"), labels = c("Temprana", "Tardia"))
test_data$Fase <- factor(test_data$Fase, levels = c("0", "1"), labels = c("Temprana", "Tardia"))


#### ---- kNN ----
# Crear un modelo de k-NN utilizando el paquete caret
knnModel <- train(Fase ~ .,
                  data = train_data_smote,
                  method = "knn",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("center", "scale"),
                  tuneLength = 30)
knnModel

plot(knnModel)

#Realizar predicciones en los datos de test con el modelo
knnmodel_prediction <- predict (knnModel, newdata = test_data)
knnmodel_prediction 

# Luego llamar a confusionMatrix
knnmodel_confusionMatrix <- confusionMatrix(knnmodel_prediction, test_data$Fase)
knnmodel_confusionMatrix

#Probabilidades
probabilities_knnmodel <- predict(knnModel, newdata = test_data, type = "prob")
probabilities_knnmodel

#Metrica F1
f1_score_knn <- as.numeric(knnmodel_confusionMatrix$byClass['F1'])
f1_score_knn
#### ---- Support Vector Machine ----
# Crear un modelo de SVM lineal utilizando el paquete caret
# parámetro C "cost" por defecto es 1, pero puedes tunearlo. Controla la flexibilidad del modelo para encontrar un equilibrio entre un margen amplio y la clasificación correcta de las muestras
svmModelLineal <- train(Fase ~.,
                        data = train_data_smote,
                        method = "svmLinear",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(C = seq(0, 2, length = 20)), #C grande lleva al sobreajuste, C pequeño al infraajuste
                        prob.model = TRUE) 
svmModelLineal

plot(svmModelLineal)

# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
svmModelLineal_prediction <- predict(svmModelLineal, newdata = test_data)
svmModelLineal_prediction

# Evaluar la precisión del modelo utilizando la matriz de confusión
svmModelLineal_confusionMatrix <- confusionMatrix(svmModelLineal_prediction, test_data$Fase)
svmModelLineal_confusionMatrix
# SVM lineal
probabilities_svm_linear <- predict(svmModelLineal, newdata = test_data, type = "prob")
probabilities_svm_linear

f1_score_svm_linear<- as.numeric(svmModelLineal_confusionMatrix$byClass['F1'])
f1_score_svm_linear


#Crear un modelo de SVM tipo kernel utilizando el paquete caret
# no hace falta tunear el parámetro C "cost" 
svmModelKernel <- train(Fase ~.,
                        data = train_data_smote,
                        method = "svmRadial",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneLength = 10,
                        prob.model = TRUE) 
svmModelKernel

plot(svmModelKernel)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
svmModelKernel_predictions <- predict(svmModelKernel, newdata = test_data)
svmModelKernel_predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
svmModelKernel_confusionMatrix <- confusionMatrix(svmModelKernel_predictions, test_data$Fase)
svmModelKernel_confusionMatrix
# SVM kernel
probabilities_svm_kernel <- predict(svmModelKernel, newdata = test_data, type = "prob")
probabilities_svm_kernel



# Crear un modelo de SVM tipo kernel polynomial utilizando el paquete caret
# no hace falta tunear el parámetro C "cost" 
svmModelKernelPolynomial <- train(Fase ~.,
                                  data = train_data_smote,
                                  method = "svmPoly",
                                  trControl = trainControl(method = "cv", number = 10),
                                  preProcess = c("center", "scale"),
                                  tuneLength = 5,
                                  prob.model = TRUE) 
svmModelKernelPolynomial

plot(svmModelKernelPolynomial)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
svmModelKernelPolynomial_predictions <- predict(svmModelKernelPolynomial, newdata = test_data)
svmModelKernelPolynomial_predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(svmModelKernelPolynomial_predictions, test_data$Fase)

# SVM kernel
probabilities_svm_kernelpol <- predict(svmModelKernelPolynomial, newdata = test_data, type = "prob")
probabilities_svm_kernelpol

#### ---- Decission Tree ----
# Crear un modelo de DT utilizando el paquete caret
dtModel <- train(Fase ~.,
                 data = train_data_smote,
                 method = "rpart",
                 trControl = trainControl(method = "cv", number = 10),
                 preProcess = c("center", "scale"),
                 tuneLength = 10)
dtModel
plot(dtModel)

fancyRpartPlot(dtModel$finalModel, type=4)


# Evaluar el modelo con el conjunto de prueba
dtModel_predictions <- predict(dtModel, newdata = test_data, type = "raw") # raw = clases
dtModel_predictions


# Evaluar la precisión del modelo utilizando la matriz de confusión
dtModel_confusionMatrix <- confusionMatrix(dtModel_predictions, test_data$Fase)
dtModel_confusionMatrix

# Obtener probabilidades
probabilities_dt <- predict(dtModel, newdata = test_data, type = "prob")
probabilities_dt

f1_score_dtModel<- as.numeric(dtModel_confusionMatrix$byClass['F1'])
f1_score_dtModel
######Random Forest####
library(caret)
install.packages("randomForest")
library(randomForest)

set.seed(123)

rf_model <- train(Fase ~ .,
                  data = train_data_smote,
                  method = "rf",
                  trControl = trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary),
                  metric = "ROC",
                  tuneLength = 5)
rf_model



# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
rf_model_predictions <- predict(rf_model, newdata = test_data)
rf_model_predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
rf_model_confusionMatrix <- confusionMatrix(rf_model_predictions, test_data$Fase)
rf_model_confusionMatrix
# SVM kernel
probabilities_rf_model <- predict(rf_model, newdata = test_data, type = "prob")
probabilities_rf_model

f1_score_rf_model<- as.numeric(rf_model_confusionMatrix$byClass['F1'])
f1_score_rf_model

#####XGBoost ######
install.packages("xgboost")
library(xgboost)

install.packages("xgboost", type = "binary")

set.seed(123)

xgb_model <- train(Fase ~ .,
                   data = train_data_smote,
                   method = "xgbTree",
                   trControl = trainControl(method = "cv",
                                            number = 10,
                                            classProbs = TRUE,
                                            summaryFunction = twoClassSummary),
                   metric = "ROC",
                   tuneLength = 5)



# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
xgb_model_predictions <- predict(xgb_model, newdata = test_data)
xgb_model_predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
xgb_model_confusionMatrix <- confusionMatrix(xgb_model_predictions, test_data$Fase)
xgb_model_confusionMatrix
# SVM kernel
probabilities_xgb_model <- predict(xgb_model, newdata = test_data, type = "prob")
probabilities_xgb_model

f1_score_xgb_model<- as.numeric(xgb_model_confusionMatrix$byClass['F1'])
f1_score_xgb_model



###### Curvas ROC ######
#### ---- Curvas ROC ----
roc_knn <- roc(test_data$Fase, probabilities_knnmodel[,"Temprana"]) # Cambia [,2] según la clase positiva
auc_knn <- auc(roc_knn)
cat("AUC k-NN:", auc_knn, "\n")

roc_svmlinear <- roc(test_data$Fase, probabilities_svm_linear[,"Temprana"]) # Cambia [,2] según la clase positiva
auc_svmlinear <- auc(roc_svmlinear)
cat("AUC svmlinear", auc_svmlinear, "\n")

roc_dt <- roc(test_data$Fase, probabilities_dt[,"Temprana"]) # Cambia [,2] según la clase positiva
auc_dt <- auc(roc_dt)
cat("AUC dt:", auc_dt, "\n")

roc_rf <- roc(test_data$Fase, probabilities_rf_model[,"Temprana"]) # Cambia [,2] según la clase positiva
auc_rf <- auc(roc_rf)
cat("AUC rf:", auc_rf, "\n")

roc_xgb_model <- roc(test_data$Fase, probabilities_xgb_model[,"Temprana"]) # Cambia [,2] según la clase positiva
auc_xgb_model <- auc(roc_xgb_model)
cat("AUC xgb_model:", auc_xgb_model, "\n")

# Curvas ROC con límites de ejes fijos (0 a 1)
plot(roc_knn, col = "blue", main = "Curvas ROC", lwd = 2, 
     legacy.axes = TRUE)

plot(roc_svmlinear, col = "purple", add = TRUE, lwd = 2, 
     legacy.axes = TRUE)
plot(roc_dt, col = "green", add = TRUE, lwd = 2, 
     legacy.axes = TRUE)
plot(roc_rf, col = "orange", add = TRUE, lwd = 2, 
     legacy.axes = TRUE)
plot(roc_xgb_model, col = "red", add = TRUE, lwd = 2, 
     legacy.axes = TRUE)


# Agregar leyenda
knn_legend <- paste("AUC k-NN:", round(auc_knn, 2))
svmlinear_legend <- paste("AUC SVM:", round(auc_svmlinear, 2))
dt_legend <- paste("AUC Decision Tree:", round(auc_dt, 2))
rf_legend <- paste("AUC RF:", round(auc_rf, 2))
xgb_model_legend <- paste("AUC XGBoost:", round(auc_xgb_model, 2))

legend("bottomright", 
       legend = c(knn_legend, svmlinear_legend, dt_legend, rf_legend, xgb_model_legend),
       col = c("blue", "purple", "green", "orange", "red"), 
       lwd = 2)


###### Precision Recall ######


# Vector binario: 1 si Temprana, 0 si Tardia
labels_temprana <- ifelse(test_data$Fase == "Temprana", 1, 0)

# Revisamos por si acaso:
table(labels_temprana)

pr_knn <- pr.curve(scores.class0 = probabilities_knnmodel[,1], weights.class0 = test_data$Fase == "Temprana", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_svm_linear <- pr.curve(scores.class0 = probabilities_svm_linear[,1], weights.class0 = test_data$Fase == "Temprana", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_dt <- pr.curve(scores.class0 = probabilities_dt[,1], weights.class0 = test_data$Fase == "Temprana", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_rf <- pr.curve(scores.class0 = probabilities_rf_model[,1], weights.class0 = test_data$Fase == "Temprana", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_xgb_model <- pr.curve(scores.class0 = probabilities_xgb_model[,1], weights.class0 = test_data$Fase == "Temprana", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)

probabilities_knnmodel
scores.class0 = probabilities_knnmodel[,1]
scores.class0
# Mostrar curva
plot(pr_knn, col = "blue", lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_svm_linear, col = "purple", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_dt, col = "red", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_rf, col = "green", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_xgb_model, col = "orange", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)


# Agregar leyenda
knn_legend <- paste("PR-Curve k-NN:", round(pr_knn$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
svm_linear_legend <- paste("PR-Curve SVM linear:", round(pr_svm_linear$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
dt_legend <- paste("PR-Curve Decission Tree:", round(pr_dt$auc.integral, 2)) 
rf_legend <- paste("PR-Curve random Forest:", round(pr_rf$auc.integral, 2)) 
xgb_model_legend <- paste("PR-Curve XGBoost:", round(pr_xgb_model$auc.integral, 2)) 

legend("bottomright", legend = c(knn_legend, svm_linear_legend, dt_legend, rf_legend, xgb_model_legend ),
       col = c("blue", "purple", "red", "green", "orange"), lwd = 2)

