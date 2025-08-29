library(HIMA)
library(caret)

# combine abiotic and biotic data in a list for HIMA
list_dat=list(cbind(abiotic_dat,Y),bio_dat)
names(list_dat)=c("abiotic","biotic")

list_dat$biotic$Sample=NULL # remove sample column


dummy <- dummyVars("~ .", data = list_dat$abiotic)

list_dat$abiotic=data.frame(predict(dummy, newdata = list_dat$abiotic))

# scale abbiotic data
list_dat$abiotic=cbind(list_dat$abiotic[,13:19],scale(list_dat$abiotic[,-c(13:19)]))

list_dat$abiotic$field_type.conv=NULL
list_dat$abiotic$season.Spring=NULL
list_dat$abiotic$depth.10=NULL

# PCA-mix of abotic variables  
library(PCAmixdata)
X1=list_dat$abiotic[,-4] # take out NDVI of X
X1$field_type.eco=as.factor(X1$field_type.eco)
X1$season.Summer=as.factor(X1$season.Summer)
X1$depth.30=as.factor(X1$depth.30)

names(X1)=c("Field type","Season","Soil depth","P","K","N","pH","EC","Soil temperature yearly",
            "Moisture yearly","Air temperature yearly","Soil temperature prior","Moisture prior",
            "Air temperature prior","Moisture Sampling")

X.quanti <- splitmix(X1)$X.quanti
X.quali <- splitmix(X1)$X.quali

# do the PCA-mix and get eigenvalues and variance for plot
pca<-PCAmix(X.quanti,X.quali,ndim=8,rename.level = T)
PC_rot=PCArot(pca,dim=8)
PC_rot$scores
PC_rot$eig

cat_coords <- PC_rot$levels
cor_quanti=PC_rot$quanti$coord



pc_data <- data.frame(
  Component = rownames(PC_rot$V), 
  Value =rbind(cor_quanti,cat_coords$coord)[, 5]
)


# Create the bar plot with rotated x-axis labels
ggplot(pc_data, aes(x = Component, y = Value)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "PC5 loadings",
    x = "Component",
    y = "Loading"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



eig_df <- data.frame(
  Component = seq_along(PC_rot$eig[,1]),
  Eigenvalue = PC_rot$eig[,1]
)

var_df <- data.frame(
  Component = seq_along(PC_rot$eig[,1]),
  Variance = cumsum(PC_rot$eig[,2])
)

# Plot
p1=ggplot(eig_df, aes(x = Component, y = Eigenvalue)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # cutoff line
  labs(title = "Eigenvalues of PCA mixture",
       x = "Principal Component",
       y = "Eigenvalue") +
  theme_minimal()


p2=ggplot(var_df, aes(x = Component, y = Variance)) +
  geom_point(size = 3) +
  geom_line() +
#  geom_hline(yintercept = sum(PC_rot$eig[1:5,2]), linetype = "dashed", color = "red") +  # cutoff line
  labs(title = "Eigenvalues of PCA mixture",
       x = "Principal Component",
       y = "Variance (%)") +
  theme_minimal()

library(patchwork)

# Combine side by side
(p1 + p2) + plot_annotation(tag_levels = 'A')  # Will label p1 as A, p2 as B


X_PC=data.frame(PC_rot$ind$coord)
names(X_PC)=paste0("dim.",1:8)
X_PC$Y=list_dat$abiotic$Y

# get correlation plot of orginal variables
list_dat$abiotic_dummy=list_dat$abiotic
names(list_dat$abiotic_dummy)=c("Field type","Season","Soil depth","NDVI","P","K","N","pH","EC","Soil temperature yearly",
                                "Moisture yearly","Air temperature yearly","Soil temperature prior","Moisture prior",
                                "Air temperature prior","Moisture Sampling")

cor_mat=cor(list_dat$abiotic_dummy)
corrplot::corrplot(cor_mat)

#### HIMA RUNS ######

# lasso run
e1e1 <- hima(Y ~dim.1+dim.2+dim.3+dim.4+dim.5,#+dim.6+dim.7+dim.8+dim.9+dim.10+dim.11,#season.Summer+field_type.eco
            data.pheno =X_PC,#list_dat$abiotic, #X_PC,
            data.M =decostand(list_dat$biotic,"rclr") ,
#            mediator.type = "compositional",
            efficient = F,
            penalty = "lasso", # Efficient HIMA does not support DBlasso
            scale = F,
            #quantile = TRUE,
            #COV = list_dat$abiotic[, c("prior.soil.moisture","P","N","yearly.moisture","N","P","K","EC","pH")],
            sigcut = 0.8
) # Disabled only for simulation data
summary(e1e1)


p_adj=p.adjust(e1e1$`p-value`, method = "fdr",n=length(e1e1$`p-value`)) # check all Q.values are within FDR (HIMA does not do it by itself)

list_dat$abiotic1=list_dat$abiotic
list_dat$biotic_1=list_dat$biotic
names(list_dat$abiotic1)=paste0("X",1:dim(list_dat$abiotic)[2])
names(list_dat$biotic_1)= paste0("M", 1:dim(list_dat$biotic)[2])

test_list1=names(list_dat$biotic_1[,which(names(list_dat$biotic) %in%  e1e1$ID[e1e1$`p_adj`<0.1])])


# MCP run
e1e1 <- hima(Y ~dim.1+dim.2+dim.3+dim.4+dim.5,#+dim.6+dim.7+dim.8+dim.9+dim.10+dim.11,#season.Summer+field_type.eco
             data.pheno =X_PC,#list_dat$abiotic, #X_PC,
             data.M =decostand(list_dat$biotic,"rclr") ,
             #            mediator.type = "compositional",
             efficient = F,
             penalty = "MCP", # Efficient HIMA does not support DBlasso
             scale = F,
             #quantile = TRUE,
             #COV = list_dat$abiotic[, c("prior.soil.moisture","P","N","yearly.moisture","N","P","K","EC","pH")],
             sigcut = 0.8
) # Disabled only for simulation data
summary(e1e1)

list_dat$abiotic1=list_dat$abiotic
list_dat$biotic_1=list_dat$biotic
names(list_dat$abiotic1)=paste0("X",1:dim(list_dat$abiotic)[2])
names(list_dat$biotic_1)= paste0("M", 1:dim(list_dat$biotic)[2])

p_adj=p.adjust(e1e1$`p-value`, method = "fdr",n=length(e1e1$`p-value`)) # check all Q.values are within FDR (HIMA does not do it by itself)

test_list2=names(list_dat$biotic_1[,which(names(list_dat$biotic) %in% e1e1$ID[p_adj<0.1])])

# DBlasso run
e1e1 <- hima(Y ~dim.1+dim.2+dim.3+dim.4+dim.5,#+dim.6+dim.7+dim.8+dim.9+dim.10+dim.11,#season.Summer+field_type.eco
             data.pheno =X_PC,#list_dat$abiotic, #X_PC,
             data.M =list_dat$biotic,#decostand(list_dat$biotic,"rclr") ,
             #            mediator.type = "compositional",
             efficient = F,
             penalty = "DBlasso", # Efficient HIMA does not support DBlasso
             scale = F,
             #quantile = TRUE,
             #COV = list_dat$abiotic[, c("prior.soil.moisture","P","N","yearly.moisture","N","P","K","EC","pH")],
             sigcut = 0.8
) # Disabled only for simulation data
summary(e1e1)

p_adj=p.adjust(e1e1$`p-value`, method = "fdr",n=length(e1e1$`p-value`)) # check all Q.values are within FDR (HIMA does not do it by itself)

list_dat$abiotic1=list_dat$abiotic
list_dat$biotic_1=list_dat$biotic
names(list_dat$abiotic1)=paste0("X",1:dim(list_dat$abiotic)[2])
names(list_dat$biotic_1)= paste0("M", 1:dim(list_dat$biotic)[2])

test_list3=names(list_dat$biotic_1[,which(names(list_dat$biotic) %in% e1e1$ID[p_adj<0.1])])

# remove dublicates
test_list_combo=unique(c(test_list1,test_list2,test_list3))

# extract names for paper
cbind(names(list_dat$biotic[,which(names(list_dat$biotic_1) %in% test_list_combo)]),
      names(list_dat$biotic_1[,which(names(list_dat$biotic_1) %in% test_list_combo)]))

# take the data for SEM in the next script
data_sem=data.frame(cbind(Y=(list_dat$abiotic$Y),X_PC,decostand(list_dat$biotic_1,"rclr")))





