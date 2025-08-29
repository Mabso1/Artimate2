# read in bio-informatics data
library(readr)
library(readxl)
library(Metrics)
library(tidyr)
library(ggpattern)
library(ggpattern)
library(patchwork)

dat = read_excel("~/90-1117990965/Report/Result/01_OTU_Taxa/OTU/otu_taxa_table.xls")

# remove enteries with out<2
dat[dat < 2] = NA

####### genus level data #######
# keep rows which has genus level
genus_idx=which(grepl("g__",dat$taxonomy)==TRUE)

dat_genus=dat[genus_idx,]

genus_taxa=dat_genus$taxonomy # save taxa

# remove counter column and taxa ID column 
dat_genus=dat_genus[,-c(1,218)]


####### species level data ######
# keep rows which has species level
species_idx=which(grepl("s__",dat$taxonomy)==TRUE)

dat_species=dat[species_idx,]

species_taxa=dat_species$taxonomy # save taxa

# remove counter column and taxa ID column
dat_species=dat_species[,-c(1,218)]

# merge duplicated species names OTU counts within each field



OTU_dublicate_merge <- function(OTU_dat,OTU_taxa,all_fields) {
  

  if(all_fields==TRUE){
    # unique taxa
    species_taxa_unique=unique(OTU_taxa)
    
    N=length(species_taxa_unique)
    K=dim(OTU_dat)[2]
    
    
    dat_species_reduced=data.frame(matrix(,nrow=N,ncol=K))
    names(dat_species_reduced)=names(OTU_taxa)
    
    for (i in 1:N ){
      idx_sum=which(species_taxa_unique[i]==OTU_taxa)
      
      dat_species_reduced[i,]=colSums(OTU_dat[idx_sum,],na.rm = T)  
      
    }
    
    # remove rows with sums =0 (taxa not present in any samples)
    idx_row_rm=which(rowSums(dat_species_reduced)==0)
    dat_species_reduced=dat_species_reduced[-idx_row_rm,]
    species_taxa_reduced=species_taxa_unique[-idx_row_rm]
    
    out=list(dat_species_reduced, species_taxa_reduced)
    
  }
  
  else{
    # find samples within each field
    
    # unique taxa
    species_taxa_unique=unique(OTU_taxa)
    
    # unique fields
    idx_con_sp=grepl("Con-Sp",names(OTU_dat))
    idx_con_su=grepl("Con-Su",names(OTU_dat))
    idx_eco_sp=grepl("Eco-Sp",names(OTU_dat))
    idx_eco_su=grepl("Eco-Su",names(OTU_dat))
    
    field_idx=list(idx_con_sp,idx_con_su,idx_eco_sp,idx_eco_su)
    field_OTU=list()
    J=length(field_idx)
    
    # get OTU and unique taxa from each field
    
    for (j in 1:J){
      N=length(species_taxa_unique)
      K=dim(OTU_dat[,field_idx[[j]]])[2]
      
      
      dat_species_reduced=data.frame(matrix(,nrow=N,ncol=K))
      names(dat_species_reduced)=names(OTU_taxa)
      
      for (i in 1:N ){
        idx_sum=which(species_taxa_unique[i]==OTU_taxa)
        
        dat_species_reduced[i,]=colSums(OTU_dat[idx_sum,field_idx[[j]]],na.rm = T)  
        
      }
      
      # remove rows with sums =0 (taxa not present in any samples)
      idx_row_rm=which(rowSums(dat_species_reduced)==0)
      dat_species_reduced=dat_species_reduced[-idx_row_rm,]
      species_taxa_reduced=species_taxa_unique[-idx_row_rm]
      
      field_OTU[[j]]=list(dat_species_reduced, species_taxa_reduced)
      
      out=field_OTU # 4 list of 2 elements 
    }
    
    }
    
    
  

  return(out)
}


# change dat_genus and genus_taxa to dat_species and species_taxa if working with species level
all_fields_OTU=OTU_dublicate_merge(dat_genus,genus_taxa,all_fields=TRUE)


####### correct the format to work with the rest of the data ######

field_dat = read_csv("C:/Users/dasth/OneDrive/Desktop/field trials/Abiotic_data.csv")

# make sure IDs have the same names in both bio and field data - choose the bio data notation

field_dat$ID=gsub("y_sp","-Eco-Sp",field_dat$ID)
field_dat$ID=gsub("x_sp","-Con-Sp",field_dat$ID)
field_dat$ID=gsub("y_su","-Eco-Su",field_dat$ID)
field_dat$ID=gsub("x_su","-Con-Su",field_dat$ID)

new_names=c()
N=length(field_dat$ID)

for (i in 1:N ){
  
  num_temp=as.character(as.numeric(strsplit(field_dat$ID,"-")[[i]][1])-1)
  dummy_string=strsplit(field_dat$ID,"-")[[i]]
  dummy_string[1]=num_temp
  dummy_string=paste(dummy_string, collapse="-")
  
  new_names[i]=dummy_string

}

correction_idx=which(grepl("9-",new_names)==TRUE)
new_names[c(correction_idx[1],correction_idx[1]+c(1,2))]=c("10-Eco-Sp", "11-Eco-Sp", "12-Eco-Sp")
new_names[c(correction_idx[2],correction_idx[2]+c(1,2))]=c("10-Con-Sp", "11-Con-Sp", "12-Con-Sp")
new_names[c(correction_idx[3],correction_idx[3]+c(1,2))]=c("10-Eco-Su",  "11-Eco-Su", "12-Eco-Su")
new_names[c(correction_idx[4],correction_idx[4]+c(1,2))]=c("10-Con-Su" , "11-Con-Su", "12-Con-Su")

# give new name
field_dat$ID=new_names

###### merge with species data ######
dat_species_reduced=all_fields_OTU[[1]]
species_taxa_reduced=all_fields_OTU[[2]]

dat_species_t=data.frame(t(dat_species_reduced)) # sample in rows and taxa in columns
names(dat_species_t)=species_taxa_reduced

# remove names which contain unidentified
idx_rm=which(grepl("unidentified",names(dat_species_t))==TRUE)
dat_species_t=dat_species_t[,-idx_rm]

# trim names of species data keeping only species/genus part
trim_names <- function(OTU_transposed,OTU_level) {
  species_names=c()
  N=length(names(dat_species_t))
  
  if(OTU_level=="Species"){
    for (i in 1:N ){
      species_names[i]=strsplit(names(dat_species_t)[i],"s__")[[1]][2]
    }
  }
  
  else if (OTU_level=="Genus"){
    for (i in 1:N ){
      temp_string=strsplit(names(dat_species_t)[i],"g__")[[1]][2]
      species_names[i]=strsplit(temp_string,";")[[1]][1]
      
    }
  }
 
  return(species_names)
  
}

species_names=trim_names(dat_species_t,"Genus")

names(dat_species_t)=species_names

######## if working with genus data press this:

# Step 1: Identify unique column names
unique_columns <- unique(names(dat_species_t))

# Step 2: Create an empty data frame to store the result
df_sum <- data.frame(matrix(ncol = length(unique_columns), nrow = nrow(dat_species_t)))
names(df_sum) <- unique_columns

# Step 3: Loop over unique column names and sum identical columns
for (col in unique_columns) {
  # Find all column indices with the same name
  cols_to_sum <- dat_species_t[, names(dat_species_t) == col, drop = FALSE]
  
  # Sum across the rows for these columns
  df_sum[[col]] <- rowSums(cols_to_sum, na.rm = TRUE)
}

# View the result
dat_species_t=df_sum
########## do not press the above if wokring with species data, go from here!


# get rid of less abundant species, use IQR cut off to see whats left
data_matrix=dat_species_t 

# IQR filteration 
library(boot)
library(vegan)
library("factoextra")

# IQR funcion 
IQR_boot <- function(data_matrix,N_runs,conf_int) {
  
  Function = function(input, index){
    Input = input[index]
    Stat = IQR(Input,na.rm=TRUE)
    return(Stat)
  }
  
  
  
  data_matrix=decostand(data_matrix,"rclr") # get scaled relative abudances
  
  # running this might take some time
  boot_list=list()
  for (i in 1:length(data_matrix)){
    Boot = boot(data_matrix[,i], Function, R=N_runs)
    mean_boot=median(Boot$t)
    boot_CI=boot.ci(Boot, conf = conf_int, type = "norm")
    lower_CI=boot_CI$normal[c(2)]
    upper_CI=boot_CI$normal[c(3)]
    boot_CI=boot.ci(Boot, conf = conf_int, type = "basic")
    lower_CI_p=boot_CI$basic[c(4)]
    upper_CI_p=boot_CI$basic[c(5)]
    
    boot_list[[i]]=c(mean_boot,lower_CI,upper_CI, lower_CI_p,upper_CI_p)
  }
  
  
  # cluster 1 boot results (rerun to above above when switching to new cluster)
  boot_matrix=do.call("rbind",boot_list)
  colnames(boot_matrix)=c("mean","CI_lower","CI_upper","CI_lower_p","CI_upper_p")
  rownames(boot_matrix)=colnames(data_matrix)
  boot_matrix=as.data.frame(boot_matrix)
  
  # identify which CIs contain 0
  zero_list=c()
  for (i in 1:length(data_matrix)){
    zero_list[i] <- 0 >= boot_matrix$CI_lower[i] & 0 <= boot_matrix$CI_upper[i]
  }
  
  
  idx0_c1_plus=which(zero_list=="FALSE")
  
  
  boot_matrix$genus_no=seq(1:length(boot_matrix[,1]))
  boot_matrix$Zero=rep(NA,length(boot_matrix[,1]))
  boot_matrix$Zero[idx0_c1_plus]="No"
  boot_matrix$Zero[-idx0_c1_plus]="Yes"
  boot_matrix$Zero=as.factor(boot_matrix$Zero)
  
  return(boot_matrix)
  
}

# IQR filt results
boot_matrix=IQR_boot(data_matrix,2000,0.95)

# extract taxa with overlaps and remove these columns
idx_rm=which(boot_matrix$Zero=="Yes")
dat_species_t_reduced=dat_species_t[,-idx_rm]


# update taxa based on IQR filteration
dat_species_t=dat_species_t_reduced

# plot the IQR filtering results
boot_matrix$ZeroInCI <- ifelse(boot_matrix$CI_lower <= 0 & boot_matrix$CI_upper >= 0, "Yes", "No")
boot_matrix$Component <- seq_len(nrow(boot_matrix))  # For x-axis
boot_matrix$Variable <- names(data_matrix)

ggplot(boot_matrix, aes(x = Variable, y = mean, color = ZeroInCI)) +
  geom_point(size = 3) +  # Mean point
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +  # CI bars
  geom_hline(yintercept = 0, linetype = "dashed") +  # Zero line
  scale_color_manual(values = c("Yes" = "red", "No" = "blue")) +
  scale_x_discrete(expand = expansion(mult = c(0.00005, 0.00005)))+
  labs(
    title = "Bootstrapped IQR and Confidence Intervals",
    x = "Species",
    y = "IQR estimate",
    color = "Zero included in CI"
  ) +
  theme_minimal() +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(10, 5, 5, 5))


ggplot(boot_matrix, aes(y = Variable, x = mean, color = ZeroInCI)) +
  geom_point(size = 3) +  # Mean point
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +  # Horizontal CI bars
  geom_vline(xintercept = 0, linetype = "dashed") +  # Zero line
  scale_color_manual(values = c("Yes" = "red", "No" = "blue")) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Bootstrapped IQR and Confidence Intervals",
    y = "Genus",
    x = "IQR Estimate",
    color = "Zero included in CI"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.margin = margin(10, 5, 5, 5)
  )


# combine the abotic and filtered biotic data
merged_data=cbind(dat_species_t,field_dat$`field_type`,field_dat$season,field_dat$depth)

K=dim(dat_species_t)[2]

names(merged_data)[(K+1):(K+3)]=c("field type","season","depth")



##### analysis of data ####
merged_data=cbind(dat_species_t,field_dat$P,field_dat$K,field_dat$N,field_dat$pH,field_dat$`field_type`,field_dat$season,field_dat$depth, field_dat$EC, field_dat$replicate,
                  field_dat$soil_temp_yearly,field_dat$soil_water_yearly,field_dat$air_temp_yearly,
                  field_dat$prev_soil_temp, field_dat$prev_soil_water,field_dat$prev_air_temp, field_dat$current_water,field_dat$NDVI_median)

K=dim(dat_species_t)[2]

field_names=c("P","K","N","pH","field_type","season","depth","EC","replicate","yearly soil temp","yearly moisture",
  "yearly air temp","prior soil temp","prior soil moisture", "prior air temp", "sample moisture","NDVI")
N=length(field_names)

names(merged_data)[(K+1):(K+N)]=field_names

# correctly encode the factors
merged_data$depth=as.factor(merged_data$depth)
merged_data$`field_type`=as.factor(merged_data$`field_type`)

merged_data$season=as.factor(merged_data$season)

# aggregate data by replicate
library(dplyr)

agg_df <- merged_data %>%
  group_by(replicate) %>%
  summarise(
    across(where(is.numeric), median),  # Aggregate numeric variables
    field_type = names(which.max(table(field_type))),
    season = names(which.max(table(season))),# Select majority class
    depth = names(which.max(table(depth)))
  )


merged_data=agg_df

merged_data=merged_data[, !colnames(merged_data) %in% "replicate"]

# correctly encode the factors
merged_data$depth=as.factor(merged_data$depth)
merged_data$`field_type`=as.factor(merged_data$`field_type`)
merged_data$season=as.factor(merged_data$season)

merged_data=data.frame(merged_data)

# select only non bio data
k=dim(merged_data)[2]
k1=dim(dat_species_t_reduced)[2]
merged_data_nobio=merged_data


levels(merged_data_nobio$field_type)=c("Conventional","Ecological") # for plotting
levels(merged_data_nobio$depth)=c("10 cm","30 cm")

# soil nurtients levels from spring to summer
df_nutrients <- data.frame(
  Nitrogen = merged_data$N,
  Phosphate = merged_data$P,
  Potassium = merged_data$K,
  season = merged_data$season,
  field_type = merged_data$field_type
) %>%
  pivot_longer(cols = c(Nitrogen, Phosphate,  Potassium), names_to = "Nutrient", values_to = "Value")

# Function to create boxplot by nutrient
plot_nutrient <- function(nutrient_name) {
  y_label <- ifelse(nutrient_name == "Nitrogen", 
                    "Concentration (g/kg)", 
                    "Concentration (mg/kg)")
  
  ggplot(subset(df_nutrients, Nutrient == nutrient_name), aes(x = field_type, y = Value, fill = season)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
    labs(
      title = paste("Concentration of", nutrient_name),
      x = "Field Type",
      y = y_label,
      fill = "Season"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

pN <- plot_nutrient("Nitrogen")
pP <- plot_nutrient("Phosphate")
pK <- plot_nutrient("Potassium")

combined_plot <- pN / pP / pK  +  plot_annotation(tag_levels = "A")# stack vertically

# Display combined plot
print(combined_plot)


# Define colors and patterns
library(ggplot2)


season_patterns <- c("Spring" = "stripe", "Summer" = "crosshatch")
field_colors <- c("Ecological" = "blue", "Conventional" = "orange")

ggplot(merged_data_nobio, aes(x = NDVI,
                              pattern = season,
                              pattern_fill = field_type,
                              color = field_type)) +   
  geom_density_pattern(
    aes(fill = "white"),       
    alpha = 0.8,
    position = "identity",
    color = "black",
    pattern_angle = 45,
    pattern_density = 0.3,
    pattern_spacing = 0.02,
    pattern_key_scale_factor = 0.6,
    size = 0.1                    
  ) +
  scale_pattern_manual(values = season_patterns) +
  scale_pattern_fill_manual(values = field_colors) +
  scale_color_manual(values = field_colors) +  
  scale_fill_identity() +
  guides(
    pattern_fill = guide_legend(
      override.aes = list(pattern = "none", fill = c("blue", "orange"))
    ),
    color = "none"  # hides border color legend to avoid duplication
  ) +
  xlim(0, 0.8) +
  labs(
    title = "NDVI Density by Season and Field Type",
    x = "NDVI",
    y = "Density",
    pattern = "Season",
    pattern_fill = "Field Type",
    color = "Field Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )







ggplot(merged_data_nobio, aes(x = season, y = NDVI, color = field_type)) +
  geom_boxplot() +  # Boxplot to show NDVI distribution
  facet_wrap(~ depth) +  # Separate plots for each depth level
  labs(title = "NDVI Distribution across Depth and Field Type",
       x = "Season",
       y = "NDVI",
       color = "Field Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

ggplot(merged_data_nobio, aes(x = season, y = NDVI, color = season)) +
  geom_boxplot() +  # Boxplot to show NDVI distribution
  facet_wrap(~ depth) +  # Separate plots for each depth level
  labs(title = "NDVI Distribution across Depth and Season",
       x = "Season",
       y = "NDVI",
       color = "Field Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


# encode to biotic and biotic variables for the next step
bio_dat=merged_data[,1:k1]
abiotic_dat=merged_data[,(k1+1):k]


Y=abiotic_dat$NDVI
abiotic_dat$NDVI=NULL
abiotic_dat$field.type=NULL

# abudance plot
library(tidyr)
library(tidyverse)

# dummy set for plotting
bio_dat1=bio_dat

names(bio_dat1) <- paste0(names(bio_dat1), " (M", seq_along(bio_dat1), ")")

bio_dat1$Sample <- 1:nrow(bio_dat1)

# Convert to long format
bio_long <- decostand(bio_dat1,"total") %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Genus",
    values_to = "Abundance"
  )

bio_long <-decostand(bio_dat1,"total") %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Genus",
    values_to = "Abundance"
  ) %>%
  group_by(Genus) %>%
  mutate(TotalAbundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Genus = fct_reorder(Genus, TotalAbundance, .desc = F)) %>%
  select(-TotalAbundance)


p=ggplot(bio_long, aes(y = Genus, x =Abundance)) +
  geom_bar(stat = "identity") +
  labs(x = "Frequency (%)", y = "Genus", fill = "Crop type") +
  ggtitle("Relative abudance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(size = 16),
        axis.text=element_text(size=12),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


p

