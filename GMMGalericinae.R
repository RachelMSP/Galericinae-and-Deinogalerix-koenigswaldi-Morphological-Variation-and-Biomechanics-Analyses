#GPA, PCA, PCA PLOTS + MESH DEFORMATIONS, ANOVAS, CAC, CAC PLOT + MESH DEFORMATIONS >> FOR EXTANT GYMNURE DATASET ONLY

rm(list=ls()) 
setwd("/Users/raquelcardenas/Desktop/prog/gmmgalericinae")

#install.packages("Matrix")
#install.packages("RRPP")
#install.packages("rgl")
#install.packages("geomorph")
#install.packages("Morpho")
#install.packages("Rvcg")


library(Matrix)
library(RRPP)
library(rgl)
library(geomorph)
library(Morpho)
library(Rvcg)

# SECTION 1- LOADING AND ASSIGNING COLORS AND NAMES
# Loading of the specimen files 
specimen_files <- list.files(pattern = "*.csv", full.names = TRUE)

specimen_data_list <- list()
specimen_names <- c()
species_codes <- c()
temp_dfs <- list()  
specimen_data_list = specimen_data_list[-46]

for (file in specimen_files) {
  data <- read.csv(file, header = FALSE, skip = 2)
  landmark_matrix <- data[, 4:6]
  
  # Storing 3D landmark matrix
  specimen_data_list[[length(specimen_data_list) + 1]] <- landmark_matrix
  
  # Extracting and storing the file names
  filename <- basename(file)
  specimen_names <- c(specimen_names, filename)
  
  # Extracting characters 5 and 6 from the filenames (These refer to the 'species code', meaning the 2 letter abbreviation given to each species) 
  species_code <- substr(filename, 5, 6)
  species_codes <- c(species_codes, species_code)
  #species_codes = species_codes[-46]
  
#  Long-format dataframe (for inspection)
  temp_dfs[[length(temp_dfs) + 1]] <- data.frame(
    filename = filename,
    x = landmark_matrix[, 1],
    y = landmark_matrix[, 2],
    z = landmark_matrix[, 3]
  )
}

# Assigning the full species names to the extracted species codes
species_factor <- factor(species_codes,
                         levels = c("Eg", "Hs", "Ns", "Om", "Pt"),    
                         labels = c("Echinosorex gymnura", "Hylomys suillus", "Neotetracus sinensis", "Otohylomys megalotis", "Podogymnura truei"))  

# Converting list to array
landmark_array <- array(unlist(specimen_data_list),
                        dim = c(22, 3, length(specimen_data_list))) #<- 22 landmarks, 3 columns (x ay and z)
#landmark_array = landmark_array[,,-46]

# Adding color factor (for plotting visualization)
species_colors <- as.factor(species_factor)
levels(species_colors) <- c("coral1", "goldenrod1", "chartreuse1", "blueviolet", "skyblue")  # Echinosorex= coral, Hylomys = yellow, Neotetracus = green, Otohylomys = purple. Podogymnura = blue


# For inspecting
dim(landmark_array)[3]       
length(specimen_names)   

# SECTION 2 - GPA RESULTS
gpa_results <- gpagen(landmark_array)
#deleteGPAM <-gpagen(landmark_array)

# SECTION 3 - MESH DEFORMATIONS ALONG PCA
ProcSym_MT = procSym(landmark_array)
which.min(ProcSym_MT$rho)

# Finding speciemen closest to the mean
min_rho_index <- which.min(ProcSym_MT$rho)
closest_landmarks <- landmark_array[, , min_rho_index]
cat(sprintf("Closest specimen is: %s\n", specimen_names[min_rho_index]))

sur.clo <- file2mesh("~/Desktop/prog/meshes_gymnures/M034HsML.stl")


# Array of the landmarks on the consensus model
mean_shape <- ProcSym_MT$mshape
mean_mesh <- tps3d(sur.clo, closest_landmarks, mean_shape)

vcgPlyWrite(mean_mesh, filename = "~/Desktop/prog/meshes_gymnures/gymnure_mean.ply")

# Mean shape
pcscores <- ProcSym_MT$PCscores
pcs <- ProcSym_MT$PCs
output_dir <- "~/Desktop/3D_Warped_Models/"
dir.create(output_dir, showWarnings = FALSE)

# Deformation of the mesh along the PC axes:
deform_mesh_along_PC <- function(pc_number, mesh, mean_shape, pcscores, pcs, output_dir) {
  # Deforming toward minimum along PC
  min_score <- min(pcscores[, pc_number])
  min_shape <- showPC(min_score, pcs[, pc_number], mean_shape)
  min_mesh <- tps3d(mesh, mean_shape, min_shape)
  
  max_score <- max(pcscores[, pc_number])
  max_shape <- showPC(max_score, pcs[, pc_number], mean_shape)
  max_mesh <- tps3d(mesh, mean_shape, max_shape)
  
  # Exporting
  min_filename <- sprintf("PC%d_min_mesh.ply", pc_number)
  max_filename <- sprintf("PC%d_max_mesh.ply", pc_number)
  
  vcgPlyWrite(min_mesh, file.path(output_dir, min_filename))
  vcgPlyWrite(max_mesh, file.path(output_dir, max_filename))
  
  cat(sprintf("Meshes for PC%d saved to %s\n", pc_number, output_dir))
}

# Applying for PC1 and PC2 
deform_mesh_along_PC(pc_number = 1, mesh = mean_mesh, mean_shape = mean_shape,
                     pcscores = pcscores, pcs = pcs, output_dir = output_dir)

deform_mesh_along_PC(pc_number = 2, mesh = mean_mesh, mean_shape = mean_shape,
                     pcscores = pcscores, pcs = pcs, output_dir = output_dir)

# SECTION 4 - PCA RESULTS AND PLOTTING
pca_results <- procSym(landmark_array)
pca_results$rotated

# Data frame
PC_data <- as.data.frame(pca_results$PCscores)

pcs <- pca_results$PCscores

# Computing eigenvalues (variance of each PC axis)
eigenvalues <- apply(pcs, 2, var)

# Computing percentage variance explained
explained_var <- eigenvalues / sum(eigenvalues) * 100

# Extracting PC1 and PC2 variance
pc1_var <- round(explained_var[1], 1)
pc2_var <- round(explained_var[2], 1)

# Plotting PCA Plot
xlim_range <- c(min(pcs[, 1]) - 0.020, max(pcs[, 1]) + 0.030)
ylim_range <- c(min(pcs[, 2]) - 0.03, max(pcs[, 2]) + 0.03)
plot(pcs[, 1], pcs[, 2],
     bg = as.character(species_colors),
     col = "black",
     pch = 21,
     #cex = log(gpa_results$Csize)*.5,
     xlab = paste0("PC1 (", pc1_var, "%)"),
     ylab = paste0("PC2 (", pc2_var, "%)"),
     main = "PCA - Galericinae (extant)",
     xlim = xlim_range,
     ylim = ylim_range)

#dev.off()


##Legend
legend("topright", legend = levels(species_factor) , col = levels(species_colors), cex = 0.80, pch = 16, title = "Galericinae Species")

# Text for visualizing individual specimens among morphoscape
text (
  PC_data[,1],
  PC_data[,2],
  labels = (specimen_names),
  pos = 3,
  col = "black",
  cex = 0.6
)


#specimen_names <- basename(specimen_files)


# SECTION 5 - Procrustes ANOVAs
# Geomorph data frame
gdf <- geomorph.data.frame(
  coords = gpa_results$coords,
  Csize = gpa_results$Csize,
  Species = species_factor
)

#CAC=CAC()

# Unique trajectory model  ANOVAs 
anova_with_log <- procD.lm(coords ~ log(Csize) * Species, data = gdf)
summary(anova_with_log)

# Common trajectory model  ANOVAs 
anova_with_log_common <- procD.lm(coords ~ log(Csize) + Species, data = gdf)
summary(anova_with_log_common)
#plot(anova_with_log)

# ANOVA between both models
anova(anova_with_log_common,anova_with_log)
summary(test)


# SECTION 6 - CAC analysis and CAC plot

cac_result <- CAC(gdf$coords, log(gdf$Csize))
scores <- cac_result$CACscores

str(cac_result$CACscores)

#Plotting
plot(log(gdf$Csize), cac_result$CACscores,
     pch = 21, col = "black", bg = as.character(species_colors),
     xlab = "log(Csize)", ylab = "Common Allometric Component",
     main = "Shape Change Along Common Allometry")

# Legend for CAC plot
legend("bottomright", legend = levels(species_factor), col = "black", pt.bg = levels(species_colors),
 pch = 21, cex = 0.7, text.font = 3, title = "Galericinae Species")

# Text to see the specimen corresponding to dots in CAC
text(log(gdf$Csize), cac_result$CACscores,
     labels = specimen_names, pos = 3, cex = 0.7)


# SECTION 7 - MESH DEFORMATIONS ALONG CAC AXIS

# Inputs
landmarks <- gdf$coords               # 3D landmark array
specimen_names <- dimnames(landmarks)[[3]]
log_csize <- log(gdf$Csize)
cac_scores <- cac_result$CACscores[, 1]  # first CAC axis
cac_vector <- cac_result$CAC[, 1]        # shape change vector for CAC1

# Finding specimen closest to the  mean
mean_log_csize <- mean(log_csize)
mean_cac <- mean(cac_scores)

combined_distance <- (log_csize - mean_log_csize)^2 + (cac_scores - mean_cac)^2
ref_index <- which.min(combined_distance)

cat(sprintf("Closest to CAC+size mean: %s\n", specimen_files[ref_index]))

ref_lm <- landmarks[, , ref_index]  # reference specimen landmarks

# Loading specimen that’s closest to the mean
ref_mesh <- file2mesh("~/Desktop/prog/meshes_gymnures/M053PtML.stl")

# CAC vector into a 3D matrix (3 × p)
k <- dim(landmarks)[1]  # 3D
p <- dim(landmarks)[2]  # number of landmarks
cac_vector_mat <- matrix(cac_vector, nrow = 1)
cac_array <- arrayspecs(cac_vector_mat, k, p)
cac_matrix <- cac_array[, , 1]

# How far to warp along CAC axis
delta_min <- min(log_csize) - mean_log_csize  # move toward small
delta_max <- max(log_csize) - mean_log_csize  # move toward large

# Warping landmarks
lm_small <- ref_lm + delta_min * cac_matrix  # negative delta
lm_large <- ref_lm + delta_max * cac_matrix  # positive delta

# Warping using TPS
mesh_small <- tps3d(ref_mesh, ref_lm, lm_small)
mesh_large <- tps3d(ref_mesh, ref_lm, lm_large)
#shade3d(mesh_small,col = "darkgrey", add = T) #to view
#shade3d(mesh_large,col = "darkgrey", add = T) #

# Exporting
output_dir <- "~/Desktop/ONLYG3D_Warped_Models_CAC"
dir.create(output_dir, showWarnings = FALSE)

vcgPlyWrite(mesh_small, file.path(output_dir, "CAC_min_mesh.ply"))
vcgPlyWrite(mesh_large, file.path(output_dir, "CAC_max_mesh.ply"))

cat("CAC-deformed meshes saved to:", output_dir, "\n")
