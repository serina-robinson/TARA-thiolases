# Install packages
install.packages("bio3d", dependencies = TRUE)
pacman::p_load("tidyverse", "bio3d", "data.table")

# Set working directory
setwd("C:/Users/tabie/OneDrive/Documents/GitHub/TARA-thiolases/")

# Read in the summary info about the homology models
summ <- fread("data/homology_models/50_OleA_troy/orig_names/summaryinfo_extended", sep = "|", fill = T, data.table = F, header = T) %>%
  janitor::clean_names() %>%
  dplyr::filter(!grepl("Resolution", resolution)) %>%
  dplyr::group_by(number_job_description) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(org = paste0(word(number_job_description, sep = "_1_", 2))) %>%
  dplyr::mutate(pdb_id = gsub("PDBTitle: ", "", v12)) %>%
  dplyr::mutate(pdb_id = gsub("crystal structure of ", "", pdb_id)) %>%
  dplyr::mutate(filnams = paste0(job_id, ".final.pdb"))
summ$pdb_id
# Read in the old homology models
fils <- data.frame(list.files("data/homology_models/50_OleA_troy/orig_names/", pattern = ".pdb"), stringsAsFactors = F)
fils
colnames(fils) <- "filnams"
?list.files
dim(summ)
comb <- summ %>%
  dplyr::left_join(fils, ., by = "filnams") %>%
  dplyr::mutate(new_filnam = paste0(number_job_description, ".pdb"))# %>%
  dplyr::filter(pdb_id == "3-oxoacyl-(acyl carrier protein)2 synthase iii, fabh (xoo4209) from xanthomonas oryzae pv.3 oryzae kacc10331")
dim(comb)
1:nrow(comb)
for(i in 1:nrow(comb)) {
  print(i)
  pdb_temp <- read.pdb(paste0("data/homology_models/50_OleA_troy/orig_names/", comb$filnams[i]))
  write.pdb(pdb_temp, file = paste0("data/homology_models/50_OleA_troy/new_names/", comb$new_filnam[i]))
}

# Read in the new homology models
fifty <- list.files("data/homology_models/50_OleA_troy/new_names/")
fifty


# Ok now let's just start with one model

# You should set this to your own working directory
pdb1 <- read.pdb("data/homology_models/50_OleA_troy/new_names/GCDEAADJ_01881_Planctomycetota_Pirellulaceae_UBA721_marine.pdb")
attributes(pdb1)
pdb1$calpha
# Now read in the 'template'
pdb2 <- read.pdb("data/homology_models/4KU5_A.pdb")

# You can also pull PDB files directly from online
pdb3 <- read.pdb("4KU5") # You can also access online PDB records


# Make a B-factor plot (NOTE - this will NOT work with homology models)
plot.bio3d(pdb2$atom$b[pdb2$calpha], sse=pdb2, typ="l", ylab="B-factor")

# Align your two sequences
# perform iterative alignment
aln <- seqaln(pdb1, pdb2, fit = TRUE)
aln$rmsd
?seqaln
# Challenge 1. Read in and rename your 50 homology models from Phyre2 using the scripts above

# Challenge 2. Align one of your models with the Xanthomonas campestris structure (4KU5) A chain.

# Challenge 3. Write the alignment to file. See bio3d tutorial for help on how to do this.
# http://thegrantlab.org/bio3d/tutorials/

# Challenge 4. Visualize the alignment in PyMOL. Alternatively, you can align them
# in PyMol as well using the 'super' command. See PyMOL tutorials for help.

# Super challenge 5. Now try to align all of them each individually against 4KU5_A.
# Hint you will need to use a for loop!

for(i in 1:nrow(comb)) {
  print(i)
  pdb_temp <- read.pdb(paste0("data/homology_models/50_OleA_troy/new_names/", comb$filnams[i]))
  
  write.pdb(pdb_temp, file = paste0("data/homology_models/50_OleA_troy/new_names_aligned/", comb$new_filnam[i]))
}


