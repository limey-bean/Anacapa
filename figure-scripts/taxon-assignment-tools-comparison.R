# Script to make comparison figure to compare CRUX and a few other 16s databases
# Data stored in Anacapa/data-for-comparisons
# First we make the 16s comparison
library(tidyverse)
library(reshape2)
library(superheat)
mock <- read_delim("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/16S/16S_actual_taxonomy.txt", delim ="\t")
# Split the taxonomy by rank
mock <- cbind(mock,colsplit(mock$`Actual Taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))

file_names <- list.files("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/16S/")
# Get rid of the mock...
file_names <- file_names[!(str_detect(file_names, "actual"))]
all_files <- lapply(file_names, function(file) read_delim(file.path("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/16S", file), delim = "\t"))
names(all_files) <- file_names

all_files <- lapply(all_files, function(x) 
  x %>% mutate(`assigned taxonomy` = str_replace(`assigned taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `assigned taxonomy` = str_replace_all(`assigned taxonomy`, "[a-z]__", ""),     # Greengenes taxonomy has p__/c__ for each level
               `assigned taxonomy` = str_replace_all(`assigned taxonomy`, "D_[1-9]__", ""),   # Silva has D_1__ etc.
               `assigned taxonomy` = str_replace_all(`assigned taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `assigned taxonomy` = str_replace_all(`assigned taxonomy`, ";$", ";NA"),       # The next few lines add NAs between ;;
               `assigned taxonomy` = str_replace(`assigned taxonomy`, ";;", ";NA;"),          # It has to be done one at a time due to some
               `assigned taxonomy` = str_replace(`assigned taxonomy`, ";;", ";NA;"),          # wonkiness in how string replacements are done
               `assigned taxonomy` = str_replace(`assigned taxonomy`, ";;", ";NA;"),
               `assigned taxonomy` = str_replace(`assigned taxonomy`, ";;", ";NA;"),
               `assigned taxonomy` = str_replace(`assigned taxonomy`, "^Bacteria;", "")))
all_files <- lapply(all_files, function(x)
  separate(x,"assigned taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE))


compare_mock_to_actual <- function(v1,v2) {
  # This function compares the contents of two vectors v1 and v2
  # If one or both of the vectors has NA, that comparison is set to true
  # Since here, we just want to compare the number of "Corrects" to the number of actual wrong calls
  # Ambigious calls (i.e. NAs) don't count as wrong calls.
  v1[v1 == "NA"] <- NA
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- "ambiguous"
  return(same)
}


comparisons <- lapply(all_files, function(each_db) 
  cbind(# seq_name = mock$seq_name, # The output should have the seq name
        actual_taxonomy = mock$`Actual Taxonomy`,  # It should have the actual taxonomy in the mock df
        assigned_taxonomy = each_db$`assigned taxonomy`,      # It should have the full assigned taxonomy
        sapply(c("phylum", "class", "order", "family", "genus", "species"),  # And it should have the output from the six comparison calls
               function(which_rank) compare_mock_to_actual(v1 = each_db[,which_rank], v2 = mock[,which_rank]))))
comparisons <- lapply(comparisons, function(each_db)
  as.data.frame(cbind(each_db, num_levels_correct = as.character(apply(each_db, 1, function(which_row) 
    sum(as.logical(which_row[c("phylum", "class", "order", "family", "genus", "species")]), na.rm = T))))))

comparisons_df <- data.frame(lapply(comparisons, function(x) as.numeric(as.character(x$num_levels_correct))))

comparisons_df <- rbind(comparisons_df, rep(0,5))
colors <- colorRampPalette(c("white", "black"))
colors(6)


group_names <- file_names %>% str_replace(., "_taxonomy.txt","") %>% str_replace(., "16S_","") %>% str_replace(.,"_","\n")
# pdf("figures/classifier-comparisons/compare_classifiers_16S.pdf", height = 12, width = 5)
superheat(comparisons_df, membership.cols = group_names, 
          heat.pal = colors(7),
          grid.vline.col = "white", grid.vline.size = 2, 
          bottom.label.text.size = 4, bottom.label.size = .15, bottom.label.col = "white", 
          pretty.order.rows = F, pretty.order.cols = F,
          left.label = "none", left.label.size = 0, legend.width = 1, legend.breaks = 0:6)
# dev.off()


#################################
# 18S-V4
#################################

# Script to make comparison figure to compare CRUX and a few other 16s databases
# Data stored in Anacapa/data-for-comparisons
# First we make the 16s comparison
library(tidyverse)
library(reshape2)
library(superheat)
mock <- read_delim("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v4/18S_V4_actual_taxonomy.txt", delim ="\t")
# Split the taxonomy by rank
mock <- cbind(mock,colsplit(mock$`Actual Taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))

file_names <- list.files("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v4/")
# Get rid of the mock...
file_names <- file_names[!(str_detect(file_names, "actual"))]
all_files <- lapply(file_names, function(file) read_delim(file.path("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v4", file), delim = "\t"))
names(all_files) <- file_names

all_files <- lapply(all_files, function(x) 
  x %>% mutate(`assigned Taxonomy` = str_replace(`assigned Taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "[a-z]__", ""),     # Greengenes Taxonomy has p__/c__ for each level
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "D_[1-9]__", ""),   # Silva has D_1__ etc.
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, ";$", ";NA"),       # The next few lines add NAs between ;;
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # It has to be done one at a time due to some
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # wonkiness in how string replacements are done
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, "^Bacteria;", "")))
all_files <- lapply(all_files, function(x)
  separate(x,"assigned Taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE))


compare_mock_to_actual <- function(v1,v2) {
  # This function compares the contents of two vectors v1 and v2
  # If one or both of the vectors has NA, that comparison is set to true
  # Since here, we just want to compare the number of "Corrects" to the number of actual wrong calls
  # Ambigious calls (i.e. NAs) don't count as wrong calls.
  v1[v1 == "NA"] <- NA
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- "ambiguous"
  return(same)
}


comparisons <- lapply(all_files, function(each_db) 
  cbind(# seq_name = mock$seq_name, # The output should have the seq name
    actual_taxonomy = mock$`Actual Taxonomy`,  # It should have the actual taxonomy in the mock df
    assigned_taxonomy = each_db$`assigned taxonomy`,      # It should have the full assigned taxonomy
    sapply(c("phylum", "class", "order", "family", "genus", "species"),  # And it should have the output from the six comparison calls
           function(which_rank) compare_mock_to_actual(v1 = each_db[,which_rank], v2 = mock[,which_rank]))))
comparisons <- lapply(comparisons, function(each_db)
  as.data.frame(cbind(each_db, num_levels_correct = as.character(apply(each_db, 1, function(which_row) 
    sum(as.logical(which_row[c("phylum", "class", "order", "family", "genus", "species")]), na.rm = T))))))

comparisons_df <- data.frame(lapply(comparisons, function(x) as.numeric(as.character(x$num_levels_correct))))

comparisons_df <- rbind(comparisons_df, rep(0,5))
colors <- colorRampPalette(c("white", "black"))
colors(6)


group_names <- file_names %>% str_replace(., "_taxonomy.txt","") %>% str_replace(., "18S_V4_","") %>% str_replace(.,"_","\n")
# pdf("figures/classifier-comparisons/compare_classifiers_16S.pdf", height = 12, width = 5)
superheat(comparisons_df, membership.cols = group_names, 
          heat.pal = colors(7),
          grid.vline.col = "white", grid.vline.size = 2, 
          bottom.label.text.size = 4, bottom.label.size = .15, bottom.label.col = "white", 
          pretty.order.rows = F, pretty.order.cols = F,
          left.label = "none", left.label.size = 0, legend.width = 1, legend.breaks = 0:6)
# dev.off()


#################################
# 18S-V4
#################################

mock <- read_delim("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v4/18S_V4_actual_taxonomy.txt", delim ="\t")
# Split the taxonomy by rank
mock <- cbind(mock,colsplit(mock$`Actual Taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))

file_names <- list.files("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v4/")
# Get rid of the mock...
file_names <- file_names[!(str_detect(file_names, "actual"))]
all_files <- lapply(file_names, function(file) read_delim(file.path("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v4/", file), delim = "\t"))
names(all_files) <- file_names

all_files <- lapply(all_files, function(x) 
  x %>% mutate(`assigned Taxonomy` = str_replace(`assigned Taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "[a-z]__", ""),     # Greengenes Taxonomy has p__/c__ for each level
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "D_[1-9]__", ""),   # Silva has D_1__ etc.
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, ";$", ";NA"),       # The next few lines add NAs between ;;
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # It has to be done one at a time due to some
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # wonkiness in how string replacements are done
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, "^Bacteria;", "")))
all_files <- lapply(all_files, function(x)
  separate(x,"assigned Taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE))


comparisons <- lapply(all_files, function(each_db) 
  cbind(# seq_name = mock$seq_name, # The output should have the seq name
    actual_taxonomy = mock$`Actual Taxonomy`,  # It should have the actual taxonomy in the mock df
    assigned_taxonomy = each_db$`assigned Taxonomy`,      # It should have the full assigned taxonomy
    sapply(c("phylum", "class", "order", "family", "genus", "species"),  # And it should have the output from the six comparison calls
           function(which_rank) compare_mock_to_actual(v1 = each_db[,which_rank], v2 = mock[,which_rank]))))
comparisons <- lapply(comparisons, function(each_db)
  as.data.frame(cbind(each_db, num_levels_correct = as.character(apply(each_db, 1, function(which_row) 
    sum(as.logical(which_row[c("phylum", "class", "order", "family", "genus", "species")]), na.rm = T))))))

comparisons_df <- data.frame(lapply(comparisons, function(x) as.numeric(as.character(x$num_levels_correct))))

comparisons_df <- rbind(comparisons_df, rep(0,5))
colors <- colorRampPalette(c("white", "black"))
colors(6)

group_names <- file_names %>% str_replace(., "_taxonomy.txt","") %>% str_replace(., "18S_V4_","") %>% str_replace(.,"_","\n")
pdf("figures/classifier-comparisons/compare_classifiers_18SV4.pdf", height = 12, width = 5)
superheat(comparisons_df, membership.cols = group_names, 
          heat.pal = colors(6), grid.hline = F,
          grid.vline.col = "white", grid.vline.size = 2, 
          bottom.label.text.size = 4, bottom.label.size = .15, bottom.label.col = "white", 
          pretty.order.rows = F, pretty.order.cols = F,
          left.label = "none", left.label.size = 0, legend.width = 1, legend.breaks = 0:6)
dev.off()

#################################
# 18S-V8-9
#################################

mock <- read_delim("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v8-9/18S_V8-9_actual_taxonomy.txt", delim ="\t")
# Split the taxonomy by rank
mock <- cbind(mock,colsplit(mock$`Actual Taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))

file_names <- list.files("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v8-9/")
# Get rid of the mock...
file_names <- file_names[!(str_detect(file_names, "actual"))]
all_files <- lapply(file_names, function(file) read_delim(file.path("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/18S_v8-9/", file), delim = "\t"))
names(all_files) <- file_names

all_files <- lapply(all_files, function(x) 
  x %>% mutate(`assigned Taxonomy` = str_replace(`assigned Taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "[a-z]__", ""),     # Greengenes Taxonomy has p__/c__ for each level
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "D_[1-9]__", ""),   # Silva has D_1__ etc.
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, ";$", ";NA"),       # The next few lines add NAs between ;;
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # It has to be done one at a time due to some
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # wonkiness in how string replacements are done
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, "^Bacteria;", "")))
all_files <- lapply(all_files, function(x)
  separate(x,"assigned Taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE))



comparisons <- lapply(all_files, function(each_db) 
  cbind(# seq_name = mock$seq_name, # The output should have the seq name
    actual_taxonomy = mock$`Actual Taxonomy`,  # It should have the actual taxonomy in the mock df
    assigned_taxonomy = each_db$`assigned Taxonomy`,      # It should have the full assigned taxonomy
    sapply(c("phylum", "class", "order", "family", "genus", "species"),  # And it should have the output from the six comparison calls
           function(which_rank) compare_mock_to_actual(v1 = each_db[,which_rank], v2 = mock[,which_rank]))))
comparisons <- lapply(comparisons, function(each_db)
  as.data.frame(cbind(each_db, num_levels_correct = as.character(apply(each_db, 1, function(which_row) 
    sum(as.logical(which_row[c("phylum", "class", "order", "family", "genus", "species")]), na.rm = T))))))

comparisons_df <- data.frame(lapply(comparisons, function(x) as.numeric(as.character(x$num_levels_correct))))

comparisons_df <- rbind(comparisons_df, rep(0,5))
colors <- colorRampPalette(c("white", "black"))
colors(6)

group_names <- file_names %>% str_replace(., "_taxonomy.txt","") %>% str_replace(., "18S_V8-9_","") %>% str_replace(.,"_","\n")
pdf("figures/classifier-comparisons/compare_classifiers_18SV8-9.pdf", height = 12, width = 5)
superheat(comparisons_df, membership.cols = group_names, 
          heat.pal = colors(6), grid.hline = F,
          grid.vline.col = "white", grid.vline.size = 2, 
          bottom.label.text.size = 4, bottom.label.size = .15, bottom.label.col = "white", 
          pretty.order.rows = F, pretty.order.cols = F,
          left.label = "none", left.label.size = 0, legend.width = 1, legend.breaks = 0:6)
dev.off()

#################################
# CO1
#################################

mock <- read_delim("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/CO1/CO1_actual_taxonomy.txt", delim ="\t")
# Split the taxonomy by rank
mock <- cbind(mock,colsplit(mock$`Actual Taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))

file_names <- list.files("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/CO1/")
# Get rid of the mock...
file_names <- file_names[!(str_detect(file_names, "actual"))]
all_files <- lapply(file_names, function(file) read_delim(file.path("data-for-figs/SangerMock_compare-blast_v_bowtie2_blca-qiime2/CO1/", file), delim = "\t"))
names(all_files) <- file_names

all_files <- lapply(all_files, function(x) 
  x %>% mutate(`assigned Taxonomy` = str_replace(`assigned Taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "[a-z]__", ""),     # Greengenes Taxonomy has p__/c__ for each level
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "D_[1-9]__", ""),   # Silva has D_1__ etc.
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `assigned Taxonomy` = str_replace_all(`assigned Taxonomy`, ";$", ";NA"),       # The next few lines add NAs between ;;
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # It has to be done one at a time due to some
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),          # wonkiness in how string replacements are done
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, ";;", ";NA;"),
               `assigned Taxonomy` = str_replace(`assigned Taxonomy`, "^Bacteria;", "")))
all_files <- lapply(all_files, function(x)
  separate(x,"assigned Taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE))


comparisons <- lapply(all_files, function(each_db) 
  cbind(# seq_name = mock$seq_name, # The output should have the seq name
    actual_taxonomy = mock$`Actual Taxonomy`,  # It should have the actual taxonomy in the mock df
    assigned_taxonomy = each_db$`assigned Taxonomy`,      # It should have the full assigned taxonomy
    sapply(c("phylum", "class", "order", "family", "genus", "species"),  # And it should have the output from the six comparison calls
           function(which_rank) compare_mock_to_actual(v1 = each_db[,which_rank], v2 = mock[,which_rank]))))
comparisons <- lapply(comparisons, function(each_db)
  as.data.frame(cbind(each_db, num_levels_correct = as.character(apply(each_db, 1, function(which_row) 
    sum(as.logical(which_row[c("phylum", "class", "order", "family", "genus", "species")]), na.rm = T))))))

comparisons_df <- data.frame(lapply(comparisons, function(x) as.numeric(as.character(x$num_levels_correct))))

comparisons_df <- rbind(comparisons_df, rep(0,5))
colors <- colorRampPalette(c("white", "black"))
colors(6)

group_names <- file_names %>% str_replace(., "_taxonomy.txt","") %>% str_replace(., "CO1_","") %>% str_replace(.,"_","\n")
pdf("figures/classifier-comparisons/compare_classifiers_CO1.pdf", height = 12, width = 5)
superheat(comparisons_df, membership.cols = group_names, 
          heat.pal = colors(6), grid.hline = F,
          grid.vline.col = "white", grid.vline.size = 2, 
          bottom.label.text.size = 4, bottom.label.size = .15, bottom.label.col = "white", 
          pretty.order.rows = F, pretty.order.cols = F,
          left.label = "none", left.label.size = 0, legend.width = 1, legend.breaks = 0:6)
dev.off()