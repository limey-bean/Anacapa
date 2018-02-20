library(tidyverse)
library(gplots)
library(reshape2)

subdir_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/16S", pattern = "b")
file_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/16S/blast/")
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/compare-blast_v_bowtie2_blca/16S/", subdir, each_file), delim = "\t")))
names(all_files) <- subdir_names

mock <- read_delim("data-for-figs/compare-blast_v_bowtie2_blca/16S/mock_3_actual_taxonomy_summary.txt", delim ="\t")
# Fix some "Not Available" calls to just read NAs
mock$`actual taxonomy` <- str_replace(mock$`actual taxonomy`, "Not Available", "NA")
mock <- cbind(mock,colsplit(mock$`actual taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))


all_files_for_heatmap <- lapply(all_files, function(each_db) lapply(each_db, function(x)
  x %>% mutate(`sum taxonomy` = str_replace(`sum taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "[a-z]__", ""),     # Greengenes taxonomy has p__/c__ for each level
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "D_[1-9]__", ""),   # Silve has D_1__ etc.
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `sum taxonomy` = str_replace_all(`sum taxonomy`, ";$", ";NA"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"))))     
# Split up the single `sum taxonomy` column into a column for each level
all_files_for_heatmap <- lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x)
  separate(x,"sum taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE)))


compare_mock_to_actual <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- TRUE
  return(same)
}

comparisons <-lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x) 
  cbind(seq_name = mock$seq_name, actual_taxonomy = mock$`actual taxonomy`, assigned_taxonomy = x$`sum taxonomy`,
        sapply(c("phylum", "class", "order", "family", "genus", "species"), 
               function(y) compare_mock_to_actual(v1 = x[,y], v2 = mock[,y])))))
comparisons <- lapply(comparisons, function(each_db) lapply(each_db, function(x) as.data.frame(cbind(x, num_levels_idd = 
                                                                                                       as.character(apply(x, 1, function(y) 
                                                                                                         sum(as.logical(y[c("phylum", "class", "order", "family", "genus", "species")]))))))))


comparisons_for_graphs <- as.data.frame(sapply(1:length(comparisons), function(x) 
  as.numeric(as.character(comparisons[[x]]$num_levels_idd))))

comparisons_df <- data.frame(lapply(comparisons, function(each_db) lapply(each_db, function(x) as.numeric(as.character(x$num_levels_idd)))))

colnames(comparisons_df) <- paste(rep(subdir_names, each = 5), file_names, sep = "_")

bowtie_calls <- comparisons_df %>% select(matches("bowtie2"))
blast_calls <- comparisons_df %>% select(matches("blast"))
for_hm_16s <- (bowtie_calls-blast_calls)

colnames(for_hm_16s) <- paste0("16S","_",c(60,70,80,90,95))

colors <- colorRampPalette(colors = c("darkred", "white", "darkblue"))
hmcol = colors(9)
# hmcol[6] <- "white"
pdf("figures/bowtie2_v_blca_16S.pdf")
heatmap.2(as.matrix(for_hm_16s), Rowv = F, 
          Colv = F,  breaks =  seq(from=-4.5, to = 4.5, by = 1), trace ="none", col = hmcol,
          colsep = 1:ncol(for_hm_16s), rowsep = 1:nrow(for_hm_16s), sepcolor = "grey", sepwidth=c(0.0001,0.0001), 
          main = "16S", cexCol = 1)#,
# ColSideColors = rep("darkblue", 5))
dev.off()


# --------------------
# 18S ----------------
library(tidyverse)
library(gplots)
library(reshape2)

subdir_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/18S_V4", pattern = "b")
file_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/18S_V4/blast/")
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/compare-blast_v_bowtie2_blca/18S_V4/", subdir, each_file), delim = "\t") %>% arrange(seq_name) ))
names(all_files) <- subdir_names

mock <- read_delim("data-for-figs/compare-blast_v_bowtie2_blca/18S_V4/18S_actual_taxonomy_summary.txt", delim ="\t")
mock <- mock %>% arrange(seq_name)
# Fix some "Not Available" calls to just read NAs
mock$`actual taxonomy` <- str_replace(mock$`sum taxonomy`, "Not Available", "NA")
mock <- cbind(mock,colsplit(mock$`sum taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))


all_files_for_heatmap <- lapply(all_files, function(each_db) lapply(each_db, function(x)
  x %>% mutate(`sum taxonomy` = str_replace(`sum taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "[a-z]__", ""),     # Greengenes taxonomy has p__/c__ for each level
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "D_[1-9]__", ""),   # Silve has D_1__ etc.
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `sum taxonomy` = str_replace_all(`sum taxonomy`, ";$", ";NA"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"))))     
# Split up the single `sum taxonomy` column into a column for each level
all_files_for_heatmap <- lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x)
  separate(x,"sum taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE)))


compare_mock_to_actual <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- TRUE
  return(same)
}

comparisons <-lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x) 
  cbind(seq_name = mock$seq_name, actual_taxonomy = mock$`actual taxonomy`, assigned_taxonomy = x$`sum taxonomy`,
        sapply(c("phylum", "class", "order", "family", "genus", "species"), 
               function(y) compare_mock_to_actual(v1 = x[,y], v2 = mock[,y])))))
comparisons <- lapply(comparisons, function(each_db) lapply(each_db, function(x) as.data.frame(cbind(x, num_levels_idd = 
                                                                                                       as.character(apply(x, 1, function(y) 
                                                                                                         sum(as.logical(y[c("phylum", "class", "order", "family", "genus", "species")]))))))))


comparisons_for_graphs <- as.data.frame(sapply(1:length(comparisons), function(x) 
  as.numeric(as.character(comparisons[[x]]$num_levels_idd))))

comparisons_df <- data.frame(lapply(comparisons, function(each_db) lapply(each_db, function(x) as.numeric(as.character(x$num_levels_idd)))))

colnames(comparisons_df) <- paste(rep(subdir_names, each = 5), file_names, sep = "_")

bowtie_calls <- comparisons_df %>% select(matches("bowtie2"))
blast_calls <- comparisons_df %>% select(matches("blast"))
for_hm <- (bowtie_calls-blast_calls)


colors <- colorRampPalette(colors = c("darkred", "white", "darkblue"))
hmcol = colors(10)
hmcol[6] <- "white"
# heatmap.2(as.matrix(for_hm), Rowv = F, 
#           Colv = F, breaks = seq(from=-4, to = 4, by = .75), trace ="none", col = hmcol,
#           colsep = 1:ncol(for_hm), rowsep = 1:nrow(for_hm), sepcolor = "grey", sepwidth=c(0.0001,0.0001), 
#           main = "18S_V4")#,
# ColSideColors = rep("darkblue", 5))

# -----------------
# 18S-V8V9
library(tidyverse)
library(gplots)
library(reshape2)

subdir_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/18S_V8-9", pattern = "b")
file_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/18S_V8-9/blast/")
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/compare-blast_v_bowtie2_blca/18S_V8-9/", subdir, each_file), delim = "\t")))
names(all_files) <- subdir_names

mock <- read_delim("data-for-figs/compare-blast_v_bowtie2_blca/18S_V4/18S_actual_taxonomy_summary.txt", delim ="\t")
# Fix some "Not Available" calls to just read NAs
mock$`sum taxonomy` <- str_replace(mock$`sum taxonomy`, "Not Available", "NA")
mock <- cbind(mock,colsplit(mock$`sum taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))


all_files_for_heatmap <- lapply(all_files, function(each_db) lapply(each_db, function(x)
  x %>% mutate(`sum taxonomy` = str_replace(`sum taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "[a-z]__", ""),     # Greengenes taxonomy has p__/c__ for each level
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "D_[1-9]__", ""),   # Silve has D_1__ etc.
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `sum taxonomy` = str_replace_all(`sum taxonomy`, ";$", ";NA"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"))))     
# Split up the single `sum taxonomy` column into a column for each level
all_files_for_heatmap <- lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x)
  separate(x,"sum taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE)))


compare_mock_to_actual <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- TRUE
  return(same)
}

comparisons <-lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x) 
  cbind(seq_name = mock$seq_name, actual_taxonomy = mock$`actual taxonomy`, assigned_taxonomy = x$`sum taxonomy`,
        sapply(c("phylum", "class", "order", "family", "genus", "species"), 
               function(y) compare_mock_to_actual(v1 = x[,y], v2 = mock[,y])))))
comparisons <- lapply(comparisons, function(each_db) lapply(each_db, function(x) as.data.frame(cbind(x, num_levels_idd = 
                                                                                                       as.character(apply(x, 1, function(y) 
                                                                                                         sum(as.logical(y[c("phylum", "class", "order", "family", "genus", "species")]))))))))


comparisons_for_graphs <- as.data.frame(sapply(1:length(comparisons), function(x) 
  as.numeric(as.character(comparisons[[x]]$num_levels_idd))))

comparisons_df <- data.frame(lapply(comparisons, function(each_db) lapply(each_db, function(x) as.numeric(as.character(x$num_levels_idd)))))

colnames(comparisons_df) <- paste(rep(subdir_names, each = 5), file_names, sep = "_")

bowtie_calls <- comparisons_df %>% select(matches("bowtie2"))
blast_calls <- comparisons_df %>% select(matches("blast"))
for_hm2 <- (bowtie_calls-blast_calls)

for_hm3 <- cbind(for_hm, for_hm2)
colors <- colorRampPalette(colors = c("darkred", "white", "darkblue"))
hmcol = colors(9)
# hmcol[6] <- "white"
colnames(for_hm3) <- paste0(rep(c("V4", "V8-9"), each = 5), "_",rep(c(60,70,80,90,95), 2))

pdf("figures/bowtie2_v_blca_18S.pdf")
heatmap.2(as.matrix(for_hm3), Rowv = F, 
          Colv = F, breaks =  seq(from=-4.5, to = 4.5, by = 1), trace ="none", col = hmcol,
          colsep = 1:ncol(for_hm3), rowsep = 1:nrow(for_hm), sepcolor = "grey", sepwidth=c(0.0001,0.0001), 
          main = "18S_V4 and 18S_V8-9",
          ColSideColors = c(rep("darkblue", 5), rep("darkgreen",5)))
dev.off()

# ---------------
# For CO1, once it's available
library(tidyverse)
library(gplots)
library(reshape2)

subdir_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/CO1", pattern = "b")
file_names <- list.files("data-for-figs/compare-blast_v_bowtie2_blca/CO1/blast/")
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/compare-blast_v_bowtie2_blca/CO1/", subdir, each_file), delim = "\t") %>% arrange(CO1_seq_number)))
names(all_files) <- subdir_names

mock <- read_delim("data-for-figs/compare-blast_v_bowtie2_blca/CO1/mock_CO1_expected.txt", delim ="\t")
mock <- mock %>% arrange(accession_numbers)
# Fix some "Not Available" calls to just read NAs
mock$`taxonomy` <- str_replace(mock$`taxonomy`, "Not Available", "NA")
mock <- cbind(mock,colsplit(mock$`taxonomy`, ";", names = c("kingdom","phylum", "class", "order", "family", "genus", "species")))
mock <- mock %>% select(-kingdom)

all_files_for_heatmap <- lapply(all_files, function(each_db) lapply(each_db, function(x)
  x %>% mutate(`sum taxonomy` = str_replace(`sum taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "[a-z]__", ""),     # Greengenes taxonomy has p__/c__ for each level
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "D_[1-9]__", ""),   # Silve has D_1__ etc.
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `sum taxonomy` = str_replace_all(`sum taxonomy`, ";$", ";NA"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"))))     
# Split up the single `sum taxonomy` column into a column for each level
all_files_for_heatmap <- lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x)
  separate(x,"sum taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE)))


compare_mock_to_actual <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- TRUE
  return(same)
}

comparisons <-lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x) 
  cbind(seq_name = mock$seq_name, actual_taxonomy = mock$`actual taxonomy`, assigned_taxonomy = x$`sum taxonomy`,
        sapply(c("phylum", "class", "order", "family", "genus", "species"), 
               function(y) compare_mock_to_actual(v1 = x[,y], v2 = mock[,y])))))
comparisons <- lapply(comparisons, function(each_db) lapply(each_db, function(x) as.data.frame(cbind(x, num_levels_idd = 
                                                                                                       as.character(apply(x, 1, function(y) 
                                                                                                         sum(as.logical(y[c("phylum", "class", "order", "family", "genus", "species")]))))))))


comparisons_for_graphs <- as.data.frame(sapply(1:length(comparisons), function(x) 
  as.numeric(as.character(comparisons[[x]]$num_levels_idd))))

comparisons_df <- data.frame(lapply(comparisons, function(each_db) lapply(each_db, function(x) as.numeric(as.character(x$num_levels_idd)))))

colnames(comparisons_df) <- paste(rep(subdir_names, each = 5), file_names, sep = "_")

bowtie_calls <- comparisons_df %>% select(matches("bowtie2"))
blast_calls <- comparisons_df %>% select(matches("blast"))
for_hm <- (bowtie_calls-blast_calls)


colors <- colorRampPalette(colors = c("darkred", "white", "darkblue"))
hmcol = colors(10)
hmcol[6] <- "white"
# pdf("figures/bowtie2_v_blca_CO1.pdf")
heatmap.2(as.matrix(for_hm), Rowv = F,
          Colv = F, breaks = seq(from=-4, to = 4, by = .75), trace ="none", col = hmcol,
          colsep = 1:ncol(for_hm), rowsep = 1:nrow(for_hm), sepcolor = "grey", sepwidth=c(0.0001,0.0001), 
          main = "CO1")#,
# ColSideColors = rep("darkblue", 5))
# dev.off()

