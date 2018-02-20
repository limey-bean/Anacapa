# Script to make comparison figure to compare CRUX and a few other 16s databases
# Data stored in Anacapa/data-for-comparisons


library(tidyverse)
library(gplots)
library(reshape2)
# NOTE!
# Using a modified version of superheat::superheat...
# Need to figure out how to make this smoother
library(superheat)

# Read in the mock dataset
mock <- read_delim("data-for-figs/compare-16s-databases/mock_3_actual_taxonomy_summary.txt", delim ="\t")
# Fix some "Not Available" calls to just read NAs
mock$`actual taxonomy` <- str_replace(mock$`actual taxonomy`, "Not Available", "NA")
mock <- cbind(mock,colsplit(mock$`actual taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))


# Import all the files to compare. The structure is slightly convoluted, so a brief explanation:
# All files will be saved in a single list, named "all_files_for_heatmap_for_heatmap"
# The list will have as many elements as there are databases (as of this writing, 7 databases)
# Each of the 7 elements is in turn a list, with 5 elements each.
# Each of the 5 is one percent_confidence cutoff (60,70,80,90,95).
subdir_names <- list.files("data-for-figs/compare-16s-databases/", pattern = "taxonomy-tables")
file_names <- list.files("data-for-figs/compare-16s-databases/filtered-16s-taxonomy-tables/")
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/compare-16s-databases/", subdir, each_file), delim = "\t")))
names(all_files) <- subdir_names

# Take care of some quirks in taxonomy paths
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


rownames(comparisons_df) <- mock$seq_name
colors <- colorRampPalette(c("white", "black"))
colors(6)

heatmap.2(as.matrix(comparisons_df), Rowv = NA, Colv = NA, col = colors(6),
          trace = "none", density.info = "none", srtCol = 45, cexCol = .5, colsep = seq(from = 5, to = 30, by = 5),
          lhei = c(1,7,1), lwid = c(0.1,10), lmat = rbind(c(0,0),c(2,1), c(3,4)), row)

# 
# group_names <- c("Blast BLCA\n80ID, 0 Over", "Blast BLCA\n90ID 10 Over", "CRUX-Blast-BLCA\n80ID 100 Returns\n0 Over",
#                  "Filtered 16S", "Greengenes", "Silva", "Unfiltered 16S")
group_names <- c("CRUX\nUnfiltered 16S",
                 "CRUX\nFiltered 16S", "Greengenes", "Silva")


pdf("figures/heatmap-16s_crux_v_greengenes_v_silva.pdf", height = 15, width = 22)
superheat.2(comparisons_df, membership.cols = rep(c(2,1,3,4), each = 5), 
            heat.pal = colors(7),
            grid.vline.col = "white", grid.vline.size = 2, bottom.label.text.size = 7,
            bottom.label.size = .15,
            legend.breaks = 0:6, bottom.label.col = "white", pretty.order.rows = F, pretty.order.cols = F,
            X.text.size = .15,X.text = as.matrix(comparisons_df), bottom.label.names = group_names,
            # yt = (colSums(comparisons_df)/nrow(comparisons_df)/6), yt.plot.type = "bar",
            # yt = barplot_of_empties, yt.plot.type = "bar",
            # yt.bar.col = "black",yt.obs.col = rep("grey", length(all_files_for_heatmap)*5), yt.point.size = 1.25, yt.num.ticks = 6,
            # yt.axis.name = "Average percentage\nof taxonomy\nassigned correctly", left.label = "none", 
            left.label.size = 0, yt.axis.size = 20, yt.axis.name.size = 20,yt.lim = c(0,1))
dev.off()

# Make the histogram include the unkown calls from Anacapa ----
all_files_for_barplot <- lapply(all_files, function(each_db) lapply(each_db, function(x)
  x %>% mutate(`sum taxonomy` = str_replace(`sum taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "[a-z]__", ""),     # Greengenes taxonomy has p__/c__ for each level
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "D_[1-9]__", ""),   # Silve has D_1__ etc.
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "^;", "empty;"),       # Replace first empty with NA
               `sum taxonomy` = str_replace_all(`sum taxonomy`, ";$", ";empty"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";empty;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";empty;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";empty;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";empty;"))))     
all_files_for_barplot <- lapply(all_files_for_barplot, function(each_db) lapply(each_db, function(x)
  separate(x,"sum taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE)))

compare_mock_to_actual_with_empties <- function(mock_taxonomy,assigned_taxonomy) {
  same <- (mock_taxonomy == assigned_taxonomy) | (is.na(mock_taxonomy) & is.na(assigned_taxonomy))
  same[is.na(same)] <- TRUE
  same[which(assigned_taxonomy=="empty")] <- "empty"
  return(same)
}

comparisons_e <-lapply(all_files_for_barplot, function(each_db) lapply(each_db, function(x) 
  cbind(seq_name = mock$seq_name, actual_taxonomy = mock$`actual taxonomy`, assigned_taxonomy = x$`sum taxonomy`,
        sapply(c("phylum", "class", "order", "family", "genus", "species"), 
               function(y) compare_mock_to_actual_with_empties(assigned_taxonomy = x[,y], mock_taxonomy = mock[,y])))))
comparisons_e <- lapply(comparisons_e, function(each_db) 
  lapply(each_db, function(x) as.data.frame(cbind(#x, 
                                                  num_levels_correct = 
                                                    apply(x, 1, function(y) 
                                                      sum(y[c("phylum", "class", "order", "family", "genus", "species")] == "TRUE")),
                                                  num_levels_wrong = 
                                                    apply(x, 1, function(y) 
                                                      sum(y[c("phylum", "class", "order", "family", "genus", "species")] == "FALSE")),
                                                  num_levels_ambig = 
                                                    apply(x, 1, function(y) 
                                                      sum(y[c("phylum", "class", "order", "family", "genus", "species")] == "empty"))
                                                  ))))
barplot_of_empties <- do.call(cbind, lapply(comparisons_e, function(each_db) sapply(each_db, function(x) 
  c(sum(x[,"num_levels_correct"]),
    sum(x[,"num_levels_wrong"]),
    sum(x[,"num_levels_ambig"])))))
colnames(barplot_of_empties) <- paste0(rep(names(comparisons_e), each = 5), 
                                       c("60", "70", "80", "90", "95"))

barplot_of_empties <- barplot_of_empties/colSums(barplot_of_empties)
barplot_of_empties <- barplot_of_empties %>% tbl_df %>% select(6:10, 1:5, 11:15, 16:20) %>% as.matrix
rownames(barplot_of_empties) <- c("correct", "wrong", "ambiguous call")
colnames(barplot_of_empties) <- NULL
# pdf("figures/16s_database_accuracy_comparisons_v1.pdf", height = 7, width = 14)
barplot(barplot_of_empties, legend = rownames(barplot_of_empties),
        args.legend = list(x = "topright", bty = "n", inset=c(0, -0.1), 
                           horiz = T),
        ylab = "Proportion of taxonomic calls done correctly",
        col = c("black", "grey", "white"))
# dev.off()
