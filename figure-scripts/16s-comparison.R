# Script to make comparison figure to compare CRUX and a few other 16s databases
# Data stored in Anacapa/data-for-comparisons


library(tidyverse)
library(reshape2)
library(superheat)

# Read in the mock dataset
mock <- read_delim("data-for-figs/Database_comparisons/compare-16s-databases/mock_3_actual_taxonomy_summary.txt", delim ="\t")
# Split the taxonomy by rank
mock <- cbind(mock,colsplit(mock$`actual taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))
# Get rid of extraneous columns
mock <- mock %>% select(-Andrew, -HMP, -Extreme)

# Import all the files to compare. The structure is slightly convoluted, so a brief explanation:
# All files will be saved in a single list, named "all_files_for_heatmap_for_heatmap"
# The list will have as many elements as there are databases (as of this writing, 7 databases)
# Each of the 7 elements is in turn a list, with 5 elements each.
# Each of the 5 is one percent_confidence cutoff (60,70,80,90,95).
subdir_names <- list.files("data-for-figs/Database_comparisons/compare-16s-databases/", pattern = "taxonomy-tables")
file_names <- list.files("data-for-figs/Database_comparisons/compare-16s-databases/filtered-16s-taxonomy-tables/")

# Read in all the files into a nested list with the data structure described above
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/Database_comparisons/compare-16s-databases/", subdir, each_file), delim = "\t") %>%
    select(-Andrew, -HMP, -Extreme))) # Get rid of the Andrew, HMP, and Extreme columns from these too, just to keep things clean
names(all_files) <- subdir_names

# Take care of some quirks in taxonomy paths
all_files_for_heatmap <- lapply(all_files, function(each_db) lapply(each_db, function(x)
  x %>% mutate(`sum taxonomy` = str_replace(`sum taxonomy`, "Not Available", "NA"), # Fix the Not Available issue, if it persists
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "[a-z]__", ""),     # Greengenes taxonomy has p__/c__ for each level
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "D_[1-9]__", ""),   # Silva has D_1__ etc.
               `sum taxonomy` = str_replace_all(`sum taxonomy`, "^;", "NA;"),       # Replace first empty with NA
               `sum taxonomy` = str_replace_all(`sum taxonomy`, ";$", ";NA"),       # The next few lines add NAs between ;;
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),          # It has to be done one at a time due to some
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),          # wonkiness in how string replacements are done
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"),
               `sum taxonomy` = str_replace(`sum taxonomy`, ";;", ";NA;"))))     



# Split up the single `sum taxonomy` column by taxonomic rank
all_files_for_heatmap <- lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(x)
  separate(x,"sum taxonomy", sep = ";", into = c("phylum", "class", "order", "family", "genus", "species"), remove = FALSE)))

compare_mock_to_actual <- function(v1,v2) {
  # This function compares the contents of two vectors v1 and v2
  # If one or both of the vectors has NA, that comparison is set to true
  # Since here, we just want to compare the number of "Corrects" to the number of actual wrong calls
  # Ambigious calls (i.e. NAs) don't count as wrong calls.
  # NOTE! Check with EC if this is the right way to think about this. The other way around
  # might be more appropriate.
  # Actually, this might just be the number correct, ambiguous is considered wrong
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- TRUE
  return(same)
}

# Compare each taxonomic rank in each file to the mock community
comparisons <- lapply(all_files_for_heatmap, function(each_db) lapply(each_db, function(which_cutoff) 
    cbind(seq_name = mock$seq_name, # The output should have the seq name
          actual_taxonomy = mock$`actual taxonomy`,  # It should have the actual taxonomy in the mock df
          assigned_taxonomy = which_cutoff$`sum taxonomy`,      # It should have the full assigned taxonomy
          sapply(c("phylum", "class", "order", "family", "genus", "species"),  # And it should have the output from the six comparison calls
                 function(which_rank) compare_mock_to_actual(v1 = which_cutoff[,which_rank], v2 = mock[,which_rank])))))

# Now we just want to count the number of TRUEs for each row.
comparisons <- lapply(comparisons, function(each_db) lapply(each_db, function(which_cutoff) 
  as.data.frame(cbind(which_cutoff, num_levels_correct = as.character(apply(which_cutoff, 1, function(which_row) 
    sum(as.logical(which_row[c("phylum", "class", "order", "family", "genus", "species")]))))))))

# Let's do some spot checks...
comparisons[[3]][[1]]$actual_taxonomy[9]
comparisons[[3]][[1]]$assigned_taxonomy[9]
comparisons[[3]][[1]]$num_levels_correct[9] == 6

comparisons_for_graphs <- as.data.frame(sapply(1:length(comparisons), function(x) 
  as.numeric(as.character(comparisons[[x]]$num_levels_correct))))

comparisons_df <- data.frame(lapply(comparisons, function(each_db) lapply(each_db, function(x) as.numeric(as.character(x$num_levels_correct)))))

colnames(comparisons_df) <- paste(rep(subdir_names, each = 5), file_names, sep = "_")


rownames(comparisons_df) <- mock$seq_name
colors <- colorRampPalette(c("white", "black"))
colors(6)

# 
# group_names <- c("Blast BLCA\n80ID, 0 Over", "Blast BLCA\n90ID 10 Over", "CRUX-Blast-BLCA\n80ID 100 Returns\n0 Over",
#                  "Filtered 16S", "Greengenes", "Silva", "Unfiltered 16S")
group_names <- c("CRUX\nUnfiltered 16S",
                 "CRUX\nFiltered 16S", "Greengenes", "Silva")


pdf("figures/heatmap-16s_crux_v_greengenes_v_silva.pdf", height = 15, width = 22)
superheat(comparisons_df, membership.cols = rep(group_names, each = 5), 
          heat.pal = colors(7),
          grid.vline.col = "white", grid.vline.size = 2, 
          bottom.label.text.size = 7, bottom.label.size = .15, bottom.label.col = "white", 
          pretty.order.rows = F, pretty.order.cols = F,
          left.label = "none", left.label.size = 0, legend.width = 1, legend.breaks = 0:6)
dev.off()

# Make the histogram include the unkown calls from Anacapa ----
# We have to make this all over again so that instead of calling things "NA"
# We can call them 'empty' to be explicit
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

# This function is better at dealing with the ambiguous alls
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

# Now we need to make a column for the number of Trues, the number of Empties, and the number of Falses
comparisons_e <- lapply(comparisons_e, function(each_db) 
  lapply(each_db, function(x) as.data.frame(cbind(#x, 
                                                  num_levels_correct = apply(x, 1, function(y) 
                                                      sum(y[c("phylum", "class", "order", "family", "genus", "species")] == "TRUE")),
                                                  num_levels_wrong = apply(x, 1, function(y) 
                                                      sum(y[c("phylum", "class", "order", "family", "genus", "species")] == "FALSE")),
                                                  num_levels_ambig = apply(x, 1, function(y) 
                                                      sum(y[c("phylum", "class", "order", "family", "genus", "species")] == "empty"))
                                                  ))))
# This do.call just converts the list above into a usable matrix to make a barplot with
barplot_of_empties <- do.call(cbind, lapply(comparisons_e, function(each_db) sapply(each_db, function(x) 
  c(sum(x[,"num_levels_correct"]),
    sum(x[,"num_levels_wrong"]),
    sum(x[,"num_levels_ambig"])))))
colnames(barplot_of_empties) <- paste0(rep(names(comparisons_e), each = 5), 
                                       c("60", "70", "80", "90", "95"))


barplot_of_empties <- barplot_of_empties/colSums(barplot_of_empties)
barplot_of_empties <- barplot_of_empties %>% tbl_df %>% select(6:10, 1:5, 11:15, 16:20) %>% as.matrix
rownames(barplot_of_empties) <- c("Correct", "Wrong", "Ambiguous")
write.csv(barplot_of_empties, "figures/accuracy_percentages_per_barcode/accuracy_16s.csv")

colnames(barplot_of_empties) <- NULL
pdf("figures/16s_database_accuracy_comparisons_v1.pdf", height = 7, width = 22)
par(xpd = T, mar = c(5.1,7.1,4.1,2.1), oma = c(0,1,0,0))
barplot(barplot_of_empties, 
        ylab = "Proportion of taxonomic calls\n(Assigned taxonomy compared to mock)",
        col = c("black", "grey", "white"), cex.axis = 2, cex.lab = 2)
legend(1,1.25, legend = rownames(barplot_of_empties), horiz = T, x.intersp = .1, bty = "n", 
       text.width = c(0,4.5,4), fill =  c("black", "grey", "white"))
dev.off()


