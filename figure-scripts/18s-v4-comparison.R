library(dplyr)
library(tidyverse)
library(superheat)
subdir_names <- list.files("data-for-figs/Database_comparisons/compare-18S-databases/18S_V4/", pattern = "coded")
file_names <- list.files("data-for-figs/Database_comparisons/compare-18S-databases/18S_V4/Silva_V4_18S_coded/")
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/Database_comparisons/compare-18S-databases/18S_V4/", subdir, each_file), delim = "\t") %>% 
    select(-mock_sanger, correct = `TRUE`, wrong = `FALSE`)))
names(all_files) <- subdir_names


# We need to extract the `correct` columns and make a matrix out of it
comparisons_df <- data.frame(lapply(all_files, function(each_db) lapply(each_db, function(x) 
  as.numeric(as.character(x$correct)))))
colnames(comparisons_df) <- paste(rep(subdir_names, each = 5), file_names, sep = "_")

group_names <- c("CRUX\nFiltered 18S-V4","CRUX\nUnfiltered 18S-V4", "Silva")


colors <- colorRampPalette(c("white", "black"))
colors(6)
pdf("figures/heatmap-18s_V4.pdf", height = 15, width = 22)
superheat(comparisons_df, membership.cols = rep(group_names, each = 5), 
          heat.pal = colors(7),
          grid.vline.col = "white", grid.vline.size = 2, 
          bottom.label.text.size = 7, bottom.label.size = .15, bottom.label.col = "white", 
          pretty.order.rows = F, pretty.order.cols = F,
          left.label = "none", left.label.size = 0, legend.width = 1, legend.breaks = 0:6)
dev.off()

#-------------

barplot_of_empties <- do.call(cbind, lapply(all_files, function(each_db) 
  sapply(each_db, function(x)  colSums(x[,2:4]))))
barplot_of_empties <- barplot_of_empties/colSums(barplot_of_empties)
rownames(barplot_of_empties) <- c("correct", "wrong", "ambiguous call")
colnames(barplot_of_empties) <- paste0(rep(c("CRUX_filtered_", "CRUX_unfiltered_", "Silva_"),each = 5), c(60,70,80,90,95))
write.csv(barplot_of_empties, "figures/accuracy_percentages_per_barcode/accuracy_18s_V4.csv")
colnames(barplot_of_empties) <- NULL


pdf("figures/18s_v4_database_accuracy_comparisons.pdf", width = 22, height = 7)
par(xpd = T, mar = c(5.1,7.1,4.1,2.1), oma = c(0,1,0,0))
barplot(barplot_of_empties, cex.axis = 2, cex.lab = 2,
        ylab = "Proportion of taxonomic calls done correctly",
        col = c("black", "grey", "white"))
legend(1,1.25, legend = rownames(barplot_of_empties), horiz = T, x.intersp = .1, bty = "n", 
       text.width = c(0,4.5,4), fill = c("black", "grey", "white"))
dev.off()
