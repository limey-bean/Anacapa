library(readr)
library(dplyr); library(reshape2); library(stringr); library(gplots)
mock <- read_delim("~/Desktop/mock_3_actual_taxonomy_summary.txt", delim ="\t")

file_names <- list.files("~/Desktop/", pattern = "16S*")

actuals <- lapply(file_names, function(x) read_delim(file = file.path("~/Desktop", x), delim = "\t"))

names(actuals) <- file_names
mock$`actual taxonomy` <- str_replace(mock$`actual taxonomy`, "Not Available", "NA")
actuals <- lapply(actuals, function(x) x %>% mutate(`sum taxonomy` = str_replace(`sum taxonomy`, "Not Available", "NA")))

mock <- cbind(mock,colsplit(mock$`actual taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))
actuals <- lapply(actuals, function(x) cbind(x, colsplit(x$`sum taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species"))))

compare_mock_to_actual <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}


comparisons <-lapply(actuals, function(x) cbind(seq_name = mock$seq_name, actual_taxonomy = mock$`actual taxonomy`, assigned_taxonomy = x$`sum taxonomy`,
                                                sapply(c("phylum", "class", "order", "family", "genus", "species"), 
                                                       function(y) compare_mock_to_actual(v1 = x[,y], v2 = mock[,y]))))

comparisons <- lapply(comparisons, function(x) as.data.frame(cbind(x, num_levels_idd = 
                                                                     as.numeric(apply(x, 1, function(y) 
                                                                       sum(as.logical(y[c("phylum", "class", "order", "family", "genus", "species")])))))))


comparisons_for_graphs <- as.data.frame(sapply(1:length(comparisons), function(x) 
  as.numeric(as.character(comparisons[[x]]$num_levels_idd))))
colnames(comparisons_for_graphs) <- names(comparisons)
rownames(comparisons_for_graphs) <- mock$seq_name
colors <- colorRampPalette(c("white", "royalblue4"))
colors(6)
heatmap.2(as.matrix(comparisons_for_graphs), Rowv = NA, Colv = NA, col = colors(6), 
          trace = "none", density.info = "none", breaks = 0:6, srtCol = 45, cexCol = .5, colsep = 3,
          ColSideColors = c(rep("red",3),rep("darkgreen",2)))
