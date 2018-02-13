library(dplyr)
library(readr)
library(stringr)
library(reshape2)
library(tidyr)
subdir_names <- list.files("data-for-figs/compare-18S-databases/18S_V4/", pattern = "taxonomy_table")
file_names <- list.files("data-for-figs/compare-18S-databases/18S_V4/Crux_filtered_18S_V4-taxonomy_table/")
all_files <- lapply(subdir_names, function(subdir) lapply(file_names, function(each_file) 
  read_delim(file.path("data-for-figs/compare-18S-databases/18S_V4/", subdir, each_file), delim = "\t")))
names(all_files) <- subdir_names

mock <- read_delim("data-for-figs/compare-18S-databases/18S_V4/18S_V4_taxonomy.txt", delim ="\t")
mock$`actual taxonomy` <- str_replace(mock$`actual taxonomy`, "Not Available", "NA")
mock <- cbind(mock,colsplit(mock$`actual taxonomy`, ";", names = c("phylum", "class", "order", "family", "genus", "species")))
mock <- mock %>% arrange(seq_name)
# ---------



all_files <- lapply(all_files, function(each_database) lapply(each_database, function(each_file) arrange(each_file, seq_name)))

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


group_names <- c("CRUX\nFiltered 18S-V4","CRUX\nUnfiltered 18S-V4", "Silva")




#-------------
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
comparisons_e
colnames(barplot_of_empties) <- paste0(rep(names(comparisons_e), each = 5), 
                                       c("60", "70", "80", "90", "95"))

barplot_of_empties <- barplot_of_empties/colSums(barplot_of_empties)
# barplot_of_empties <- barplot_of_empties %>% tbl_df %>% select(6:10, 1:5, 11:15, 16:20) %>% as.matrix
rownames(barplot_of_empties) <- c("correct", "wrong", "ambiguous call")

# Add the Silva stuff
v418s_file_list <- list.files("data-for-figs/compare-18S-databases/18S_V4/Silva_V4_18S_coded/")
silva_performance <- lapply(v418s_file_list, function(x) 
  read.table(paste0("data-for-figs/compare-18S-databases/18S_V4/Silva_V4_18S_coded/",x), header = T, sep = "\t",stringsAsFactors = F))
names(silva_performance) <- v418s_file_list
v418s_silv_bar <- sapply(silva_performance, function(x) colSums(x[,3:5]))
v418s_silv_bar <- v418s_silv_bar/colSums(v418s_silv_bar)
rownames(v418s_silv_bar) <- c("correct", "wrong", "ambiguous call")
barplot_of_empties <- cbind(barplot_of_empties,v418s_silv_bar )
colnames(barplot_of_empties) <- NULL

# pdf("figures/18s_v4_database_accuracy_comparisons_v1.pdf", height = 7, width = 14)
barplot(barplot_of_empties, legend = rownames(barplot_of_empties),
        args.legend = list(x = "topright", bty = "n", inset=c(0, -0.1), 
                           horiz = T),
        ylab = "Proportion of taxonomic calls done correctly",
        col = c("black", "grey", "white"), main = "18S-V4")
# dev.off()


comparisons_df <- cbind(comparisons_df, sapply(silva_performance, function(x) x[,"TRUE."]))


# pdf("figures/heatmap-18s_v4_heatmap.pdf", height = 15, width = 22)
superheat(comparisons_df, membership.cols = rep(c("CRUX\n18S-V4 Filtered","CRUX\n18S-V4 Unfiltered","Silva"), each = 5), 
          heat.pal = colors(7),
          grid.vline.col = "white", grid.vline.size = 2, bottom.label.text.size = 7,
          bottom.label.size = .15,
          legend.breaks = 0:6, bottom.label.col = "white", pretty.order.rows = F, pretty.order.cols = F,
          X.text.size = 0,X.text = as.matrix(comparisons_df),# bottom.label.names = group_names,
          # yt = (colSums(comparisons_df)/nrow(comparisons_df)/6), yt.plot.type = "bar",
          # yt = barplot_of_empties, yt.plot.type = "bar",
          # yt.bar.col = "black",yt.obs.col = rep("grey", length(all_files_for_heatmap)*5), yt.point.size = 1.25, yt.num.ticks = 6,
          # yt.axis.name = "Average percentage\nof taxonomy\nassigned correctly", 
          left.label = "none", 
          left.label.size = 0, yt.axis.size = 20, yt.axis.name.size = 20,yt.lim = c(0,1))
# dev.off()


# Try to redo all of this cleanly
# to_plot_ggbar <- barplot_of_empties %>% t %>% data.frame() %>% tibble::rownames_to_column() %>% gather(.,db, percentage, -rowname) 

# With new V4 data --------------
v418s_file_list_1 <- list.files("data-for-figs/compare-18S-databases/18S_V4/CRUX_filtered_V4_18S_coded/")
crux_f_performance <- lapply(v418s_file_list_1, function(x) 
  read.table(paste0("data-for-figs/compare-18S-databases/18S_V4/CRUX_filtered_V4_18S_coded/",x), header = T, sep = "\t",stringsAsFactors = F))
names(crux_f_performance) <- v418s_file_list_1
v418s_cruxf_bar <- sapply(crux_f_performance, function(x) colSums(x[,3:5]))
v418s_cruxf_bar <- v418s_cruxf_bar/colSums(v418s_cruxf_bar)
rownames(v418s_cruxf_bar) <- c("correct", "wrong", "ambiguous call")


v418s_file_list_2 <- list.files("data-for-figs/compare-18S-databases/18S_V4/CRUX_unfiltered_V4_18S_coded/")
crux_uf_performance <- lapply(v418s_file_list_2, function(x) 
  read.table(paste0("data-for-figs/compare-18S-databases/18S_V4/CRUX_unfiltered_V4_18S_coded/",x), header = T, sep = "\t",stringsAsFactors = F))
# names(crux_uf_performance) <- v418s_file_list_2
v418s_cruxuf_bar <- sapply(crux_uf_performance, function(x) colSums(x[,3:5]))
v418s_cruxuf_bar <- v418s_cruxuf_bar/colSums(v418s_cruxuf_bar)
rownames(v418s_cruxuf_bar) <- c("correct", "wrong", "ambiguous call")

for_full_bar <- cbind(v418s_cruxf_bar,v418s_cruxuf_bar,v418s_silv_bar)

pdf("figures/18s_v4_database_accuracy_comparisons.pdf")
barplot(for_full_bar, legend = rownames(barplot_of_empties),
        args.legend = list(x = "topright", bty = "n", inset=c(0, -0.1), 
                           horiz = T),
        ylab = "Proportion of taxonomic calls done correctly",
        col = c("black", "grey", "white"), main = "18S-V4")
dev.off()


for_heatmap <- sapply(crux_f_performance, function(x) x[,"TRUE."])
for_heatmap <- cbind(for_heatmap,sapply(crux_uf_performance, function(x) x[,"TRUE."]))
for_heatmap <- cbind(for_heatmap,sapply(silva_performance, function(x) x[,"TRUE."]))
pdf("figures/heatmap-18s_v4.pdf")
superheat(for_heatmap, membership.cols = rep(c("CRUX 18S V4\nFiltered","CRUX 18S V4\nUnfiltered",
                                               "Silva 18S V4"), each = 5), 
          heat.pal = colors(7),
          grid.vline.col = "white", grid.vline.size = 2, bottom.label.text.size = 7,
          bottom.label.size = .15,
          legend.breaks = 0:6, bottom.label.col = "white", pretty.order.rows = F, pretty.order.cols = F,
          X.text.size = 0,#X.text = as.matrix(comparisons_df),# bottom.label.names = group_names,
          # yt = (colSums(comparisons_df)/nrow(comparisons_df)/6), yt.plot.type = "bar",
          # yt = barplot_of_empties, yt.plot.type = "bar",
          # yt.bar.col = "black",yt.obs.col = rep("grey", length(all_files_for_heatmap)*5), yt.point.size = 1.25, yt.num.ticks = 6,
          # yt.axis.name = "Average percentage\nof taxonomy\nassigned correctly", 
          left.label = "none", 
          left.label.size = 0, yt.axis.size = 20, yt.axis.name.size = 20,yt.lim = c(0,1))
dev.off()