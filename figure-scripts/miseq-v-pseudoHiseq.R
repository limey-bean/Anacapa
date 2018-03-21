# Script to make MiSeq vs Pseudo-HiSeq comparisons figure
# GSK 20 March 2018

allCO1 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/allCO1.txt", header = T)
all16S <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/all_16S.txt", header = T)
mc11_v4 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC11_v4_all.txt", header = T)
mc11_v8 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC11_v89_all.txt", header = T)
mc21_v4 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC21_v4_all.txt", header = T)
mc21_v8 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC21_v89_all.txt", header = T)
mc51_v4 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC51_v4_all.txt", header = T)
mc51_v8 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC51_v89_all.txt", header = T)
mc61_v4 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC61_v4_all.txt", header = T)
mc61_v8 <- read.table("data-for-figs/Mock_Mi_vs_Hi/data/MC61_v89_all.txt", header = T)

co1_hi <- allCO1$Hi
co1_mi <- allCO1$Mi
S16_hi <- all16S$hi+1
S16_mi <- all16S$mi+1
mc11_v4_hi <- mc11_v4$hi
mc11_v4_mi <- mc11_v4$mi
mc11_v8_hi <- mc11_v8$hi
mc11_v8_mi <- mc11_v8$mi

mc21_v4_hi <- mc21_v4$hi
mc21_v4_mi <- mc21_v4$mi
mc21_v8_hi <- mc21_v8$hi
mc21_v8_mi <- mc21_v8$mi

mc51_v4_hi <- mc51_v4$hi
mc51_v4_mi <- mc51_v4$mi
mc51_v8_hi <- mc51_v8$hi
mc51_v8_mi <- mc51_v8$mi

mc61_v4_hi <- mc61_v4$hi
mc61_v4_mi <- mc61_v4$mi
mc61_v8_hi <- mc61_v8$hi
mc61_v8_mi <- mc61_v8$mi

# Write a function that creates all the figures
compare_hi_mi <- function(hi, mi, log_ = F, title) {
  fit <- lm(hi~mi)

  greybg = ggplot2::alpha("black", .2)
  if(log_ == T){
    plot(mi,hi,col="black",
         cex=3,pch = 21, bg = greybg, log = "xy", main = title, xlab = "", ylab = "")
    mtext("Note logged axes", adj = 1, cex = 0.8)
  } else {
    plot(mi,hi,col="black",
         cex=3,pch = 21, bg = greybg, main = title, xlab = "", ylab = "")
  }
  
  abline(0,1, col="grey25", lty=2)
  abline(fit, col ="black", untf = T)
  cf <- round(coef(fit), 2) 
  eq <- paste0("y = ", 
               cf[2], " x",
               ifelse(sign(cf[1])==1, " + ", " - "), abs(cf[1]),"\n",
               "R^2 = ", round(summary(fit)$adj.r.squared, 2), 
               ", p = ", round(summary(fit)$coefficients[2,4], 6),
               ", N = ", nrow(fit$model))
  
  legend("topleft", legend = eq, bty = "n", inset = c(0, -.45), xpd = T)
}

pdf("~/Desktop/test.pdf", height = 10, width = 7, dpi = 300)
par(mfrow = c(5,2), oma = c(1,1,3,0), xpd = F)
compare_hi_mi(co1_hi, co1_mi, title = "CO1")
compare_hi_mi(S16_hi, S16_mi, log_ = T, title = "16S")
compare_hi_mi(mc11_v4_hi, mc11_v4_mi, title = "MC 11 V4")
compare_hi_mi(mc11_v8_hi, mc11_v8_mi, title = "MC 11 V8-9")

compare_hi_mi(mc21_v4_hi, mc21_v4_mi, title = "MC 21 V4")
compare_hi_mi(mc21_v8_hi, mc21_v8_mi, title = "MC 21 V8-9")

compare_hi_mi(mc51_v4_hi, mc51_v4_mi, title = "MC 51 V4")
compare_hi_mi(mc51_v8_hi, mc51_v8_mi, title = "MC 51 V8-9")

compare_hi_mi(mc61_v4_hi, mc61_v4_mi, title = "MC 61 V4")
compare_hi_mi(mc61_v8_hi, mc61_v8_mi, title = "MC 61 V8-9")
mtext(expression(bold("Taxonomy compared at the Species level")), outer = T, cex = 1.25)
mtext(side = 1, "Relative Abundance in MiSeq length (%)", outer = TRUE, cex = 1.25, line = -1)
mtext(side = 2, "Relative Abundance in HiSeq length (%)", outer = TRUE, cex = 1.25, line = -2)
dev.off()
##################################
library(tidyverse)

summarize_to_genus <- function(df) {
  df <- df %>% separate(taxon, into = c("phylum", "class", "order", "family", "genus", "species"), sep = ";")
  df <- df %>% group_by(phylum, class, order, family, genus) %>% summarize_if(is.numeric, sum)
  return(df)
}

CO1 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/allCO1_sp_max.txt", delim = "\t")
S16 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/all_16S_taxon.txt", delim = "\t")
mc11_v4 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC11_v4_all_sp_max.txt", delim = "\t")
mc11_v8 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC11_v89_all_sp_max.txt", delim = "\t")
mc21_v4 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC21_v4_all_sp_max.txt", delim = "\t")
mc21_v8 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC21_v89_all_sp_max.txt", delim = "\t")
mc51_v4 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC51_v4_all_sp_max.txt", delim = "\t")
mc51_v8 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC51_v89_all_sp_max.txt", delim = "\t")
mc61_v4 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC61_v4_all_sp_max.txt", delim = "\t")
mc61_v8 <- read_delim("data-for-figs/Mock_Mi_vs_Hi/data/MC61_v89_all_sp_max.txt", delim = "\t")

CO1_g <- summarize_to_genus(CO1)
S16_g <- summarize_to_genus(S16)

mc11_v4_g <- summarize_to_genus(mc11_v4)
mc11_v8_g <- summarize_to_genus(mc11_v8)

mc21_v4_g <- summarize_to_genus(mc21_v4)
mc21_v8_g <- summarize_to_genus(mc21_v8)

mc51_v4_g <- summarize_to_genus(mc51_v4)
mc51_v8_g <- summarize_to_genus(mc51_v8)

mc61_v4_g <- summarize_to_genus(mc61_v4)
mc61_v8_g <- summarize_to_genus(mc61_v8)

#################3

co1_hi_g <- CO1_g$Hi
co1_mi_g <- CO1_g$Mi
S16_hi_g <- S16_g$hi+1
S16_mi_g <- S16_g$mi+1
mc11_v4_hi_g <- mc11_v4_g$hi
mc11_v4_mi_g <- mc11_v4_g$mi
mc11_v8_hi_g <- mc11_v8_g$hi
mc11_v8_mi_g <- mc11_v8_g$mi

mc21_v4_hi_g <- mc21_v4_g$hi
mc21_v4_mi_g <- mc21_v4_g$mi
mc21_v8_hi_g <- mc21_v8_g$hi
mc21_v8_mi_g <- mc21_v8_g$mi

mc51_v4_hi_g <- mc51_v4_g$hi
mc51_v4_mi_g <- mc51_v4_g$mi
mc51_v8_hi_g <- mc51_v8_g$hi
mc51_v8_mi_g <- mc51_v8_g$mi

mc61_v4_hi_g <- mc61_v4_g$hi
mc61_v4_mi_g <- mc61_v4_g$mi
mc61_v8_hi_g <- mc61_v8_g$hi
mc61_v8_mi_g <- mc61_v8_g$mi


par(mfrow = c(5,2), oma = c(0,0,3,0))
compare_hi_mi(co1_hi_g, co1_mi_g, title = "CO1")
compare_hi_mi(S16_hi_g, S16_mi_g, log_ = T, title = "16S")
compare_hi_mi(mc11_v4_hi_g, mc11_v4_mi_g, title = "MC 11 V4")
compare_hi_mi(mc11_v8_hi_g, mc11_v8_mi_g, title = "MC 11 V8-9")

compare_hi_mi(mc21_v4_hi_g, mc21_v4_mi_g, title = "MC 21 V4")
compare_hi_mi(mc21_v8_hi_g, mc21_v8_mi_g, title = "MC 21 V8-9")

compare_hi_mi(mc51_v4_hi_g, mc51_v4_mi_g, title = "MC 51 V4")
compare_hi_mi(mc51_v8_hi_g, mc51_v8_mi_g, title = "MC 51 V8-9")

compare_hi_mi(mc61_v4_hi_g, mc61_v4_mi_g, title = "MC 61 V4")
compare_hi_mi(mc61_v8_hi_g, mc61_v8_mi_g, title = "MC 61 V8-9")
mtext(expression(bold("Taxonomy compared at the Genus level")), outer = T, cex = 1.25)
mtext(side = 1, "Relative Abundance in MiSeq length (%)", outer = TRUE, cex = 1.25, line = -1)
mtext(side = 2, "Relative Abundance in HiSeq length (%)", outer = TRUE, cex = 1.25, line = -2)

generate_table <- function(hi, mi, title) {
  fit <- lm(hi~mi)
  summ <- summary(fit)
  
  cf <- round(coef(fit), 5) 
  eq <- paste0("y = ", 
               cf[2], " x",
               ifelse(sign(cf[1])==1, " + ", " - "), abs(cf[1]))
  
  rse <- summ$sigma
  df <- summ$df[2]
  adj_R <- summ$adj.r.squared
  fstat <- summ$fstatistic[1]
  pval <- summ$coefficients[2,4]
 
  return(c(title, eq, rse,df, adj_R, fstat, pval)) 
}

a <- generate_table(co1_hi_g, co1_mi_g, title = "CO1 - genus")
b <- generate_table(S16_hi_g, S16_mi_g,  title = "16S - genus")
c <- generate_table(mc11_v4_hi_g, mc11_v4_mi_g, title = "MC 11 V4 - genus")
d <- generate_table(mc11_v8_hi_g, mc11_v8_mi_g, title = "MC 11 V8-9 - genus")

e <- generate_table(mc21_v4_hi_g, mc21_v4_mi_g, title = "MC 21 V4 - genus")
f <- generate_table(mc21_v8_hi_g, mc21_v8_mi_g, title = "MC 21 V8-9 - genus")

h <- generate_table(mc51_v4_hi_g, mc51_v4_mi_g, title = "MC 51 V4 - genus")
i <- generate_table(mc51_v8_hi_g, mc51_v8_mi_g, title = "MC 51 V8-9 - genus")

j <- generate_table(mc61_v4_hi_g, mc61_v4_mi_g, title = "MC 61 V4 - genus")
k <- generate_table(mc61_v8_hi_g, mc61_v8_mi_g, title = "MC 61 V8-9 - genus")


aa <- generate_table(co1_hi, co1_mi, title = "CO1 - species")
bb <- generate_table(S16_hi, S16_mi,  title = "16S - species")
cc <- generate_table(mc11_v4_hi, mc11_v4_mi, title = "MC 11 V4 - species")
dd <- generate_table(mc11_v8_hi, mc11_v8_mi, title = "MC 11 V8-9 - species")

ee <- generate_table(mc21_v4_hi, mc21_v4_mi, title = "MC 21 V4 - species")
ff <- generate_table(mc21_v8_hi, mc21_v8_mi, title = "MC 21 V8-9 - species")

hh <- generate_table(mc51_v4_hi, mc51_v4_mi, title = "MC 51 V4 - species")
ii <- generate_table(mc51_v8_hi, mc51_v8_mi, title = "MC 51 V8-9 - species")

jj <- generate_table(mc61_v4_hi, mc61_v4_mi, title = "MC 61 V4 - species")
kk <- generate_table(mc61_v8_hi, mc61_v8_mi, title = "MC 61 V8-9 - species")

to_save <- rbind(a,aa,b,bb,c,cc,d,dd,e,ee,f,ff,h,hh,i,ii,j,jj,k,kk)
colnames(to_save) <- c("Dataset", "linear regression eq", "Residual SE", "DF", "Adj. R Squared", "F statistic", "P Value")
write.csv(to_save, file = "figures-and-tables-for-the-Github/MiSeq-v-PseudoHiSeq-SupTable.csv")
