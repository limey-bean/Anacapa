#Use the function pairwise.adonis() with following arguments
#
#x = community table
#
#factors = a column or vector with all factors to be tested pairwise
#
#sim.function = which function to calculate the similarity matrix. eg 'daisy' or 'vegdist' default is 'vegdist'
# NOTE that if you wnat to use daisy, you need to install library 'cluster'
#
#sim.method = similarity method from daisy or vegdist: default is 'bray' 
#
#p.adjust.m = the p.value correction method, one of the methods supported by p.adjust(); default is 'bonferroni'
#
#The function will return a table with the pairwise factors, F-values, R^2, p.value and adjusted p.value
#
# load in runnig R session with: 
# source('pairwise.adonis.txt')
#

# example:
# data(iris)
# pairwise.adonis(iris[,1:4],iris$Species)
#
#[1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
#                    pairs    F.Model        R2 p.value p.adjusted sig
#1    setosa vs versicolor  552.51266 0.8493496   0.001      0.003   *
#2     setosa vs virginica 1001.54509 0.9108722   0.001      0.003   *
#3 versicolor vs virginica   91.82959 0.4837475   0.001      0.003   *

# similarity euclidean from vegdist and holm correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')

# similarity manhattan from daisy and bonferroni correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='daisy',sim.method='manhattan',p.adjust.m='bonferroni')
##############################


##start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

## end copy here
