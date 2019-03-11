#!/usr/bin/python2
import sys
import re

file= sys.argv[1] 
outfile = sys.argv[2]


#file='/Users/limeybean/Desktop/test_blac.out.txt'
#tempfile='/Users/limeybean/Desktop/test_blac.out.temp.txt'
f= open(outfile,"w+")

# make a temp file that removes everything but taxonomic rank and taxon name

Rank_names = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
nodice = "Unclassified"

with open(file, "r") as my_file:
  for line in my_file:
  	  line = re.sub('Unclassified', 'Unclassified\\t', line)
  	  str = line.split('\t')
  	  un = str[1]
  	  if un == 'Unclassified':
  	  	print >> f, str[0], '\t', str[1]
  	  else:
 	  	rank = str[1]
  	  	exclusions = '|'.join(Rank_names)
  	  	rank = re.sub(exclusions, '', rank)
  	  	rank = re.sub('^:', '', rank)
  	  	rank = re.sub(';$', '', rank)
  	  	rank_id = rank.split(';:')
  	  	Rank_dict = {}
  	  	for i in range(len(Rank_names)):
  	  		Rank_dict[Rank_names[i]] = rank_id[i]
  	  	bcc = str[2]
  	  	exclusions = '|'.join(Rank_names)
  	  	bcc = re.sub(exclusions, '', bcc)
  	  	bcc = re.sub('^:', '', bcc)
  	  	bcc = re.sub(';$', '', bcc)
  	  	bcc_id = bcc.split(';:')
  	  	bcc_dict = {}
  	  	for i in range(len(Rank_names)):
  	  		bcc_dict[Rank_names[i]] = bcc_id[i]
  	  	print >> f, str[0], '\t', 'superkingdom:', Rank_dict['superkingdom'], ';', bcc_dict['superkingdom'], ';' , 'phylum:', Rank_dict['phylum'], ';', bcc_dict['phylum'], ';' , 'class:', Rank_dict['class'], ';', bcc_dict['class'], ';' , 'order:', Rank_dict['order'], ';', bcc_dict['order'], ';' , 'family:', Rank_dict['family'], ';', bcc_dict['family'], ';' , 'genus:', Rank_dict['genus'], ';', bcc_dict['genus'], ';' , 'species:', Rank_dict['species'], ';', bcc_dict['species'], ';'
       
#REF_GI_265678402	superkingdom:Bacteria;100.0;phylum:Actinobacteria;100.0;class:Actinobacteria;100.0;order:Micromonosporales;66.2638888889;family:Micromonosporaceae;66.2638888889;genus:Micromonospora;22.2078671329;species:Catenulispora_rubra;18.1111111111;



f.close()



