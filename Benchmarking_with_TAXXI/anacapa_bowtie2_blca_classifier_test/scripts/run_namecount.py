#!/usr/bin/python2

# script take a taxonomy file and counts the occurrences of taxonomic rank
# usage e.g.:
# python2 ~/taxxi/py/run_namecount.py ~/Downloads/test.txt ~/taxxi/blca_recreate/namecounts test

import re
import string 
import sys
import os


file_to_namecount= sys.argv[1] 
outdir = sys.argv[2]
name = sys.argv[3]

tempfile = outdir + "/temp.txt" 


file=file_to_namecount
f= open(tempfile,"w+")

# make a temp file that removes everything but taxonomic rank and taxon name

with open(file, "r") as my_file:
  for line in my_file:
      str = line.split('tax=')[-1]
      newstr = str.replace(";\n", ",")
      newstr1 = newstr.replace(";", ",")
      f.write(newstr1.replace(",", "\n"))

f.close()

# count occurrences of taxonomic rank and taxon name printing it to file


f= open(outdir + "/" + name + "_namecount","w+")

frequency = {}
document_text = open(tempfile, 'r')
text_string = document_text.read()
match_pattern = re.findall('[a-z]?:?[A-Z]?[^\s]+', text_string)


 
for word in match_pattern:
    count = frequency.get(word,0)
    frequency[word] = count + 1
    
    
     
frequency_list = frequency.keys()

 
 
for words in frequency_list:
    print >> f, frequency[words],'\t', words
    
    
document_text.close()
f.close()

os.remove(tempfile)