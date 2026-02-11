#Script Title: Concatipede Testing
#Authors: Marcus Barela, Ariel Goldsztejn,Matteo Vecchi and Matthieau Bruneaux
#Date:4.XI.2025


install.packages("concatipede")
install.packages("tidyverse")
install.packages("devtools")
devtools::install_github("tardipede/concatipede")

library(concatipede)
library(tidyverse)
library(devtools)
library(DT)

# Save the path to the initial directory for later clean-up

setwd("G:/Shared drives/bioinf")

old_dir <- getwd()

# Create a directory to put the fasta files for this example if you dont have one

dir.create("concatipede_test")

# Set it as the working directory

setwd("concatipede_test")

#THIS SECTION IS FOR USE WITH THE PACKAGE EXAMPLE FILES ONLY, IF YOU PLAN TO USE YOUR OWN EXCEL SHEETS OR FASTA FILES, DO NOT USE THIS SECTION OF CODE
# Copy the example fasta files shipped with the package into that directory


example_files = list.files(system.file("extdata", package="concatipede"), full.names = TRUE)

file.copy(from = example_files, to = getwd())

#END EXAMPLE CODE SECTION HERE, FROM HERE ON OUT YOU CAN USE AS YOU SEE FIT

#DO NOT NEED THIS SECTION IF AUTOMATCHING BINS BY SEQUENCE NAMES- this is the same as using the sequence names as the reference for the bins you will create fo concatenation in the next step
find_fasta()

concatipede_prepare(out = "seqnames")
#END OF OPTIONAL CODE 

#USE THIS FOR IF YOU ARE AUTOMATCHING AS EXPLAINED PREVIOUSLY

find_fasta() %>% #finds fasta files in the current dir
  concatipede_prepare() %>%                  #prepares the concatenation file
  
  write_xl("template.xlsx") %>%              #sets the excel file to be used to create the finished binned excel file
  
  auto_match_seqs() %>%                      #automatches the sequences in the input excel to create bins based on the genbank names in the cells for the individuals
  
  write_xl("template_automatched.xlsx") %>%  #sets the output excel file
  
  concatipede() %>%                          #does the actual concatenation 
  
  write_fasta("merged-seqs.fasta") %>%       #creates a new FASTA
  image(cex=0.3)                             #creates the concatenation as an image in the plot area



#Use this if not automatching or if you want to just continue on with your current dataset after automatching

concatipede(filename = "template_automatched.xlsx", out = "Macrobiotidae_4genes", excel.sheet = 1)

old_par <- par(mar=c(1,12,1,1))
image(concatipede(filename = "template_automatched.xlsx",out="Macrobiotidae_4genes",excel.sheet = 1),cex=0.5)

#HERE make a second sheet in the workbook you are currently using to create he concatenation that has just the names and only one locus

concatipede(filename = "template_automatched.xlsx",out="Macrobiotus_4genes",excel.sheet = 2) #the sheet number listed will be the sheet number in the workbook used to make the concatentation

par(mar=c(1,12,1,1))
image(concatipede(filename = "template_automatched.xlsx",out="Macrobiotus_4genes",excel.sheet = 2),cex=0.5)

#same thing but for a sheet #3

concatipede(filename = "template_automatched.xlsx",out="Macrobiotidae_COI",excel.sheet = 3)

par(mar=c(1,12,1,1))
image(concatipede(filename = "template_automatched.xlsx",out="Macrobiotidae_ITS2",excel.sheet = 3),cex=0.5)

#this is where your genbank accession table will be created from the accession  numbers and tags from your correspondence table
genbank.table = get_genbank_table(filename = "template_automatched.xlsx", excel.sheet = 1)

#prints the table, makes it interactive with arguments

get_genbank_table(filename = "template_automatched.xlsx", excel.sheet = 1) %>% datatable(extensions = 'Buttons',
                                                                                           options = list(dom = 'Blfrtip',
                                                                                                          rownames = F,
                                                                                                          buttons = c('excel'),
                                                                                                          scrollX=TRUE,
                                                                                                          lengthMenu = list(c(nrow(.),25,50,-1),
                                                                                                                            c(nrow(.),25,50,"All"))))
#rename the sequence marker names for simplicity or whatever you want

rename_sequences(filename = "template_automatched.xlsx", excel.sheet = 1, marker_names = c("COI","ITS2","LSU","SSU"))

#back to the main wd

setwd(old_dir)

# Invisible clean-up
par(old_par)
options(old_options)

