# load up jackalope
# read fasta file
# using read_fasta

# read vcf file
# using create_haplotypes



library(jackalope)
library(tidyverse)



# strains that we want to work with
strains <- c("Y382", "Y383", "Y385", "Y466", "Y467", "Y644",
             "Y818", "Y821", "Y858", "Y893", "Y1092", "MR1")

# get the reference file
ref_file <-  paste0("~/Documents/Stanford Courses/Lucas Lab/ncbi_dataset/",
       "ncbi_dataset/data/GCA_003401565.2/GCA_003401565.2_MR1_a10_genomic.fna")

ref <- read_fasta(ref_file, cut_names=FALSE)


# str_replace("input string", "remove condition eg X.*Y", "replace")

# str_replace("Y989-TAGGCATG-TCTACTCT_S43_2", "-.*_", "_")
# str_replace("MR1_1", "-.*_", "_")

# make new names
# replace old names with new names
# remove any names in old haplotypes vector that we dont want


# reformat the reference file into the format we want
new_names <- ref$chrom_names() |>
    str_split(" ") |>
    map_chr(\(x) x[6]) |>
    str_remove_all(",")
    
ref$set_names(new_names)


# get haplotypes from vcf_file

vcf_file <- paste0("~/Documents/Stanford Courses/Lucas Lab/strain-seq/",
        "Dhami_et_al_final_SNPset_noMT.vcf")

haplotypes <- create_haplotypes(ref, haps_vcf(vcf_file))

# reformat the names into the format we want
new_hap_names <- haplotypes$hap_names() |>
    map_chr(\(x) str_replace(x, "-.*_", "_"))

haplotypes$set_names(new_hap_names)


# z <- letters[1:4]
# z[z %in% c("a", "b")]

# create temporary vector of all haplotype names
nobueno_haplotypes <- haplotypes$hap_names()

# add _1 and _2 to every element in strains to compare to our haplotype names
new_strains <- c(paste0(strains, "_1"), paste0(strains, "_2"))
    
# filter out the haplotypes we do want to leave only the ones to be removed
nobueno_haplotypes <- nobueno_haplotypes[!(nobueno_haplotypes %in% new_strains)]

# paste0 combines strings

# c(paste0("_1"), paste0("_2"))

# remove all the haplotypes we don't want
haplotypes$rm_haps(nobueno_haplotypes)

# to see which names are left. should be the same as new_strains
result <- haplotypes$hap_names()
