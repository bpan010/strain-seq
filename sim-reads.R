# load up jackalope
# read fasta file
# using read_fasta

# read vcf file
# using create_haplotypes



library(jackalope)
library(tidyverse)



ref_file <-  paste0("~/Documents/Stanford Courses/Lucas Lab/ncbi_dataset/",
       "ncbi_dataset/data/GCA_003401565.2/GCA_003401565.2_MR1_a10_genomic.fna")

ref <- read_fasta(ref_file, cut_names=FALSE)




new_names <- ref$chrom_names() |>
    str_split(" ") |>
    map_chr(\(x) x[6]) |>
    str_remove_all(",")
    

ref$set_names(new_names)



vcf_file <- paste0("~/Documents/Stanford Courses/Lucas Lab/strain-seq/",
        "Dhami_et_al_final_SNPset_noMT.vcf")

mutations <- create_haplotypes(ref, haps_vcf(vcf_file))

mutations
