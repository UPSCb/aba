library(tidyverse)
library(stringr)

raw_data_dir <- "/mnt/picea/storage/data/aspseq/rbhalerao/ABA/RNA-Seq/IGA/rishi/ABAproject"
fq_files <- list.files(raw_data_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

data.frame(FileName = basename(fq_files), FileLocation = dirname(fq_files),
           stringsAsFactors = FALSE) %>% 
  mutate(SampleName = str_extract(FileName, "\\d+w-(abi1|t89)_\\d+"),
         time = as.integer(sub("^(\\d+)w.+$", "\\1", SampleName)),
         genotype = str_extract(SampleName, "abi1|t89"),
         genotype = ifelse(genotype == "abi1", "abi1 mutant", "T89"),
         SampleDescription = paste(genotype, "after", time, "weeks of short day treatment"),
         SequencingDate = "2014-01-01",
         ExperimentTitle = "Short day treatment of ABA insensitive hybrid aspen") %>% 
  select(ExperimentTitle, SampleName, SampleDescription, SequencingDate, FileName, FileLocation) %>% 
  write_csv("~niklasm/git/UPSCb/projects/aba/doc/ena_submission.csv")
