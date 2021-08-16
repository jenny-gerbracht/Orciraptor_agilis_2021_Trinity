library(dplyr)
library(stringr)

#Simplify tip labels and add attributes for coloring in figree
setwd("~/transfer/Jenny/Paper_Trinity/Module_8/Input")
randomnames <- read.delim("mafft_input_randomNames.txt")

setwd("/srv/Jenny/EukProt/")
eukprot_anno <- read.delim("EukProt_included_data_sets.v02.2020_06_30.txt",
                           fileEncoding = "UTF16LE")


orc <- filter(randomnames, str_detect(Original.seqname, "Orciraptor"))
euk <- filter(randomnames, str_detect(Original.seqname, "^\\.\\/"))
out <- randomnames[!randomnames$Original.seqname %in% c(euk$Original.seqname, orc$Original.seqname),]


orc$attribute <- c("Orciraptor")


euk$EukProt_ID <- str_match(euk$Original.seqname, "EP\\d{5}")
eukprot_anno <- select(eukprot_anno, EukProt_ID, Supergroup_UniEuk)
euk %>% 
  left_join(eukprot_anno, by = "EukProt_ID") -> euk
euk$Original.seqname<- paste0(str_split(euk$Original.seqname, "\\_", simplify = TRUE)[,2],
                              "_",
                              str_split(euk$Original.seqname, "\\_", simplify = TRUE)[,3])
euk %>% 
  select(c(1, 2, 4)) %>% 
  setNames(., c("Assigned.random.name", "Original.seqname", "attribute")) -> euk


out$Original.seqname <- str_match(out$Original.seqname, "\\[(.*?)\\]")[,2]
out$attribute <- c("Bacteria")
out[c(101:105),]$attribute <- c("Outgroup")

randonnames_edit <- rbind(orc, euk, out)
randonnames_edit$Original.seqname <- make.unique(randonnames_edit$Original.seqname, sep = "_")
colnames(randonnames_edit) <- c("Assigned random name", "Original seqname", "attribute")

setwd("~/transfer/Jenny/Paper_Trinity/Module_8/Input")
write.table(select(randonnames_edit, c(1, 2)),
            file = "mafft_input_randomNames_edited.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

write.table(select(randonnames_edit, c(2, 3)),
            file = "annotation.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

