---
title: "R Notebook"
output:
  pdf_document: default
  html_document:
    df_print: paged
---

---
title: "R Notebook"
output:
  html_document:
    df_print: paged
---

---
title: "R_diversity_analysis"
author: "Alexander W. Gofton"
date: "26/07/2021"
output: html_document
---

## Install and load packages
```{r}
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggthemes)
library(ggpubr)
library(splitstackshape)
#library(devtools)
#devtools::install_github("microbiome/microbiome")
library(microbiome)
#devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
#devtools::install_github("MadsAlbertsen/ampvis2")
library(ampvis2)
library(RColorBrewer)
```
Check working dir
```{r}
getwd()
```

### Calculate tpm from mapping stats for each contig
```{r}
file_list <- list.files(path="tpm_files", full.names = TRUE, pattern = "*.txt")

for (i in 1:length(file_list)){
  # Read in the bt2 mapping stats file and take cols 1-0
  tpm <- read.csv(file_list[i], header=T, sep="\t") %>% select(X.rname:meanmapq)
  # Get length of contig in Kb
  tpm$kb <- tpm[,3]/1000
  # Divide num_reads / contig length(kb) to give RPK
  tpm$RPK <- tpm[,4]/tpm[,10]
  # Sum all RPK and divide by 1 mil to give tpm scaling factor
  tpm_scaling_factor <- sum(tpm$RPK)/1000000
  # Divide RPK by scaling factor
  tpm$TPM <- tpm[,11]/tpm_scaling_factor
  # Get relivant cols and write file
  tpm <- select(tpm, c(X.rname:meanmapq, TPM))
  # Save file
  sampleID <- read.table(text = file_list[i], sep = ".", as.is = TRUE)$V1
  outfile <- paste(SampleID, '_tpm.csv', sep="")
  #write.csv(tpm, file=outfile, quote=FALSE, row.names=FALSE)
}

# Files are saved so now clear environment
rm(list=ls())
gc()
```

### Join read to class (from MEGAN) files and tpm files
```{r}
# list have corresponding files
tpm_files <- list.files(path = "tpm_files", full.names = TRUE, pattern = "*.csv")
r2c_files <- list.files(path = "r2c_files", full.names = TRUE, pattern = "*.txt")
out_files <- list.files(path = "r2c_files", pattern = "*.txt")

for (i in 1:length(r2c_files)){
  # Import files and format cols
  tpm <- read.csv(tpm_files[i], header = TRUE) %>%
         select(X.rname, TPM) %>% 
         rename(Contig_ID = X.rname)
  r2c <- read.csv(r2c_files[i], header = FALSE, sep = "\t", 
                  col.names = c("Contig_ID", "Level", "Tax_path")) %>%
         select(Contig_ID, Tax_path)
  # Join r2c and TPM files, keeping only contigs with an r2c taxonomy path
  r2c_tpm <- left_join(r2c, tpm, by="Contig_ID")
  # Separate tax path into cols and rename cols
  r2c_tpm <- cSplit(r2c_tpm, "Tax_path", sep = "; ", stripWhite = TRUE) %>%
             rename(Domain = Tax_path_1, Kingdom = Tax_path_2, Phylum = Tax_path_3,
                    Class = Tax_path_4, Order = Tax_path_5, Family = Tax_path_6, 
                    Genus = Tax_path_7, Species = Tax_path_8)
  # Get sample ID from file name (out_files)
  SampleID <- read.table(text = out_files[i], sep = ".", as.is = TRUE)$V1
  # Group and summarize by identical taxonomy, summing TPM
  col_r2c_tpm <- r2c_tpm %>% 
                 group_by(Domain, Kingdom, Phylum, 
                          Class, Order, Family, 
                          Genus, Species) %>%
                 summarize(Sum_TPM = sum(TPM)) %>%
                 ungroup()
  # Rename TPM col with sampleID
  names(col_r2c_tpm)[names(col_r2c_tpm) == "Sum_TPM"] <- SampleID
  # Create new character string to be used as output file name
  r2c_tpm_sampleID <- paste('r2c_tpm_files/', SampleID, '_r2c_tpm.csv', sep="")
  # Save file
  #write.csv(col_r2c_tpm, file = r2c_tpm_sampleID, quote = FALSE, row.names = FALSE)
}

# Files are saved so clear the environment
rm(list=ls())
gc()
```
### Join all r2c_tpm files together to creat an OTU + taxonomy file
```{r}
## Cols to merge by
#tax_cols <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## List files
#r2c_tpm_files <- list.files(path = "r2c_tpm_files", full.names = TRUE, pattern = "*.csv")
#
## Initial join
#df1 <- read.csv(r2c_tpm_files[1])
#df2 <- read.csv(r2c_tpm_files[2])
#merged_r2c_tpm <- full_join(df1, df2, by = tax_cols)
#rm(df1, df2)
#
## Loop through other files, adding one by one, so that memory isn't exceeded 
#for (i in 3:152){
#  dfn <- read.csv(r2c_tpm_files[i])
#  merged_r2c_tpm <- full_join(merged_r2c_tpm, dfn, by = tax_cols)
#  rm(dfn)
#}
#
## Convert NA values in count cols with 0s
#merged_r2c_tpm <- replace_na(merged_r2c_tpm, list(AR1=0, AR2=0, AR3=0, AR4=0, AR5=0, AR6=0, AR7=0, BH1=0, K2=0, K3=0, K4=0, K7=0, K10=0, K11=0, K12=0, K16=0, K18=0, K19=0, K20=0, K21=0, K22=0, K23=0, K24=0, #K26=0, K27=0, K29=0, K30=0, K31=0, K32=0, K33=0, K34=0, K35=0, K36=0, K37=0, K39=0, K40=0, K41=0, K42=0, K43=0, K44=0, K45=0, K46=0, K47=0, K48=0, K63=0, K64=0, K65=0, K66=0, K67=0, K68=0, K69=0, K71=0, K72=0, #K73=0, K74=0, K75=0, K76=0, K77=0, K78=0, K80=0, K82=0, K84=0, K85=0, K95=0, K96=0, K99=0, K100=0, K101=0, K102=0, K103=0, K109=0, K110=0, K117=0, K118=0, K119=0, K120=0, K122=0, K126=0, K127=0, K128=0, K130=0, #LR1=0, LR2=0, LR3=0, LR4=0, PR1=0, PR2=0, PR3=0, PR5=0, PR9=0, PR10=0, WW4=0, WW5=0, WW6=0, WW7=0, WW8=0, WW9=0, WW11=0, WW12=0, WW13=0, WW14=0, WW15=0, WW16=0, WW17=0, WW18=0, WW21=0, WW22=0, PR426_B=0, #PR428_B=0, PR429_B=0, PR430_B=0, PR431_B=0, PR432_B=0, PR433_B=0, PR434_B=0, PR435_B=0, PR436_B=0, PR437_B=0, PR438_B=0, PR439_B=0, PR441_B=0, LR2_B=0, LR6_B=0, LR7_B=0, LR8_B=0, LR9_B=0, LR10_B=0, LR11_B=0, #LR12_B=0, LR13_B=0, LR14_B=0, CR4_B=0, CR7_B=0, CR8_B=0, CR9_B=0, M1_B=0, M2_B=0, M4_B=0, P2_B=0, P3_B=0, P4_B=0, P6_B=0, R3a_B=0, R3b_B=0, SD1_B=0, SD2_B=0, SD3_B=0, SD4_B=0, SD5_B=0, SD6_B=0, SD7_B=0, #KP1_B=0))
#
## Save file
##write.csv(merged_r2c_tpm, file = "r2c_tpm_merged.csv", row.names = FALSE)
```
"r2c_tpm_merged.csv" was then edited in excel to properly format and curate taxonomy names.

### Importing OTU+taxonomy table and metadata table for ampzis2
```{r}
# Import file & collapse by taxonomy, rounding up because phyloseq and ampsiz2 expect integers
#av2_OTU_tab <- read.csv("r2c_tpm_merged_formatted.csv", header = TRUE) %>%
#               group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
#               summarize(across(.cols = AR1:WW9, sum)) %>%
#               relocate(Kingdom:Species, .after = WW9) %>%
#               ungroup() %>% 
#               mutate_at(1:152, ceiling)
## Count number of cols
#nrows <- nrow(av2_OTU_tab)
## Make OTUn sequence
#prefix <- "OTU"
#suffix <- as.vector(seq(1:nrows))  # make sure this number=num rows
#OTU_IDs <- paste(prefix, suffix, sep="_")
## Add OTUn column
#av2_OTU_tab <- mutate(av2_OTU_tab, OTU = OTU_IDs, .before = AR1)
## Write file
#write.csv(av2_OTU_tab, file = "av2_OTU_tax_table.csv", quote = FALSE, row.names = FALSE)
```
######################################################################
Import av2 OTU table
```{r}
av2_OTU_tab <- read.csv("av2_OTU_tax_table2.csv", header = TRUE)
```

Import av2 metadata file
```{r}
av2_metadata <- read.csv("metadata.csv", header = TRUE)
```

### Create phyloseq-formatted OTU table, taxonomy table, and metadata table
```{r}
# Metadata is the same as for av2
ps_metadata <- read.csv("metadata.csv", header = TRUE)
# Phyloseq required tax and otu counts to be in separate files, with OTU_ID in the first col of each.
# We can select these columns from the av2 OTU table and save them as separate dfs

ps_otu <- select(av2_OTU_tab, OTU, AR1:WW9)
write_csv(ps_otu, "ps_OTU_table.csv")
ps_tax <- select(av2_OTU_tab, OTU, Kingdom:Species) %>% 
  mutate_if(is.character, str_replace_all, pattern = "k__", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "p__", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "c__", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "o__", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "f__", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "g__", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "s__", replacement = "")
write_csv(ps_tax, "ps_tax_table.csv")
```

### Create ampviz2 object
```{r}
av2 <- amp_load(otutable=av2_OTU_tab, metadata=av2_metadata)
```

Create ampviz2 object and filter unassigned and host contigs from ampziv2 object
```{r}
av2 <- amp_load(otutable=av2_OTU_tab, metadata=av2_metadata)
unass_tax <- c("k__Unassigned", "p__Unassigned", "c__Unassigned", "o__Unassigned", "f__Unassigned","g__Unassigned", 
               "s__Unassigned", "k__unclassified Bacteria",	"p__unclassified Bacteria",	"c__unclassified Bacteria",	
               "o__unclassified Bacteria",	"f__unclassified Bacteria",	"g__unclassified Bacteria",	"s__unclassified Bacteria", 
               "k__unclassified Eukaryota",	"p__unclassified Eukaryota",	"c__unclassified Eukaryota",	"o__unclassified Eukaryota", 
               "f__unclassified Eukaryota",	"g__unclassified Eukaryota",	"s__unclassified Eukaryota")
av2_unass_filt <- amp_subset_taxa(av2, tax_vector = unass_tax, normalise = FALSE, remove = TRUE)
host_tax <- c("c__Mammalia", "f__Ixodidae")
av2_filt <- amp_subset_taxa(av2_unass_filt, tax_vector = host_tax, normalise = FALSE, remove = TRUE)
```

### Create phyloseq object
```{r}
ps_otus <- read.csv("ps_OTU_table.csv", row.names = "OTU")
ps_otu_mat <- as.matrix(ps_otus)
ps_otu_mat <- otu_table(ps_otu_mat, taxa_are_rows = TRUE)
# Tax file headers needs some formatting for phyloseq
ps_tax <- read.csv("ps_tax_table.csv", row.names = "OTU") %>% 
          rename(Domain = Kingdom)
# Convert to matricies
ps_tax_mat <- as.matrix(ps_tax)
ps_metadata <- read.csv("metadata.csv", row.names = "Sample_ID")
ps_metadata = sample_data(data.frame(ps_metadata))
# Import as ps objects
psOTU <- otu_table(ps_otu_mat, taxa_are_rows = TRUE)
psTAX <- tax_table(ps_tax_mat)
ps = phyloseq(psOTU, psTAX, ps_metadata)
```

Filtering unassigned and host taxonomies from phyloseq object
```{r}
# Filter out unassigned and host contigs
ps_filt <- subset_taxa(ps, Domain != "k__Unassigned")
ps_filt <- subset_taxa(ps_filt, Class != "c__Mammalia")
ps_filt <- subset_taxa(ps_filt, Family != "f__Ixodidae")
ps_filt <- subset_taxa(ps_filt, Domain != "k__unclassified Eukaryota")
ps_filt <- subset_taxa(ps_filt, Domain != "k__unclassified Bacteria")

# Species == I. holocyclus only
ps_filt_Ihol <- subset_samples(ps_filt, Sample_species == "Ixodes holocyclus")
# Species == H. bancrofti only
ps_filt_Hban <- subset_samples(ps_filt, Sample_species == "Haemaphysalis bancrofti")
ps_filt_Hban_Kioloa <- subset_samples(ps_filt_Hban, Site == "Kioloa")
# Ticks only
ps_filt_ticks <- subset_samples(ps_filt, Type == "Tick")
# Bloods only
ps_filt_bloods <- subset_samples(ps_filt, Type == "Blood")
```

Ticks alpha rarefaction with ampzis 2
```{r}
av2_ticks_filt <- amp_subset_samples(av2_filt, Type == "Tick", removeAbsents = TRUE)
arare_ticks_filt <- amp_rarecurve(av2_ticks_filt, stepsize = 100, color_by = "Lifestage", 
                                  facet_by = c("Site", "Sample_species"), facet_scales = "free") + 
  ylab(label = "Taxonomic units") + 
  xlab(label = "Rarefaction depth (number of subsampled contigs)") + 
  theme_bw() + 
  scale_color_brewer(palette = "Paired")

arare_ticks_filt
ggsave("arare_ticks_filt.pdf", path="./plots/Final_figs", width=10, height=10)
```

Blood samples alpha rarefaction (on host + unassigned filtered data)
```{r}
av2_bloods_filt <- amp_subset_samples(av2_filt, Type == "Blood", removeAbsents = TRUE)
arare_bloods_filt <- amp_rarecurve(av2_bloods_filt, stepsize = 100, color_by = "Sample_species",
                                   facet_by = "Site", facet_scales = "free") + 
  ylab(label = "Taxonomic units") + 
  xlab(label = "Rarefaction depth (number of subsampled contigs)") + 
  theme_bw() + 
  scale_color_brewer(palette = "Paired")

arare_bloods_filt
ggsave("arare_bloods_filt.pdf", path = "./plots/alpha_rarefaction", width = 10, height=10)
```
Arrange tick and blood alpha rarefaction together
```{r}
x <- ggarrange(arare_ticks_filt, arare_bloods_filt, nrow = 2, ncol = 1, labels = c("A", "B"))
x

ggsave("Figure S1 Alpha rarefaction.pdf", path = "./plots/Final_figs", height = 22, width = 17, units = "cm")
```
```{r}
distrib <- plot_read_distribution(ps_filt, groups = "Sample_species", plot.type = "density") + 
  xlab("TPM-normalised contig abundance") + ylab("Density")


distrib

```mbh










### Top taxa boxplots before and after filtering hosts
```{r}
# Ticks with host contigs
av2_ticks <- amp_subset_samples(av2_unass_filt, Type == "Tick", removeAbsents = TRUE)
ticks_top_tax_bp <- amp_boxplot(av2_ticks, group_by = "Type", sort_by = "median", plot_type = "boxplot", 
                                point_size = 0.8, tax_aggregate = "Phylum", tax_add = "Kingdom", 
                                tax_show = 10, tax_empty = "remove", normalise = TRUE, plot_flip = TRUE) + 
  theme_bw() + 
  scale_color_manual(values = c("grey25")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.95, size = 8), axis.title.x = element_blank()) + 
  theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))

ticks_top_tax_bp
# Ticks without host contigs
av2_ticks <- amp_subset_samples(av2_filt, Type == "Tick", removeAbsents = TRUE)
ticks_top_tax_bp_nohost <- amp_boxplot(av2_ticks, group_by = "Type", sort_by = "median", plot_type = "boxplot", 
                                point_size = 0.8, tax_aggregate = "Phylum", tax_add = "Kingdom", 
                                tax_show = 10, tax_empty = "remove", normalise = TRUE, plot_flip = TRUE) + 
  theme_bw() + 
  scale_color_manual(values = c("grey25")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.95, size = 8)) + 
  theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))

ticks_top_tax_bp_nohost
# Bloods with host contigs
av2_blood <- amp_subset_samples(av2_unass_filt, Type == "Blood", removeAbsents = TRUE)
blood_top_tax_bp <- amp_boxplot(av2_blood, group_by = "Type", sort_by = "median", plot_type = "boxplot", 
                                point_size = 0.8, tax_aggregate = "Phylum", tax_add = "Kingdom", tax_show = 10, 
                                tax_empty = "remove", normalise = TRUE,plot_flip = TRUE) + 
  theme_bw() + 
  scale_color_manual(values = c("grey25")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.95, size = 8), 
        axis.title.x = element_blank()) +
  theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))

blood_top_tax_bp
# Bloods without host contigs
av2_blood <- amp_subset_samples(av2_filt, Type == "Blood", removeAbsents = TRUE)
blood_top_tax_bp_nohost <- amp_boxplot(av2_blood, group_by = "Type", sort_by = "median", plot_type = "boxplot", 
                                point_size = 0.8, tax_aggregate = "Phylum", tax_add = "Kingdom", tax_show = 10, 
                                tax_empty = "remove", normalise = TRUE,plot_flip = TRUE) + 
  theme_bw() + 
  scale_color_manual(values = c("grey25")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.95, size = 8), axis.title.x = element_blank()) +
  theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))

blood_top_tax_bp_nohost
# Putting the plots together
together <- ggarrange(ticks_top_tax_bp, ticks_top_tax_bp_nohost, 
                      blood_top_tax_bp, blood_top_tax_bp_nohost,
                      nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
together
ggsave("Top10_phyla_host_not_filtered.pdf", path = "./plots", height = 20, width = 25, units = "cm")
```

##Taxonomy heat maps
Top phyla with and without host contigs (as above)
```{r}
# all with host
av2_hm_all_unass_filt_txt <- amp_heatmap(av2_unass_filt, facet_by = "Type", group_by = "Type", 
                                           tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                                            normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 1, 
                                            plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                                            plot_legendbreaks = c(0.1, 1, 10), plot_na = FALSE, textmap = TRUE)
write.csv(av2_hm_all_unass_filt_txt, file = "av2_hm_all_unass_filt.csv")
# all without host
av2_hm_all_host_filt_txt <- amp_heatmap(av2_filt, facet_by = "Type", group_by = "Type", 
                                           tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                                            normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 1, 
                                            plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                                            plot_legendbreaks = c(0.1, 1, 10), plot_na = FALSE, textmap = TRUE)
write.csv(av2_hm_all_host_filt_txt, file = "av2_hm_all_host_filt.csv")

```
all heatmap (family-level)
```{r}
hm_rel_txt <- amp_heatmap(av2_filt, facet_by = "Sample_species", group_by = "Site", 
                                tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 100, 
                                normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 1, 
                                plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                                plot_legendbreaks = c(0.1, 1, 10), plot_na = FALSE, textmap = TRUE)
write.csv(hm_rel_txt, file = "hm_all_rel_txt_ticks.csv")

hm_all_rel <- amp_heatmap(av2_filt, facet_by = "Sample_species", group_by = "Site", 
                            tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 50, 
                            normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                            plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                            plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(strip.text.x = element_text(size = 10, face = "italic")) + 
  scale_y_discrete(labels = c("Remaining taxa (444)", "Bacteria; Bacteroidetes; Sphingobacteriaceae", "Bacteria; Actinobacteria; Kineosporiaceae", "Bacteria; Proteobacteria; Rhodanobacteraceae", "Bacteria; Proteobacteria; Acetobacteraceae", "Bacteria; Proteobacteria; Beijerinckiaceae", "Bacteria; Proteobacteria; Betaproteobacteria", "Bacteria; Bacteroidetes; Cytophagaceae", "Bacteria; Proteobacteria; Enterobacterales", "Bacteria; Proteobacteria; Bradyrhizobiaceae", "Bacteria; Actinobacteria; Corynebacteriales", "Bacteria; Actinobacteria; Streptomycetaceae", "Bacteria; Bacteroidetes; Hymenobacteraceae", "Eukaryota; Nematoda; Angiostrongylidae", "Viruses; Anelloviridae", "Bacteria; Bacteroidetes; Muribaculaceae", "Bacteria; Proteobacteria; Gammaproteobacteria", "Pararvirae; Artverviricota; Retroviridae", "Bacteria; Actinobacteria; Nocardiaceae", "Bacteria; Proteobacteria; Rhizobiaceae", "Bacteria; Proteobacteria; Proteobacteria", "Bacteria; Actinobacteria; Microbacteriaceae", "Bacteria; Actinobacteria; Propionibacteriaceae", "Bacteria; Proteobacteria; Enterobacteriaceae", "Bacteria; Proteobacteria; Burkholderiales", "Bacteria; Actinobacteria; Pseudonocardiaceae", "Bacteria; Proteobacteria; Yersiniaceae", "Bacteria; Proteobacteria; Alcaligeceae", "Bacteria; Actinobacteria; Nocardioidaceae", "Bacteria; Proteobacteria; Coxiellaceae", "Bacteria; Proteobacteria; Oxalobacteraceae", "Bacteria; Proteobacteria; Alphaproteobacteria", "Bacteria; Proteobacteria; Rhizobiales", "Bacteria; Proteobacteria; Methylobacteriaceae", "Eukaryota; Apicomplexa; Piroplasmida", "Eukaryota; Apicomplexa; Babesiidae", "Bacteria; Proteobacteria; Rickettsiaceae", "Bacteria; Actinobacteria; Actinobacteria", "Bacteria; Proteobacteria; Sphingomodaceae", "Bacteria; Proteobacteria; Bartonellaceae", "Eukaryota; Apicomplexa; Hepatozoidae", "Bacteria; Proteobacteria; Anaplasmataceae", "Bacteria; Proteobacteria; Comamodaceae", "Bacteria; Proteobacteria; Pseudomodaceae", "Bacteria; Actinobacteria; Mycobacteriaceae", "Bacteria; Proteobacteria; Burkholderiaceae", "Eukaryota; Euglenozoa; Trypanosomatidae", "Bacteria; Proteobacteria; Francisellaceae", "Bacteria; Tenericutes; Mycoplasmataceae", "Eukaryota; Apicomplexa; Theileriidae", "Bacteria; Proteobacteria; Candidatus Midichloriaceae"))

hm_all_rel
ggsave("Figure 3 Heatmap.pdf", path="./plots/Final_figs", width = 43, height = 20, units = "cm")
```
all heatt map (no site grouping)
```{r}
hm_rel_nogroup_txt <- amp_heatmap(av2_filt, facet_by = "Sample_species",
                                tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 3000, 
                                normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 1, 
                                plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                                plot_legendbreaks = c(0.1, 1, 10), plot_na = FALSE, textmap = TRUE)
write.csv(hm_rel_nogroup_txt, file = "hm_all_rel_txt_ticks.csv")

hm_all_nogroup_rel <- amp_heatmap(av2_filt, facet_by = "Sample_species", group_by = "Sample_species", 
                            tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 50, 
                            normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                            plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                            plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(strip.text.x = element_text(size = 10, face = "italic"))
  #scale_y_discrete(labels = c())

hm_all_nogroup_rel
ggsave("All_heatmap_nogroup.pdf", path="./plots/heatmaps", width = 43, height = 20, units = "cm")
```
Ticks heat map
```{r}
hm_ticks_rel_txt <- amp_heatmap(av2_ticks_filt, facet_by = "Sample_species", group_by = "Site", 
                                tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                                normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 1, 
                                plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                                plot_legendbreaks = c(0.1, 1, 10), plot_na = FALSE, textmap = TRUE)
write.csv(hm_ticks_rel_txt, file = "hm_rel_txt_ticks.csv")

hm_ticks_rel <- amp_heatmap(av2_ticks_filt, facet_by = "Sample_species", group_by = "Site", 
                            tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                            normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                            plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                            plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(strip.text.x = element_text(size = 10, face = "italic"))
  #scale_y_discrete(labels = c())

hm_ticks_rel
ggsave("Tick_heatmap.pdf", path="./plots/heatmaps", width = 30, height = 30, units = "cm")
```

Bloods (all samples) heat map
```{r}
hm_bloods_rel_txt <- amp_heatmap(av2_bloods_filt, facet_by = "Sample_species", group_by = "Site", 
                                 tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                                 normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 1, 
                                 plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                                 plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE, textmap = TRUE)
write.csv(hm_bloods_rel_txt, file = "hm_rel_txt_bloods.csv")

hm_bloods_rel <- amp_heatmap(av2_bloods_filt, facet_by = "Sample_species", group_by = "Site", 
                             tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                             normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                             plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                             plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(strip.text.x = element_text(size = 10, face = "italic"))
  #scale_y_discrete(labels = c())

hm_bloods_rel
ggsave("Blood_heatmap.pdf", path="./plots/Final_figs", width = 30, height = 30, units = "cm")
```

Ixodes holocyclus heatmap (Lifestage x site)
Family
```{r}
av2_Ihol <- amp_subset_samples(av2_filt, Sample_species == "Ixodes holocyclus", removeAbsents = TRUE)

hm_Ihol_rel_txt <- amp_heatmap(av2_Ihol, facet_by = "Site", group_by = "Lifestage", 
                               tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 100, 
                               normalise = TRUE, plot_values = TRUE, plot_values_size = 2, round = 1, 
                               plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                               plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE, textmap = TRUE)
write.csv(hm_Ihol_rel_txt, file = "hm_rel_txt_Ihol.csv")

hm_Ihol_rel <- amp_heatmap(av2_Ihol, facet_by = "Site", group_by = "Lifestage", 
                           tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 50, 
                           normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                           plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                           plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(strip.text.x = element_text(size = 10)) + 
  scale_y_discrete(labels = c("Remaining taxa (366)", "Bacteria; Actinobacteria; Gordoniaceae", "Bacteria; Bacteroidetes; Bacteroidetes", "Bacteria; Bacteroidetes; Flavobacteriaceae", "Eukaryota; Nematoda; Rhabditida", "Bacteria; Actinobacteria; Micrococcales", "Bacteria; Actinobacteria; Geodermatophilaceae", "Bacteria; Proteobacteria; Erwiniaceae", "Bacteria; Proteobacteria; Bradyrhizobiaceae", "Bacteria; Proteobacteria; Enterobacteriaceae", "Bacteria; Bacteroidetes; Cytophagaceae", "Bacteria; Proteobacteria; Alcaligeceae", "Bacteria; Proteobacteria; Beijerinckiaceae", "Bacteria; Proteobacteria; Xanthomodaceae", "Bacteria; Bacteroidetes; Sphingobacteriaceae", "Bacteria; Proteobacteria; Acetobacteraceae", "Bacteria; Actinobacteria; Nocardiaceae", "Bacteria; Actinobacteria; Streptomycetaceae", "Bacteria; Proteobacteria; Rhizobiaceae", "Bacteria; Actinobacteria; Kineosporiaceae", "Eukaryota; Nematoda; Mermithidae", "Bacteria; Proteobacteria; Betaproteobacteria", "Orthorvirae; Negarviricota; Iroviridae", "Bacteria; Actinobacteria; Corynebacteriales", "Bacteria; Proteobacteria; Rhodanobacteraceae", "Bacteria; Proteobacteria; Enterobacterales", "Bacteria; Bacteroidetes; Hymenobacteraceae", "Bacteria; Proteobacteria; Proteobacteria", "Bacteria; Proteobacteria; Burkholderiales", "Bacteria; Proteobacteria; Caulobacteraceae", "Bacteria; Actinobacteria; Pseudonocardiaceae", "Bacteria; Proteobacteria; Gammaproteobacteria", "Bacteria; Actinobacteria; Microbacteriaceae", "Bacteria; Proteobacteria; Oxalobacteraceae", "Bacteria; Proteobacteria; Alphaproteobacteria", "Bacteria; Actinobacteria; Nocardioidaceae", "Bacteria; Proteobacteria; Rhizobiales", "Bacteria; Proteobacteria; Comamodaceae", "Bacteria; Proteobacteria; Coxiellaceae", "Bacteria; Proteobacteria; Burkholderiaceae", "Bacteria; Proteobacteria; Yersiniaceae", "Bacteria; Proteobacteria; Anaplasmataceae", "Bacteria; Proteobacteria; Methylobacteriaceae", "Bacteria; Proteobacteria; Sphingomodaceae", "Bacteria; Actinobacteria; Actinobacteria", "Eukaryota; Apicomplexa; Piroplasmida", "Bacteria; Proteobacteria; Pseudomodaceae", "Bacteria; Actinobacteria; Mycobacteriaceae", "Eukaryota; Euglenozoa; Trypanosomatidae", "Eukaryota; Apicomplexa; Theileriidae", "Bacteria; Proteobacteria; Candidatus Midichloriaceae"))

hm_Ihol_rel
ggsave("Ihol_heatmap.pdf", path="./plots/Final_figs", width = 25, height = 30, units = "cm")
```

H. bancrofti (Lifestage x site)
```{r}
av2_Hban <- amp_subset_samples(av2_filt, Sample_species == "Haemaphysalis bancrofti", removeAbsents = TRUE)

hm_Hban_rel_txt <- amp_heatmap(av2_Hban, group_by = "Sample_species", 
                               tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 50, 
                               normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                               plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                               plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE, textmap = TRUE)
write.csv(hm_Hban_rel_txt, file = "hm_rel_txt_Hban.csv")

hm_Hban_rel_txt <- amp_heatmap(av2_Hban, facet_by = "Site", group_by = "Lifestage", 
                               tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 50, 
                               normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                               plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                               plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(strip.text.x = element_text(size = 10, face = "italic"))
  #scale_y_discrete(labels = c("Remaining taxa (358)", "Bacteria; Proteobacteria; Enterobacterales", "Bacteria; Bacteroidetes; #Chitinophagaceae", "Bacteria; Proteobacteria; Sphingomodales", "Bacteria; Planctomycetes; Isosphaeraceae", "Bacteria; Armatimodetes; #Fimbriimodaceae", "Bacteria; Proteobacteria; Gammaproteobacteria", "Bacteria; Proteobacteria; Lichenihabitantaceae", "Bacteria; #Actinobacteria; Corynebacteriales", "Bacteria; Actinobacteria; Micrococcales", "Bacteria; Actinobacteria; Nocardiaceae", "Bacteria; #Proteobacteria; Aurantimodaceae", "Bacteria; Acidobacteria; Acidobacteriaceae", "Bacteria; Actinobacteria; Micromonosporaceae", "Bacteria; #Proteobacteria; Hyphomicrobiaceae", "Eukaryota; Apicomplexa; Eucoccidiorida", "Bacteria; Actinobacteria; Frankiales", "Bacteria; #Bacteroidetes; Sphingobacteriaceae", "Bacteria; Proteobacteria; Rhizobiaceae", "Bacteria; Proteobacteria; Betaproteobacteria", "Bacteria; #Proteobacteria; Bradyrhizobiaceae", "Bacteria; Proteobacteria; Caulobacteraceae", "Bacteria; Proteobacteria; Anaplasmataceae", "Bacteria; #Proteobacteria; Rhodobacteraceae", "Bacteria; Actinobacteria; Kineosporiaceae", "Bacteria; Proteobacteria; Pseudomodaceae", "Bacteria; #Proteobacteria; Burkholderiaceae", "Bacteria; Proteobacteria; Oxalobacteraceae", "Bacteria; Actinobacteria; Geodermatophilaceae", "Bacteria; #Proteobacteria; Proteobacteria", "Bacteria; Proteobacteria; Burkholderiales", "Bacteria; Proteobacteria; Acetobacteraceae", "Bacteria; #Bacteroidetes; Hymenobacteraceae", "Bacteria; Proteobacteria; Beijerinckiaceae", "Bacteria; Bacteroidetes; Cytophagaceae", "Bacteria; #Actinobacteria; Microbacteriaceae", "Bacteria; Proteobacteria; Coxiellaceae", "Bacteria; Proteobacteria; Alphaproteobacteria", "Bacteria; #Actinobacteria; Nocardioidaceae", "Bacteria; Proteobacteria; Comamodaceae", "Bacteria; Proteobacteria; Rhizobiales", "Bacteria; #Actinobacteria; Pseudonocardiaceae", "Eukaryota; Apicomplexa; Babesiidae", "Eukaryota; Apicomplexa; Hepatozoidae", "Bacteria; Actinobacteria; #Actinobacteria", "Bacteria; Proteobacteria; Methylobacteriaceae", "Bacteria; Actinobacteria; Mycobacteriaceae", "Bacteria; Proteobacteria; #Sphingomodaceae", "Eukaryota; Euglenozoa; Trypanosomatidae", "Bacteria; Proteobacteria; Rickettsiaceae", "Bacteria; Proteobacteria; #Francisellaceae"))

hm_Hban_rel_txt
#ggsave("Hban_heatmap.pdf", path="./plots/Final_figs", width = 25, height = 30, units = "cm")

```

Ixodes tichosuri only
```{r}
av2_Itri <- amp_subset_samples(av2_filt, Sample_species == "Ixodes holocyclus", removeAbsents = TRUE)

hm_Itri_rel_txt <- amp_heatmap(av2_Itri, group_by = "Lifestage", 
                               tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                               normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                               plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                               plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE, textmap = TRUE)
write.csv(hm_Itri_rel_txt, file = "hm_rel_txt_Itri.csv")

hm_Itri_rel_txt <- amp_heatmap(av2_Itri, group_by = "Sample_species", 
                               tax_aggregate = "Family", tax_add = c("Kingdom", "Phylum"), tax_show = 30, 
                               normalise =  TRUE, plot_values = TRUE, plot_values_size = 2, round = 2, 
                               plot_colorscale = "log10", showRemainingTaxa = TRUE, 
                               plot_legendbreaks = c(0.1, 1, 10), plot_na = TRUE) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=10)) + 
  theme(strip.text.x = element_text(size = 10, face = "italic"))
  #scale_y_discrete(labels = c())

hm_Itri_rel_txt
ggsave("Itri_heatmap.pdf", path="./plots/heatmaps", width = 25, height = 30, units = "cm")
```


### Alpha diversity plots
Ixodes holocyclus lifestages vs sites
```{r}
Ihol_aplot_chao1 <- plot_diversity_stats(ps_filt_Ihol, 
                                         label.format = "p.format", 
                                         group="Species_Lifestage_Site", 
                                         index = "chao1", 
                                         group.colors = c("#E69F00", 
                                                          "#56B4E9", 
                                                          "#009E73", 
                                                          "#E69F00", 
                                                          "#56B4E9", 
                                                          "#009E73"), 
                                         group.order = c("Ixodes holocyclus; Female; Kioloa", 
                                                         "Ixodes holocyclus; Male; Kioloa", 
                                                         "Ixodes holocyclus; Nymph; Kioloa", 
                                                         "Ixodes holocyclus; Female; Sydney", 
                                                         "Ixodes holocyclus; Male; Sydney", 
                                                         "Ixodes holocyclus; Nymph; Sydney")) + 
  scale_x_discrete(labels = c("Ixodes holocyclus; Female; Kioloa" = "Female",
                              "Ixodes holocyclus; Male; Kioloa" = "Male", 
                              "Ixodes holocyclus; Nymph; Kioloa" = "Nymph", 
                              "Ixodes holocyclus; Female; Sydney" = "Female", 
                              "Ixodes holocyclus; Male; Sydney" = "Male", 
                              "Ixodes holocyclus; Nymph; Sydney" = "Nymph")) + 
 # facet_wrap(~Site, scales="free") + 
  theme_bw() +  
  ylab(label = "chao1") + 
  xlab(label = "")
  #ylim(0, 525) 

Ihol_aplot_shannon <- plot_diversity_stats(ps_filt_Ihol, 
                                           label.format = "p.format", 
                                           group="Species_Lifestage_Site", 
                                           index = "diversity_shannon", 
                                           group.colors = c("#E69F00", 
                                                            "#56B4E9", 
                                                            "#009E73", 
                                                            "#E69F00", 
                                                            "#56B4E9",
                                                            "#009E73"), 
                                           group.order = c("Ixodes holocyclus; Female; Kioloa", 
                                                           "Ixodes holocyclus; Male; Kioloa", 
                                                           "Ixodes holocyclus; Nymph; Kioloa", 
                                                           "Ixodes holocyclus; Female; Sydney", 
                                                           "Ixodes holocyclus; Male; Sydney", 
                                                           "Ixodes holocyclus; Nymph; Sydney")) + 
  scale_x_discrete(labels = c("Ixodes holocyclus; Female; Kioloa" = "Female", 
                              "Ixodes holocyclus; Male; Kioloa" = "Male", 
                              "Ixodes holocyclus; Nymph; Kioloa" = "Nymph", 
                              "Ixodes holocyclus; Female; Sydney" = "Female", 
                              "Ixodes holocyclus; Male; Sydney" = "Male", 
                              "Ixodes holocyclus; Nymph; Sydney" = "Nymph")) + 
  #facet_wrap(~Site, scales="free") + 
  theme_bw() + 
  ylab(label = "shannon") + 
  xlab(label = "")
  #ylim(0, 3)

Iho_aplot_chao1_shan <- ggarrange(Ihol_aplot_chao1, 
                                  Ihol_aplot_shannon, 
                                  ncol = 1, 
                                  nrow = 2, 
                                  labels = c("A", ""), align = "h")
Iho_aplot_chao1_shan
ggsave("I_holocyclus_alpha_div_chao1_shan.pdf", 
       path="./plots/alpha_diversity", 
       width = 30, height = 30, units = "cm")
```

Haemaphysalis bancrofti lifestages vs sites
```{r}
Hban_aplot_chao1 <- plot_diversity_stats(ps_filt_Hban_Kioloa, 
                                         label.format = "p.format", 
                                         group = "Species_Lifestage_Site", 
                                         index = "chao1", 
                                         group.colors = c("#E69F00", 
                                                          "#56B4E9", 
                                                          "#009E73"), 
                                         group.order = c("Haemaphysalis bancrofti; Female; Kioloa", 
                                                         "Haemaphysalis bancrofti; Male; Kioloa", 
                                                         "Haemaphysalis bancrofti; Nymph; Kioloa")) + 
  scale_x_discrete(labels = c("Haemaphysalis bancrofti; Female; Kioloa" = "Female", 
                              "Haemaphysalis bancrofti; Male; Kioloa" = "Male", 
                              "Haemaphysalis bancrofti; Nymph; Kioloa" = "Nymph")) + 
  facet_wrap(~Site, scales = "free") + 
  theme_bw() + 
  ylab(label = "chao1") + 
  xlab(label = "") + 
  ylim(0,800)

Hban_aplot_shan <- plot_diversity_stats(ps_filt_Hban_Kioloa, 
                                        label.format = "p.format", 
                                         group = "Species_Lifestage_Site", 
                                         index = "diversity_shannon", 
                                         group.colors = c("#E69F00", 
                                                          "#56B4E9", 
                                                          "#009E73"), 
                                         group.order = c("Haemaphysalis bancrofti; Female; Kioloa", 
                                                         "Haemaphysalis bancrofti; Male; Kioloa", 
                                                         "Haemaphysalis bancrofti; Nymph; Kioloa")) + 
  scale_x_discrete(labels = c("Haemaphysalis bancrofti; Female; Kioloa" = "Female", 
                              "Haemaphysalis bancrofti; Male; Kioloa" = "Male", 
                              "Haemaphysalis bancrofti; Nymph; Kioloa" = "Nymph")) + 
  facet_wrap(~Site, scales = "free") + 
  theme_bw() + 
  ylab(label = "shannon") + 
  xlab(label = "") + 
  ylim(0,3.7)

Hban_aplot_chao1_shan <- ggarrange(Hban_aplot_chao1, 
                                   Hban_aplot_shan, 
                                   ncol = 1, 
                                   nrow = 2, 
                                   labels = c("B", ""), align = "h")
Hban_aplot_chao1_shan
ggsave("H_bancrofti_alpha_div_chao1_shan.pdf", 
       path="./plots/alpha_diversity", 
       width = 15, height = 30, units = "cm")
```

Wildlife bloods, by species x site
```{r}
Blood_aplot_chao1 <- plot_diversity_stats(ps_filt_bloods, 
                                          label.format = "p.format", 
                                          group = "Sample_species", 
                                          index = "chao1", 
                                          group.colors = c("#D55E00", 
                                                           "#CC79A7", 
                                                           "#0072B2", 
                                                           "#E69F00", 
                                                           "#009E73")) + 
  theme_bw() + 
  ylab(label = "chao1") + 
  xlab(label = "") + 
  #ylim(0,310) + 
  theme(axis.text.x = element_text(face = "italic"))


Blood_aplot_shan <- plot_diversity_stats(ps_filt_bloods, 
                                         label.format = "p.format", 
                                         group = "Sample_species", 
                                         index = "diversity_shannon", 
                                         group.colors = c("#D55E00", 
                                                          "#CC79A7", 
                                                          "#0072B2", 
                                                          "#E69F00", 
                                                          "#009E73")) + 
  theme_bw() + 
  ylab(label = "shannon") + 
  xlab(label = "") + 
  #ylim(0,0.21) + 
  theme(axis.text.x = element_text(face = "italic"))
Blood_aplot_chao1_shan <- ggarrange(Blood_aplot_chao1, 
                                    Blood_aplot_shan, 
                                    ncol = 1, 
                                    nrow = 2, 
                                    labels = c("C", ""))
Blood_aplot_chao1_shan
ggsave("Bloods_alpha_div_chao1_shan.pdf", 
       path="./plots/Final_figs", 
       width = 15, height = 15, units = "cm")
```

Put ticks alpha diversity together
```{r}
library("cowplot")
y <- ggdraw() + 
  draw_plot(Iho_aplot_chao1_shan, x = 0, y = 0.5, width = 1, height = 0.5) + 
  draw_plot(Hban_aplot_chao1_shan, x = 0, y = 0, width = 0.5, height = 0.5) + 
  draw_plot(Blood_aplot_chao1_shan, x = 0.5, y = 0, width = 0.5, height = 0.5)
y

ggsave("Ticks_alpha_div.pdf", path = "./plots/Final_figs", height = 18, width = 25, units = "cm")
```


### Ordination plots with ampvis2

All samples PCoA
```{r}
av2_all_ord <- amp_ordinate(av2_filt, filter_species = 0,type = "PCoA", 
                            distmeasure = "jsd", 
                            x_axis = 1, y_axis = 2,print_caption = FALSE, sample_color_by = "Site", 
                            sample_shape_by = "Sample_species", opacity = 0.7, sample_colorframe = "Sample_species")+
  theme_bw() + 
  scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values=c(16, 17, 15, 0, 1, 6, 5))
av2_all_ord    
ggsave("all_samples_PCA.pdf", path = "./plots/Final_figs", width = 25, height = 20, units = "cm")
```
All ord without endosymbiotns
```{r}
endo_tax <- c("f__Candidatus Midichloriaceae", "g__Candidatus Midichloria", "f__Francisellaceae", "g__Francisella", "f__Rickettsiaceae", 
              "g__Rickettsia")
av2_noendo <- amp_subset_taxa(av2_unass_filt, tax_vector = endo_tax, normalise = FALSE, remove = TRUE)

av2_no_endo_ord <- amp_ordinate(av2_noendo, filter_species = 0,type = "PCoA", 
                            distmeasure = "jsd", 
                            x_axis = 1, y_axis = 2,print_caption = FALSE, sample_color_by = "Site", 
                            sample_shape_by = "Sample_species", opacity = 0.7, sample_colorframe = "Sample_species")+
  theme_bw() + 
  scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values=c(16, 17, 15, 0, 1, 6, 5))
av2_no_endo_ord    
ggsave("all_samples_no_endo_PCA.pdf", path = "./plots/Final_figs", width = 25, height = 20, units = "cm")
```




Ixodes holocyclus only PCoA
```{r}
av2_Ihol <- amp_subset_samples(av2_filt, Sample_species == "Ixodes holocyclus", removeAbsents = TRUE)
av2_Ihol_ord <- amp_ordinate(av2_Ihol, filter_species = 0, type = "PCoA", 
                             distmeasure = "jsd", 
                             x_axis = 1, y_axis = 2, sample_color_by = "Site", 
                             sample_shape_by = "Lifestage", opacity = 0.7, sample_colorframe = "Lifestage") + 
  theme_bw() + 
  scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values=c(16, 17, 15))
av2_Ihol_ord    
ggsave("Ixodes_holocyclus_PCoA.pdf", path = "./plots/PCoA", width = 25, height = 25, units = "cm")
```

H. bancrofti only PCoA
```{r}
av2_Hban<- amp_subset_samples(av2_filt, Sample_species == "Haemaphysalis bancrofti", removeAbsents = TRUE)
av2_Hban_ord <- amp_ordinate(av2_Hban, filter_species = 0, type = "PCoA", 
                             distmeasure = "jsd", 
                             x_axis = 1, y_axis = 2, sample_color_by = "Site", 
                             sample_shape_by = "Lifestage", opacity = 0.7, sample_colorframe = "Lifestage") + 
  theme_bw() + 
  scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values=c(16, 17, 15))

av2_Hban_ord    
ggsave("Haemaphysalis_bancrofti_PCoA.pdf", path = "./plots/PCoA", width = 25, height = 25, units = "cm")
```

PCoA animal blood samples
```{r}
av2_bloods_ord <- amp_ordinate(av2_bloods_filt, filter_species = 0, type = "PCoA", 
                               distmeasure = "jsd", 
                               x_axis = 1, y_axis = 2, print_caption = FALSE, 
                               sample_color_by = "Site", sample_shape_by = "Sample_species", 
                               opacity = 0.7, sample_colorframe = "Species_Site") + 
  theme_bw() + 
  scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values=c(5, 1, 2, 0))
av2_bloods_ord    
ggsave("Wildlife_bloods_PCoA.pdf", path = "./plots/PCA", width = 25, height = 25, units = "cm")
```
Join PCoA for Ihol and Hban
```{r}
library("cowplot")
y <- ggarrange(av2_Ihol_ord, av2_Hban_ord, ncol = 1, nrow = 2, labels = c("A", "B"))
y

ggsave("Figure S3 PcoA plots.pdf", path = "./plots/Final_figs", height = 21, width = 28, units = "cm")
```




```{r}
source("taxa_summary.R", local = T)
```
### Prevalence vs. total Count scatter plots

Ixodes holocyclus
```{r}
fm_Ihol = fast_melt(ps_filt_Ihol)
prevdt = fm_Ihol[, list(Prevalence = sum(count > 0), 
                        TotalCounts = sum(count)), 
                   by = TaxaID]
addPhylum = unique(copy(fm_Ihol[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
# Top 10 phyla
showPhyla = prevdt[, sum(TotalCounts), by = Phylum][order(-V1)][1:10]$Phylum
setkey(prevdt, Phylum)
# Make plots for phyla
top_phyla1 = ggplot(prevdt[showPhyla], 
                    mapping = aes(Prevalence, 
                                  TotalCounts, 
                                  color = Phylum)) + 
  geom_point(size = 1.5, alpha = 0.75) + 
  scale_y_log10() + 
  theme_bw() + 
  scale_color_brewer(palette = "Spectral") + 
  ylab("Abundance") + 
  labs(title = "Ixodes holocyclys, all samples")
top_phyla1
ggsave("top_phyla_Ihol.pdf", plot = top_phyla1, path = "./plots", width = 15, height = 10, units = "cm")
```

Ixodes holocyclus Kioloa vs Sydney
```{r}
# Kioloa
ps_filt_Ihol_Kioloa <- subset_samples(ps_filt_Ihol, Site == "Kioloa")
fm_Ihol_Kioloa = fast_melt(ps_filt_Ihol_Kioloa)
prevdt = fm_Ihol_Kioloa[, list(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]
addPhylum = unique(copy(fm_Ihol_Kioloa[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
# Top 10 phyla
showPhyla = prevdt[, sum(TotalCounts), by = Phylum][order(-V1)][1:10]$Phylum
setkey(prevdt, Phylum)
# Make plots for phyla
top_phyla_Ihol_Kioloa = ggplot(prevdt[showPhyla], 
                               mapping = aes(Prevalence, 
                                             TotalCounts, 
                                             color = Phylum)) + 
  geom_point(size = 1.5, alpha = 0.75) + 
  scale_y_log10() + 
  theme_bw() + 
  scale_color_brewer(palette = "Spectral") + 
  ylab("Abundance") + 
  labs(title = "Ixodes holocyclys, Kioloa")

# Sydney
ps_filt_Ihol_Syd <- subset_samples(ps_filt_Ihol, Site == "Sydney")
fm_Ihol_Syd = fast_melt(ps_filt_Ihol_Syd)
prevdt = fm_Ihol_Syd[, list(Prevalence = sum(count > 0), 
                            TotalCounts = sum(count)), 
                     by = TaxaID]
addPhylum = unique(copy(fm_Ihol_Syd[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
# Top 10 phyla
showPhyla = prevdt[, sum(TotalCounts), by = Phylum][order(-V1)][1:10]$Phylum
setkey(prevdt, Phylum)
# Make plots for phyla
prevdt[showPhyla]
top_phyla_Ihol_Sydney = ggplot(prevdt[showPhyla], 
                               mapping = aes(Prevalence, 
                                             TotalCounts, 
                                             color = Phylum)) + 
  geom_point(size = 1.5, alpha = 0.75) + 
  scale_y_log10() + 
  theme_bw() + 
  scale_color_brewer(palette = "Spectral") + 
  ylab("Abundance") + 
  labs(title = "Ixodes holocyclys, Sydney")

top_phyla_Ihol_Kioloa
top_phyla_Ihol_Sydney

ggsave("top_phyla_Ihol_Kioloa.pdf", 
       plot = top_phyla_Ihol_Kioloa, 
       path = "./plots", 
       width = 15, height = 10, units = "cm")
ggsave("top_phyla_Ihol_Sydney.pdf", 
       plot = top_phyla_Ihol_Sydney, 
       path = "./plots", 
       width = 15, height = 10, units = "cm")

z <- ggarrange(top_phyla_Ihol_Kioloa, top_phyla_Ihol_Sydney)
z
```

Ixodes holocyclus Families
```{r}
# Get top 10 most prevalent families
mdt_Ihol = fast_melt(ps_filt_Ihol)
prevdt_Ihol = mdt_Ihol[, list(Prevalence = sum(count > 0),
                              TotalCounts = sum(count)),
                              by = TaxaID]
addFamily_Ihol = unique(copy(mdt_Ihol[, list(TaxaID), Family]))
setkey(prevdt_Ihol, TaxaID)
setkey(addFamily_Ihol, TaxaID)
prevdt_Ihol <- addFamily_Ihol[prevdt_Ihol]
showFamily_Ihol = prevdt_Ihol[, sum(TotalCounts), by = Family] [order(-V1)][1:10]$Family
setkey(prevdt_Ihol, Family)

# Make abundance vs prevalence plot
top_families_Ihol <- ggplot(prevdt_Ihol[showFamily_Ihol], 
                            mapping = aes(Prevalence, TotalCounts, color = Family)) + 
  geom_point(size = 1.5, alpha = 0.75) + 
  scale_y_log10() + 
  theme_bw() + 
  ylab("Abundance") + 
  labs(title = "Ixodes holocyclys, Families") #+ 
  #theme(legend.position = "none")
  
top_families_Ihol
x <- prevdt_Ihol[showFamily_Ihol]

```

H. bancrofti
```{r}
fm_Hban = fast_melt(ps_filt_Hban)
prevdt = fm_Hban[, list(Prevalence = sum(count > 0), 
                          TotalCounts = sum(count)), 
                   by = TaxaID]
addPhylum = unique(copy(fm_Hban[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
# Top 10 phyla
showPhyla = prevdt[, sum(TotalCounts), by = Phylum][order(-V1)][1:10]$Phylum
setkey(prevdt, Phylum)
# Make plots for phyla
top_phyla1 = ggplot(prevdt[showPhyla], 
                    mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
             geom_point(size = 1.5, alpha = 0.75) + 
             scale_y_log10() + 
             theme_bw() + 
             scale_color_brewer(palette = "Spectral") + 
             ylab("Abundance") + 
             labs(title = "Overall")
top_phyla1
ggsave("top_phyla_Ihol.pdf", plot = top_phyla1, path = "./plots", width = 15, height = 10, units = "cm")
```


### Plotting tick-borne micobes

```{r}
ps_Fam <- aggregate_taxa(ps_filt, "Family")
ps_Fam_df <- psmelt(ps_Fam)
write_csv(ps_Fam_df, "ps_Family_aggregated.csv")

cPal = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
         "#D55E00", "#CC79A7")
```

Family_Anaplasmataceae
```{r}
plot_taxa_Ana <- plot_listed_taxa(ps_Fam, 
                              "Anaplasmataceae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Ana
ggsave("Anaplasmataceae_abundance_plots.pdf", 
       path="./plots/abundance_plots", 
       height=20, width=20, units="cm")
```

Family_Borreliaceae
```{r}
plot_taxa_Bor <- plot_listed_taxa(ps_Fam, 
                              "Borreliaceae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Bor
ggsave("Borreliaceae_abundance_plots.pdf", 
       path="./plots/abundance_plots", 
       height=20, width=20, units="cm")
```

Family_Rickettsiaceae
```{r}
plot_taxa_Rick <- plot_listed_taxa(ps_Fam, 
                              "Rickettsiaceae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Rick
ggsave("Rickettsiaceae_abundance_plots.pdf", 
       path="./plots/abundance_plots", 
       height=20, width=20, units="cm")
```
Family_Midichloriaceae
```{r}
plot_taxa_Midichloria <- plot_listed_taxa(ps_Fam, 
                              "Candidatus Midichloriaceae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Midichloria
ggsave("Midichloriaceae_abundance_plots.pdf", path="./plots/abundance_plots", height=20, width=20, units="cm")
```
Family_Francicellaceae
```{r}
plot_taxa_Francisellaceae <- plot_listed_taxa(ps_Fam, 
                              "Francisellaceae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Francisellaceae
ggsave("Francisellaceae_abundance_plots.pdf", path="./plots/abundance_plots", height=20, width=20, units="cm")
```
Family Babesiidae
```{r}
plot_taxa_Babesiidae <- plot_listed_taxa(ps_Fam, 
                              "Babesiidae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Babesiidae
ggsave("Babesiidae_abundance_plots.pdf", path="./plots/abundance_plots", height=20, width=20, units="cm")
```
Family Theileridae
```{r}
plot_taxa_Theileriidae <- plot_listed_taxa(ps_Fam, 
                              "Theileriidae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Theileriidae
ggsave("Theileriidae_ticks_abundance_plots.pdf", path="./plots/abundance_plots", height=20, width=20, units="cm")
```
Family Trypanosomatidae
```{r}
plot_taxa_Trypanosomatidae <- plot_listed_taxa(ps_Fam, 
                              "Trypanosomatidae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Trypanosomatidae
ggsave("Trypanosomatidae_ticks_abundance_plots.pdf", path="./plots/abundance_plots", height=20, width=20, units="cm")
```
Family f__Bartonellaceae
```{r}
plot_taxa_Bartonellaceae <- plot_listed_taxa(ps_Fam, 
                              "Bartonellaceae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Bartonellaceae
ggsave("Bartonellaceae_ticks_abundance_plots.pdf", 
       path="./plots/abundance_plots", 
       height=20, width=20, units="cm")
```
Family Hepatozoidae
```{r}
plot_taxa_Hepatozoidae <- plot_listed_taxa(ps_Fam, 
                              "Hepatozoidae", 
                              group = "Sample_species", 
                              group.colors = cPal, 
                              add.violin = FALSE,
                              dot.opacity = 0.7,
                              box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 25, 
                                   hjust = 1, 
                                   vjust = 0.95, 
                                   face = "italic", 
                                   size = 8), 
        axis.title.x = element_blank())
plot_taxa_Hepatozoidae
ggsave("Hepatozoidae_ticks_abundance_plots.pdf", 
       path="./plots/abundance_plots", 
       height=20, width=20, units="cm")
```

Putting them all togeher
```{r}
abund_plot <- ggarrange(plot_taxa_Midichloria, 
                        plot_taxa_Francisellaceae, 
                        plot_taxa_Rick, 
                        plot_taxa_Ana,
                        plot_taxa_Bartonellaceae,
                        plot_taxa_Hepatozoidae, 
                        plot_taxa_Babesiidae, 
                        plot_taxa_Theileriidae, 
                        plot_taxa_Trypanosomatidae)
abund_plot
ggsave("TBP_abundance_plot.pdf", 
       path = "./plots/abundance_plots", 
       height = 20, width = 20, units = "cm")
```
Viruses
```{r}
# Rhabdoviridae
plot_taxa_Rhabdoviridae <- plot_listed_taxa(ps_Fam, "Rhabdoviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Rhabdoviridae
# Partitiviridae
plot_taxa_Partitiviridae <- plot_listed_taxa(ps_Fam, "Partitiviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Partitiviridae
# Reoviridae
plot_taxa_Reoviridae <- plot_listed_taxa(ps_Fam, "Reoviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Reoviridae
# Virgaviridae
plot_taxa_Virgaviridae <- plot_listed_taxa(ps_Fam, "Virgaviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Virgaviridae
# Phenuiviridae
plot_taxa_Phenuiviridae <- plot_listed_taxa(ps_Fam, "Phenuiviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Phenuiviridae
# Chuviridae
plot_taxa_Chuviridae <- plot_listed_taxa(ps_Fam, "Chuviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Chuviridae
# Iflaviridae
plot_taxa_Iflaviridae <- plot_listed_taxa(ps_Fam, "Iflaviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Iflaviridae
# Circoviridae 
plot_taxa_Circoviridae <- plot_listed_taxa(ps_Fam, "Circoviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Circoviridae
# Anelloviridae  
plot_taxa_Anelloviridae <- plot_listed_taxa(ps_Fam, "Anelloviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Anelloviridae  
# Flaviviridae   
plot_taxa_Flaviviridae <- plot_listed_taxa(ps_Fam, "Flaviviridae", 
  group = "Sample_species", group.colors = cPal, add.violin = FALSE,
  dot.opacity = 0.7, box.opacity = 0.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.95, face = "italic", size = 7), 
                      axis.title.x = element_blank())
plot_taxa_Flaviviridae 
```

Put them all together
```{r}
viral_abund_plots <- ggarrange(plot_taxa_Rhabdoviridae, 
                               plot_taxa_Partitiviridae, 
                               plot_taxa_Reoviridae, 
                               plot_taxa_Virgaviridae, 
                               plot_taxa_Phenuiviridae, 
                               plot_taxa_Chuviridae, 
                               plot_taxa_Iflaviridae, 
                               plot_taxa_Circoviridae, 
                               plot_taxa_Anelloviridae, 
                               plot_taxa_Flaviviridae)
viral_abund_plots
ggsave("Viruses_abundance_plot.pdf", 
       path = "./plots/abundance_plots", 
       height = 20, width = 20, units = "cm")
```

z
```{r}
ps_unass <- subset_taxa(ps, Domain != "k__Unassigned")

y <- plot_taxa_boxplot(ps_Fam, 
                       taxonomic.level = "Family", 
                       top.otu = 10, 
                       keep.other = FALSE, 
                       group = "Sample_species", 
                       group.colors = cPal, 
                       add.violin = TRUE,
                       violin.opacity = 0.5, 
                       dot.opacity = 0.8, 
                       dot.size = 1.5)
                       
y
ggsave("taxa_composition_plot_Family.pdf", path = "plots", 
       height = 30, width = 30, units = "cm")
```





























