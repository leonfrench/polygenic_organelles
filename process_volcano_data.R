library(here)
library(readr)
library(magrittr)
library(tibble)
library(tidyr)
library(dplyr)
#volcano plot data grabbed from https://organelles.czbiohub.org/

full_matrix <- tibble()
for (filename in list.files(here("data/organelles.czbiohub/volcano_data"))) {
  print(filename)
  single_volcano <- read_csv(here("data/organelles.czbiohub/volcano_data", filename))
  col_name <- gsub("volcano_data_", "", filename)
  col_name <- gsub(".csv", "", col_name)
  print(col_name)
  single_volcano %<>% rename(gene_symbol = `Gene name canonical`)
  single_volcano %<>% mutate(sign_log_p = sign(enrichment) * `-log10(p-value)`)
  single_volcano %>% arrange(sign_log_p)
  single_volcano %<>% select(gene_symbol, sign_log_p)
  single_volcano %<>% mutate(source = col_name)
  #single_volcano %>% rename(!!col_name := sign_log_p)
  full_matrix <- bind_rows(full_matrix, single_volcano)
}
full_matrix
full_matrix %<>% group_by(source) %>% mutate(rank = rank(sign_log_p))
full_matrix %>% arrange(sign_log_p)
full_matrix %>% filter(grepl("[(]", gene_symbol))
full_matrix %<>% mutate(gene_symbol = gsub(" [(](.*)[)]","", gene_symbol))
full_matrix %>% filter(grepl("[(]", gene_symbol))

#deal with dupes introduced after removing parantheses
wide_matrix <- full_matrix %>% group_by(gene_symbol, source) %>% summarize(rank = mean(rank)) %>%
  pivot_wider(names_from = source, id_cols = c(gene_symbol), values_from = rank)

dim(wide_matrix)
wide_matrix %>% write_csv(here("data/organelles.czbiohub/processed_ranks_markers.csv"))

cell_meta_data <- tibble(cell_type = colnames(wide_matrix %>% select(-gene_symbol)))
cell_meta_data %<>% mutate(Organelle = gsub(".* [(](.*)[)]","\\1", cell_type))
cell_meta_data %<>% mutate(Marker = gsub(" [(](.*)[)]","", cell_type)) 
cell_meta_data %>% write_csv(here("data/organelles.czbiohub/processed_ranks_markers_metadata.csv"))

full_matrix %<>% inner_join(cell_meta_data %>% rename(source=cell_type))
full_matrix %>% select(Organelle, Marker) %>% distinct() %>% 
  group_by(Organelle) %>% summarize(Markers = paste(sort(Marker), collapse=",")) %>% rename(cell_type=Organelle) %>%
  write_csv(here("data/organelles.czbiohub/processed_ranks_organelle_metadata.csv"))
full_matrix %<>% group_by(Organelle, gene_symbol) %>% summarize(rank = mean(rank))
full_matrix %>% tail()
#re-rank
full_matrix %<>% group_by(Organelle) %>% mutate(rank = rank(rank))
full_matrix %>% arrange(rank)

wide_matrix_organelle <- full_matrix %>% pivot_wider(names_from = Organelle, id_cols = c(gene_symbol), values_from = rank)

wide_matrix_organelle %>% write_csv(here("data/organelles.czbiohub/processed_ranks_organelle.csv"))

