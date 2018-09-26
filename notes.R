raw_files <- list.files("raw/", pattern = "*.txt") # get all raw data from txt fil
raw_files
counts <- lapply(paste0("raw/",raw_files), 
                   function(x) {
                     readr::read_tsv(x, # read tsv
                                     col_names = c("gene_id", "count"), 
                                     col_types = readr::cols(
                                       gene_id = readr::col_character(),
                                       count = readr::col_integer()
                                       )
                                     )
                     }  # it seems like to determine the colummns characters
                   ) %>% 
  purrr::reduce(dplyr::left_join, by = 'gene_id') # AND combine all from left join by gene_id, what is purrr::reduce?
#ok, after we tried the one without "purrr::reduce", the data is still seperated lists
colnames(counts) <- c("gene_id", gsub("_accepted_hits.bam.read_sorted.bam.sam.count.txt", "", raw_files) )
# this one is to add column names by using the row files, what is gsub? let me have a look
#OK, it is like 'grep' but for this case it takes the stuff from the determined files
counts
extra_data <- readr::read_csv("raw/Mo_CMcm_reads_counts.csv",skip = 1,
                              col_names = c("gene_id", "CM_1", "CM_2", "CM_3"), 
                              col_types = readr::cols(
                                       gene_id = readr::col_character(),
                                       CM_1 = readr::col_integer(),
                                       CM_2 = readr::col_integer(),
                                       CM_3 = readr::col_integer()
                                      )
                              )
#let me compare with the previous one. what is skip?
extra_data
