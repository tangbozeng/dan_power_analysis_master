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
#let me compare with the previous one. what is skip? need to find out later.
extra_data
counts <- dplyr::left_join(counts, extra_data, by = "gene_id") # again, combine it by adding from letf join.
summary(counts) # Some hilariously high max count values, haha, ofcs, that is the total reads number.

# It reshape the counts matrix into a combined data.frame.
counts_long <- counts %>% 
  reshape2::melt(id.vars = "gene_id", value.name = "count", variable.name=c("sample") )
# need to find out what is resharpe2

# this is the plotting method. 
counts_long %>%
  ggplot2::ggplot() +
  ggplot2::aes(count) +
  ggplot2::geom_density() +
  ggplot2::facet_wrap( . ~ sample)
# make plot to have a look of the distribution.what if we try something else?
counts_long %>%
  ggplot2::ggplot() +
  ggplot2::aes(count) +
  ggplot2::geom_density() +
  ggplot2::facet_wrap( . ~ sample) +
  ggplot2::scale_x_log10() # add log10.
m <- counts[,2:ncol(counts)]
mxs <- apply(m, MARGIN=c(2), max)
which(m == mxs, arr.ind = TRUE) # he asked whether high count genes same or not?
#
##      row col
#convert data.frame to matrix
# check missing data
m <- as.matrix(counts[,2:ncol(counts)])
# copy it
n <- m
# set 0 counts to NA
n[n < 1 ] <- NA
# How many genes with full data  in all samples
sum(rowSums(m) == 0)
ggplot2::ggplot(counts_long) + 
  ggplot2::aes(x = sample, y = gene_id) +
  ggplot2::geom_tile(ggplot2::aes(fill = log10(count)))
readr::write_csv(counts, "cleaned/counts_wide.csv")
readr::write_csv(counts_long, "cleaned/counts_long.csv")
counts <- readr::read_csv("cleaned/counts_wide.csv")
count_matrix <- as.matrix(counts[,2:ncol(counts)])
colnames(count_matrix) <- colnames(counts)[2:ncol(counts)]
rownames(count_matrix) <- counts$gene_id # important
params <- estimateParam(
  countData = count_matrix,
  Distribution = "NB",
  RNAseq = "bulk",
  normalisation = "TMM" # edgeR method, can be others
)
plotParam(params, annot = T)
sum(complete.cases(n))
colnames(counts)
#define some metadata for the sample
smp_id <- c(1,1,2,2,3,3,4,4,5,5,5) # the columns each sample type is in
smp_names <- c("guy_11", "pmk1_t00", "T24_NEG", "T24_POS", "CM") #the basename of each column
#normalise using edgeR
library(edgeR)

norm_factors <- calcNormFactors(count_matrix, method = "TMM")
dge_list <- DGEList(counts = count_matrix, norm.factors = norm_factors)
norm_counts <- cpm(dge_list)

## calculate the mean counts for each sample based on the groupings above
smp_means <- lapply(1:5, function(x, norm_counts, smp_id){
  m <- rowMeans(norm_counts[,smp_id == x])
  tibble::tibble(gene_id  = names(m), mean_count = m )
}, norm_counts, smp_id)

names(smp_means) <- smp_names

#change the list into a tibble
smp_means <-  purrr::reduce(smp_means, dplyr::left_join, by = 'gene_id')
colnames(smp_means) <- c("gene_id", smp_names)
smp_means

# get all non redundant combinations of samples
pairs <- combn(smp_names, 2 )

## define a function that given the data can calculate the  log2 fold change
## for a pair of samples
get_fc <- function(p, df = NULL){
  d <- df[,p]
  d <- log2( d[,1] / d[,2] )[,1]
  tibble::tibble(gene_id = df$gene_id, log2_fc = d )
}

## apply the function to all pairs
log2_fcs <- apply(pairs, MARGIN = 2, get_fc, smp_means)
pair_names <- apply(pairs, MARGIN = 2, function(x) paste0(x[1], "-", x[2]))
names(log2_fcs) <- pair_names

## again, change the list back to a tibble
log2_fcs <-  purrr::reduce(log2_fcs, dplyr::left_join, by = 'gene_id')
colnames(log2_fcs) <- c("gene_id", pair_names)

## convert from wide to long format
log2_fcs_long <-  tibble::as.tibble(
  reshape2::melt( log2_fcs, 
                  id.vars = "gene_id", 
                  value.name = "log2_fc", 
                  variable.name=c("comparison") 
  )
)
                    
                    library(magrittr)


## define a function that given a threshold and data of each genes
## log2 fc tells the proportion of genes with absolute log2fc over that threshold

find_over <- function(threshold, df){
  df %>% 
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      n = dplyr::n_distinct(gene_id), 
      abs_log2_fc_threshold = threshold,  
      propn_over = (sum(abs(log2_fc) >= threshold, na.rm = TRUE) / n)
    )
}

## apply the function to thresholds at 0.5 intervals
threshold_info <- dplyr::bind_rows( lapply(seq(0.5, 3, 0.5), find_over, log2_fcs_long) ) %>%
  dplyr::select( -n )

##plot

threshold_info %>%
  ggplot2::ggplot() + 
  ggplot2::aes(abs_log2_fc_threshold, propn_over) +
  ggplot2::geom_line(ggplot2::aes(colour = comparison))
                    
##plot only interesting comparisons
threshold_info %>%
  dplyr::filter(comparison %in% c("guy_11-pmk1_t00", "T24_NEG-T24_POS")) %>%
  ggplot2::ggplot() + 
  ggplot2::aes(abs_log2_fc_threshold, propn_over) +
  ggplot2::geom_line(ggplot2::aes(colour = comparison))
                    
                    ggplot2::ggplot(log2_fcs_long) + 
  ggplot2::aes(log2_fc) +
  ggplot2::geom_density()
                    
