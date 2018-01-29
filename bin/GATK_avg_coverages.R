#!/usr/bin/env Rscript
# This script will parse the '.sample_interval_summary' files output from GATK DepthOfCoverage
# and give a formatted table of coverage stats across samples

# USAGE: $ find example_input/ -type f -name "*.sample_interval_summary" | xargs ./GATK_avg_coverages.R


# ~~~~~ FUNCTIONS ~~~~~ # 
source("snsxt_tools.R")

read.sample_interval_summary.file <- function(file){
    # read the GATK DepthOfCoverage sample_interval_summary file
    # coverage_df <- read.delim(file = coverage_file, header = TRUE, sep = ',')
    df <- read.delim(file = file, header = TRUE, sep = ',', check.names = FALSE)
    return(df)
}

sampleID_from_sample_interval_summary_file <- function(file){
    # get the sample ID for the file from the 4th header field
    header_line <- readLines(file)
    header_field <- unlist(strsplit(x = header_line, split = ','))[4]
    sampleID <- gsub(pattern = '_total_cvg$', replacement = '', x = header_field)
    return(sampleID)
}

build_all_coverages_df <- function(coverage_files){
    # aggregate all the coverage files into a single df
    # empty df to hold all the coverages
    all_coverages_df <- data.frame()
    
    for(coverage_file in coverage_files){
        # load file
        tsprintf("Reading from coverage file: %s", coverage_file)
        coverage_df <- read.sample_interval_summary.file(coverage_file)
        # get sampleID from field in file header
        sampleID <- sampleID_from_sample_interval_summary_file(coverage_file)
        # save the Targets as rownames
        rownames(coverage_df) <- coverage_df[["Target"]]
        # save just the 'mean_cvg' column
        coverage_df <- coverage_df[5]
        # reset column name to sampleID
        colnames(coverage_df)[1] <- sampleID
        
        # load the data into the overall df
        if(nrow(all_coverages_df) == 0){
            all_coverages_df <- coverage_df
        } else {
            all_coverages_df <- cbind(all_coverages_df, coverage_df)
        }
    }
    return(all_coverages_df)
}

# ~~~~~ RUN ~~~~~ #
# get script args
args <- commandArgs(trailingOnly = TRUE)  

print(args)
coverage_files <- args


# ~~~~~~~ IMPORT AVERAGE COVERAGE PER REGION PER SAMPLE ~~~~~~~ #
all_coverages_df <- build_all_coverages_df(coverage_files)
avg_file <- 'average_coverage_per_sample.tsv'
tsprintf("Writing sample averages to file: %s", avg_file)
write.table(x = all_coverages_df, sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA, 
            file = avg_file)

# ~~~~~~~ CALCULATE AVERAGE OF AVG'S PER REGION ~~~~~~~ #
region_coverages_df <- as.data.frame(rowMeans(all_coverages_df))
colnames(region_coverages_df) <- "average_coverage"
region_avg_file <- 'average_coverage_per_region.tsv'
tsprintf("Writing region averages to file: %s", region_avg_file)
write.table(x = region_coverages_df, sep = '\t', quote = FALSE, row.names = TRUE, col.names = FALSE, file = region_avg_file)



# ~~~~~~~ CREATE BED FOR REGIONS WITH LOW COVERAGE ~~~~~~~ #
low_cutoff <- 50
tsprintf("Finding regions with coverage below %s", low_cutoff)
low_regions <- region_coverages_df[region_coverages_df["average_coverage"] < low_cutoff, , drop = FALSE]
tsprintf("Number of regions found: %s", nrow(low_regions))


low_BED_file <- sprintf('regions_coverage_below_%s.bed', low_cutoff)
try_to_save_BED(df = low_regions, output_file = low_BED_file)

tsprintf("Finding regions with coverage below 0")
zero_regions <- region_coverages_df[region_coverages_df["average_coverage"] == 0, , drop = FALSE]
tsprintf("Number of regions found: %s", nrow(zero_regions))

zero_BED_file <- "regions_with_coverage_0.bed"
try_to_save_BED(df = zero_regions, output_file = zero_BED_file)


save.image()
