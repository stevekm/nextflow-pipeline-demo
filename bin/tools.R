#!/usr/bin/env Rscript

# functions to use in snsxt task R scripts

# ~~~~~ PACKAGES ~~~~~ # 
library("optparse")
library("tools")


timestamp <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
logfile <- file.path(".", sprintf("report_log.%s.txt", timestamp))

# default value, overwrite in scripts after 'source' to enable logfile output
hard_log <- FALSE

# ~~~~~ FUNCTIONS ~~~~~ # 
tsprintf <- function(fmt, ...){
    # print a formatted message with timestamp
    # base message
    m <- sprintf(fmt, ...)
    # message with timestamp
    tm <- sprintf('[%s] %s', format(Sys.time(), "%H:%M:%S"), m)
    # emit message
    message(tm)
    # add to log file
    if(isTRUE(hard_log)) cat(sprintf("%s\n", tm), file = logfile, append = TRUE)
}

msprintf <- function(fmt, ...) {
    message(sprintf(fmt, ...))
}

mycat <- function(text){
    # print formatted text in Rmd
    cat(gsub(pattern = "\n", replacement = "  \n", x = text))
}


make_output_filename <- function(prefix, suffix, sep = '_'){
    # create an output filename, testing for conditional length of prefix
    if(nchar(prefix) > 0){
        return(paste(prefix, suffix, sep = sep))
    } else {
        return(suffix)
    }
}

chrom_regions2df <- function(regions){
    # split the regions into chrom coordinates for BED files
    # regions <- c("chr1:236998847-236998987", "chr1:237001714-237001899")
    regions_df <- as.data.frame(do.call(rbind, strsplit(regions, ':')))
    regions_df <- cbind(regions_df[1],
                        as.data.frame(do.call(rbind, strsplit(as.character(regions_df$V2), '-'))))
    colnames(regions_df) <- c("chrom", "start", "stop")
    return(regions_df)
}


write_BED <- function(df, output_file){
    write.table(x = df, file = output_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}


try_to_save_BED <- function(df, output_file) {
    if(nrow(df) > 0){
        df_bed <- chrom_regions2df(rownames(df))
        tsprintf("Writing regions to file: %s", output_file)
        write_BED(df = df_bed, output_file = output_file)
        
    } else {
        tsprintf("No regions present, making empty file: %s", output_file)
        file.create(output_file)
    }
}

# ~~~~~ RUN ~~~~~ # 
tsprintf("loaded tools.R...")