#!/usr/bin/env Rscript

# install.packages("ncdf4")
# install.packages("dplyr")
# install.packages("reshape2")
# install.packages("RPostgreSQL")
# install.packages("optparse")
# install.packages("R.utils")
suppressPackageStartupMessages(library(ncdf4))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RPostgreSQL))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(raster))


## process command line parameters
option_list = list(
  make_option(c("-i", "--input"), 
              type = "character", 
              default = NULL, 
              help = "name of input .zip file or folder", 
              metavar="character"), 
  make_option(c("-o", "--output"), 
              type = "character", 
              default = "output", 
              help = "name of subfolder for temporary output [default = %default]", 
              metavar = "character"), 
  make_option(c("-p", "--pivot"), 
              type = "character", 
              default = "TRUE", 
              help = "flag for pivoting data frame to long format [default = %default]", 
              metavar = "character"), 
  make_option(c("-s", "--step"), 
              type = "integer", 
              default = "1000000", 
              help = "step size for pivot chunks [default = %default]", 
              metavar = "integer"), 
  make_option(c("-f", "--filter"), 
              type = "character", 
              default = "FALSE", 
              help = "flag for filtering out NAs (pivoted data will always be filtered) [default = %default]", 
              metavar = "character"), 
  make_option(c("-v", "--variable"), 
              type = "character", 
              default = NULL, 
              help = "when filtering: name of variable that NAs should be filtered from", 
              metavar = "character"), 
  make_option(c("-g", "--geotiff"), 
              type = "character", 
              default = "FALSE", 
              help = "flag for exporting a GeoTIFF raster image based on defined variable (-v / --variable) [default = %default]", 
              metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##---
#opt$input <- "C:/Users/kgreger/Downloads/S3A_OL_2_WFR____20170507T103435_20170507T103735_20170507T124202_0179_017_222_2340_MAR_O_NR_002.SEN3.zip"
#opt$output <- "output"


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied: name of input .zip file or folder.\n", 
       call. = FALSE)
} else if (opt$filter == "TRUE" & is.null(opt$variable)) {
  print_help(opt_parser)
  stop("When filtering is activated a variable name to filter NAs from must be supplied.\n", 
       call. = FALSE)
} else if (opt$geotiff == "TRUE" & is.null(opt$variable)) {
  print_help(opt_parser)
  stop("When export to GeoTIFF is activated a variable name to generate the GeoTIFF from must be supplied.\n", 
       call. = FALSE)
}


## convert boolean command line parameters to actual boolean
pivot <- if (opt$pivot == "TRUE") TRUE else FALSE
filter.nas <- if (opt$filter == "TRUE") TRUE else FALSE
geotiff <- if (opt$geotiff == "TRUE") TRUE else FALSE
input.type <- if (tolower(substr(opt$input, 
                                 nchar(opt$input) - 3, 
                                 nchar(opt$input))) == ".zip") "zip" else "folder"

## create absolute paths from relative paths
opt$input <- getAbsolutePath(opt$input)
path.to.extract.from <- if (substr(dirname(opt$input), 
                                   nchar(dirname(opt$input)) - 3, 
                                   nchar(dirname(opt$input))) == "/") dirname(opt$input) else paste0(dirname(opt$input), "/")
file.to.extract <- gsub("\\.zip$", "", basename(opt$input))
opt$output <- getAbsolutePath(opt$output)
path.to.extract.to <- if (substr(opt$output, 
                                 nchar(opt$output), 
                                 nchar(opt$output)) == "/") substr(opt$output, 
                                                                   1, 
                                                                   nchar(opt$output) - 1) else opt$output

## convert integer command line parameter to actual integer
step.size <- as.integer(opt$step)


## fork for zip vs folder
if (input.type == "zip") {
  ## extract zip file
  unzip(paste0(path.to.extract.from, 
               file.to.extract, 
               ".zip"), 
        exdir = path.to.extract.to, 
        junkpaths = TRUE)
  input.path <- path.to.extract.to
} else {
  input.path <- path.to.extract.from
}


## check for filter variable
## extract zip file
if (opt$filter == "TRUE" | opt$geotiff == "TRUE") {
  nc.files <- gsub("\\.nc", "", list.files(input.path, pattern = ".nc"))
  if (tolower(opt$filter) %in% nc.files) {
    variable.to.filter.null.values.by <- opt$filter
  } else {
    stop("Supplied filter variable was not found in input file.\n", 
         call. = FALSE)
  }
}


## connect to PostgreSQL
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, 
                 dbname = "netcdf",
                 host = "localhost", 
                 port = 5432,
                 user = "postgres", 
                 password = "postgres")


## load geo_coordinates.nc
nc <- nc_open(paste0(input.path, 
                     "/geo_coordinates.nc"))
## extract correct dimensions for this nc file
nc.dims <- c(nc$dim[[1]]$len, nc$dim[[2]]$len)
## generate ids
id <- c(matrix(1:(nc.dims[1] * nc.dims[2]), nrow = nc.dims[1]))
df <- data.frame(id)
## clean up
nc_close(nc)
rm(nc)

stop()


## loop through contents of nc file and extract all variables
files <- list.files(input.path, pattern = ".nc")
for(f in files) {
  ## load nc file
  nc <- nc_open(paste0(input.path, 
                       "/", 
                       f), 
                suppress_dimvals = TRUE)
  
  ## check for valid dimension
  if (nc$ndims > 1) {
    if (nc$dim[[1]]$len != nc.dims[1] | nc$dim[[2]]$len != nc.dims[2]) { next }
  }
  
  cat(paste("\nProcessing layer:", f))
  
  ## loop through ll variables in nc file
  for(v in 1:nc$nvars) {
    ## do stuff
    variable <- data.frame(c(ncvar_get(nc, nc$var[[v]]$name)))
    colnames(variable) <- nc$var[[v]]$name
    df <- cbind(df, variable)
  }
  
  ## clean up
  nc_close(nc)
  rm(nc, variable)
  
}


## fork execution if pivot or not
if (pivot) {
  ## pivot df to long format
  for (i in seq(0, nrow(df), step.size)) {
    cat(paste0("Writing chunk ", i / step.size + 1, "\n"))
    df.temp <- melt(df[i:(step.size + i), ], id = c("id", "latitude", "longitude")) %>% 
      filter(!is.na(.$value))
    if (i == 0) {
      ## write data to PostgreSQL
      dbWriteTable(con,
                   "data_long",
                   value = df.temp,
                   append = FALSE,
                   row.names = FALSE)
    } else {
      ## write data to PostgreSQL
      dbWriteTable(con,
                   "data_long",
                   value = df.temp,
                   append = TRUE,
                   row.names = FALSE)
    }
    rm(df.temp)
  }
} else {
  if (filter.nas) {
    ## filter for NULL values in defined variable
    cat(paste0("Filtering data by '", variable.to.filter.null.values.by, "'\n"))
    df <- df[!is.na(get(paste0("df$", variable.to.filter.null.values.by))), ]
  }
  
  ## write data to PostgreSQL
  cat("Writing data\n")
  dbWriteTable(con,
               "data",
               value = df,
               append = FALSE,
               row.names = FALSE)
}


## export to GeoTIFF
if (geotiff) {
  cat("Writing GeoTIFF\n")
  nc <- nc_open(paste0(input.path, 
                       "/", 
                       variable.to.filter.null.values.by, 
                       ".nc"))
  r <- ncvar_get(nc, nc$var[["CHL_NN"]]$name) %>% 
    t(.) %>% 
    raster(.,
           xmn = min(df$longitude), 
           xmx = max(df$longitude), 
           ymn = min(df$latitude), 
           ymx = max(df$latitude), 
           crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  #plot(r)
  writeRaster(r, 
              paste0(path.to.extract.from, 
                     "/", 
                     variable.to.filter.null.values.by, 
                     ".tif"), 
              format = "GTiff", 
              overwrite = TRUE)
}



## clean up
dbDisconnect(con)
if (input.type == "zip") unlink(path.to.extract.to, 
                                recursive = TRUE)
rm(df)
cat("Done.\n")
