# install.packages("ncdf4")
# install.packages("dplyr")
# install.packages("reshape2")
# install.packages("RPostgreSQL")
library(ncdf4)
library(dplyr)
library(reshape2)
library(RPostgreSQL)

path.to.extract.from <- "C:/Users/kgreger/Downloads/"
file.to.extract <- "S3A_OL_2_WFR____20170507T103435_20170507T103735_20170507T124202_0179_017_222_2340_MAR_O_NR_002.SEN3"
path.to.extract.to <- "C:/Users/kgreger/Downloads/output"
step.size <- 1000000  ## step size for pivot chunks
pivot <- TRUE  ## should data be pivoted to long format?
filter <- TRUE  ## should NAs be filtered? (pivoted data will always be filtered)
variable.to.filter.null.values.by <- "CHL_NN"  ## if filtering, which variable should NAs be filtered from
export.geotiff <- TRUE


##--- build fork for zip vs folder

## extract zip file
unzip(paste0(path.to.extract.from, 
             file.to.extract, 
             ".zip"), 
      exdir = path.to.extract.to, 
      junkpaths = TRUE)

## connect to PostgreSQL
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, 
                 dbname = "netcdf",
                 host = "localhost", 
                 port = 5432,
                 user = "postgres", 
                 password = "postgres")


## load geo_coordinates.nc
nc <- nc_open(paste0(path.to.extract.to, 
                     "/geo_coordinates.nc"))
## extract correct dimensions for this nc file
nc.dims <- c(nc$dim[[1]]$len, nc$dim[[2]]$len)
## generate ids
id <- c(matrix(1:(nc.dims[1] * nc.dims[2]), nrow = nc.dims[1]))
df <- data.frame(id)
## clean up
nc_close(nc)
rm(nc)


## loop through contents of nc file and extract all variables
files <- list.files(path.to.extract.to, pattern = ".nc")
for(f in files) {
  ## load nc file
  nc <- nc_open(paste0(path.to.extract.to, 
                       "/", 
                       f), 
                suppress_dimvals = TRUE)
  
  ## check for valid dimension
  if (nc$ndims > 1) {
    if (nc$dim[[1]]$len != nc.dims[1] | nc$dim[[2]]$len != nc.dims[2]) { next }
  }
  
  cat(paste("\nProcessing:", f))
  
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
    cat(paste0(i, "\n"))
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
    df <- df[!is.na(get(paste0("df$", variable.to.filter.null.values.by))), ]
  }
  
  ## write data to PostgreSQL
  dbWriteTable(con,
               "data",
               value = df,
               append = FALSE,
               row.names = FALSE)
}



## clean up
dbDisconnect(con)
unlink(path.to.extract.to, 
       recursive = TRUE)
rm(df)

