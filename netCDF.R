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

## prepare list of files in nc file
nc.files <- unzip(paste0(path.to.extract.from, 
                         file.to.extract, 
                         ".zip"), 
                  list = TRUE) %>% 
  select(Name) %>% 
  filter(substr(Name, nchar(Name) - 2, nchar(Name)) == ".nc") %>% 
  filter(substr(Name, nchar(Name) - 18, nchar(Name)) != "/geo_coordinates.nc")

## connect to PostgreSQL
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, 
                 dbname = "netcdf",
                 host = "localhost", 
                 port = 5432,
                 user = "postgres", 
                 password = "postgres")


## --- this is where the magic happens ---


## load geo_coordinates.nc
unzip(paste0(path.to.extract.from, 
             file.to.extract, 
             ".zip"), 
      files = paste0(file.to.extract, 
                     "/geo_coordinates.nc"), 
      exdir = path.to.extract.to)
nc <- nc_open(paste0(path.to.extract.to, 
              "/", 
              file.to.extract, 
              "/geo_coordinates.nc"))
## extract correct dimensions for this nc file
nc.dims <- c(nc$dim[[1]]$len, nc$dim[[2]]$len)
## extract id and geo coordinates from nc file
id <- c(matrix(1:(nc.dims[1] * nc.dims[2]), nrow = nc.dims[1]))
lat <- c(ncvar_get(nc, "latitude"))
lon <- c(ncvar_get(nc, "longitude"))
df <- data.frame(id, lon, lat)
## clean up
nc_close(nc)
unlink(path.to.extract.to, 
       recursive = TRUE)
rm(nc, lat, lon) #df



## loop through contents of nc file and extract all variables
for(i in 1:length(nc.files[, 1])) {
  cat(paste0("\n", nc.files[i, 1]))
  ## unzip single nc file
  unzip(paste0(path.to.extract.from, 
               file.to.extract, 
               ".zip"), 
        files = nc.files[i, 1], 
        exdir = path.to.extract.to)
  ## load nc file
  nc <- nc_open(paste0(path.to.extract.to, 
                       "/", 
                       nc.files[i, 1]), 
                suppress_dimvals = TRUE)
  ## do stuff (extract, pivot, ...)
  ## check for valid dimension
  if (nc$ndims > 1) {
    if (nc$dim[[1]]$len != nc.dims[1] | nc$dim[[2]]$len != nc.dims[2]) { next }
  }
  
  ## loop through ll variables in nc file
  for(v in 1:nc$nvars) {
    ## do stuff
    variable <- data.frame(c(ncvar_get(nc, nc$var[[v]]$name)))
    colnames(variable) <- nc$var[[v]]$name
    df <- cbind(df, variable)
  }
  
  ## clean up
  nc_close(nc)
  unlink(path.to.extract.to, 
         recursive = TRUE)
  rm(nc, variable)
  
}

## write data to PostgreSQL
dbWriteTable(con,
             "data",
             value = df,
             append = FALSE,
             row.names = FALSE)

## clean up
rm(df)
