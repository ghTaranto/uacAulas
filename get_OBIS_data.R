# AUTHOR: Gerald H. Taranto

# install_github("robis")
# install.packages(data.table)
library(robis)
library(data.table)

# Download spp of interest
spp_of_interest <- c("Mora moro", "Helicolenus dactylopterus")

# Save date of download
date_string <- gsub("-", "_", Sys.Date())  

# DOWNLOAD AND SAVE DATA
lspp <- list()

for(spp in spp_of_interest){
  lspp[[spp]] <- robis::occurrence(spp)
  # data.table::fwrite(sp_, gsub(" ", "_", paste0(spp, "_", date_string, ".csv")))
}

as.data.table(lspp$`Mora moro`)


lapply(lspp, function(x) as.data.table(x)[, keep_col, with=F])


# SELECT COLUMN AND RBIND CSVs
keep_col <- c("scientificName", "decimalLatitude", "decimalLongitude", "minimumDepthInMeters", "maximumDepthInMeters", "year", "basisOfRecord", "institutionCode")




lspp <- lapply(lspp, function(x) x[, keep_col])

sppDT <- do.call(rbind, lspp)



spp  <- list.files(pattern = "\\.csv")
obis <- do.call("rbind", lapply(spp, fread, select = keep_col))


# SAVE OBIS TABLE
fwrite(obis, paste0("obisAllSpp", date_string, ".csv"))


do.call(rbind, lspp)
