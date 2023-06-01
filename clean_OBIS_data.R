library(sf)
library(data.table)
library(ggplot2)
library(ggspatial)


# azoresShape <- st_read("shapeFiles/Azores_wgs84.shp")
azoresEEZ   <- st_read("shapeFiles/AzoresEEZ.shp")
world       <- st_read("shapeFiles/Mapa_Mundo_wgs84.shp")


sppDB_path <- list.files(pattern = "obisAllSpp")
sppDB <- fread(sppDB_path[2])

sppSpatial <- st_as_sf(sppDB, coords = c("decimalLongitude","decimalLatitude"))
st_crs(sppSpatial) <- 4326


# azoresEEZ_crs <- st_crs(azoresEEZ, parameters = T)
# world_crs     <- st_crs(world, parameters = T)
# sppSpatial    <- st_crs(sppSpatial, parameters = T)

# azoresEEZ_crs$epsg
# world_crs$epsg
# sppSpatial$epsg


# N <- nrow(sppSpatial)


# plot( st_geometry(world), col = "grey90", border = "grey10")
# plot( st_geometry(sppSpatial), pch = 16, col = '#ff000040', add = TRUE)
# plot( st_geometry(azoresEEZ), border = "cyan4", lwd=3, add = TRUE)


# Clean spatial
sppAzores <- st_filter(sppSpatial, azoresEEZ)


# plot( st_geometry(world), col = "grey90", border = "grey10")
# plot( st_geometry(sppAzores), pch = 16, col = '#ff000040', add = TRUE)
# plot( st_geometry(azoresEEZ), border = "cyan4", lwd=3, add = TRUE)


table( sppAzores$institutionCode ) 


dev.off()

# plot
# plot( sppAzores["scientificName"], 
#       # xlim = st_bbox(azoresEEZ)[c(1, 3)],
#       # ylim = st_bbox(azoresEEZ)[c(2, 4)],
#       pch = 16, cex = 0.51)
#       # key.width = lcm(7), 
#       # axes = TRUE)
# plot( st_geometry(azoresEEZ), border = "cyan4", lwd=2, add = TRUE)
# 
# plot( st_geometry(azoresShape), col = "grey10")
# 
# 
# legend("bottom")
# 
# legend("bottom", 
#        legend=c("sppAzores"),
#        col=c("black","green","red"),
#        lwd=2,horiz=TRUE)


# CHANGE GGPLOT FONTS
# library(showtext)
# 
# showtext_auto(TRUE)
# 
# fl <- font_files()
# fl[grepl(pattern = "geo", fl$family, ignore.case = TRUE), ]
# 
# font_add("geo_reg", "georgia.ttf")
# font_add("geo_bold","georgiab.ttf")
# font_add("geo_ita", "georgiai.ttf")
# font_add("geo_bita", "georgiaz.ttf")
# 
# 
# sppAzores$scientificName <- as.factor(sppAzores$scientificName)
# 
# ggplot() +
#   geom_sf(data = azoresEEZ, col = "black", size = 1, fill = "transparent", legend = FALSE) + 
#   geom_sf(data = sppAzores, aes(col = scientificName, shape = scientificName), size = 0.80, alpha = 0.5, legend = FALSE) + 
#   # theme_classic() + 
#   theme(#panel.border = element_rect(fill = NA), 
#         legend.position="none", 
#         strip.text.x = element_text(face = "italic", size = 20, family = "geo_bita", color = "white"),
#         axis.text = element_text(face = "italic", size = 12, family = "geo_reg"),
#         panel.grid = element_line(linetype = 2, color = "grey85", size = 0.5),
#         panel.border = element_rect(color = "black", fill = NA, size = 1),
#         panel.background = element_rect(fill = "#f6eee3"), 
#         strip.background = element_rect(color = "black", size = 1, fill = "black")) + 
#   facet_grid(. ~ scientificName)
# 

# DATA CLEAN temporal 
range(sppAzores$year, na.rm = TRUE)
table(sppAzores$scientificName)

sppAzores <- sppAzores[!is.na(sppAzores$year), ]
table(sppAzores$scientificName)

sppAzores <- sppAzores[sppAzores$year >= 2000, ]
table(sppAzores$scientificName)


# dtm 
library(terra)
dtm <- rast("rasterFiles/emodnetDTM.tif")
plot(dtm)


dtm1000 <- dtm
dtm1000[which(values(dtm1000 < -1000))] <- NA
plot(dtm1000)


# prepare models
# sppAzores$scientificName <- gsub(" ", "_", sppAzores$scientificName) # substitute space in species name with _


ext_depth <- extract(dtm1000, sppAzores, cells = TRUE, ID = FALSE)

sppAzores$depth <- ext_depth[,1]
sppAzores$rcell <- ext_depth[,2]


sppAzores <- sppAzores[!is.na(sppAzores$depth), ]


sppAzoresDT <- as.data.table(sppAzores)
sppAzoresDT <- unique(sppAzoresDT, by = c('scientificName','rcell'))
table(sppAzoresDT$scientificName)



sppAzores[!duplicated()]


sppAzores[!dup]
sppAzores[!duplicated(sppAzores[c("scientificName", "rcell")]),]

table(sppAzores$scientificName)


l <- list() 
for(s in unique(sppAzores$scientificName)){
  
  modelDT         <- as.data.table(values(dtm1000))
  modelDT$cell    <- 1:nrow(modelDT)
  modelDT         <- modelDT[!is.na(modelDT$emodnetDTM), ]
  
  modelDT$presAbs[!is.na(modelDT$emodnetDTM)] <- 0
  
  ext <- as.data.table(extract(dtm1000, sppAzores[sppAzores$scientificName == s, ], cells=TRUE, ID = FALSE))
  
  modelDT$presAbs[ modelDT$cell %in% ext$cell ] <- 1
  
  r <- dtm1000
  r[] <- NA
  
  r[modelDT$cell] <- modelDT$presAbs
  plot(r, type = "classes")
  title(s)
  
  modelDT$spp <- s
  
  l[[s]] <- modelDT
}
  

modelDT <- do.call(rbind, l)

modelDT[presAbs == 1, mean(emodnetDTM), by=spp]


sppAzoresDT

ggplot(sppAzoresDT, aes(x = scientificName, y=depth, fill = scientificName)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="grey30", alpha=0.2) + 
  ylab("Depth (m)") + xlab("") + 
  theme_bw() +
  theme(legend.position="none",
        ) +
  scale_fill_manual(values = c("bisque2", "azure2")) # check http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
                                        


sppAzoresDT[scientificName == "Helicolenus dactylopterus", depth]



t.test(sppAzoresDT[scientificName == "Helicolenus dactylopterus", depth], 
       sppAzoresDT[scientificName == "Mora moro", depth])


       







# Generating data
set.seed(2021)
d <- as.data.frame(cbind(rnorm(1:20, 500, 50), c(rep(0, 10), rep(1, 10))))
h <- modelDT[spp=="Helicolenus_dactylopterus" & presAbs == 1, emodnetDTM]
m <- modelDT[spp=="Mora_moro" & presAbs == 1, emodnetDTM]

#Difference in means
original <- diff(tapply(h, m, mean))
mean(outcome[treatment==1])-mean(outcome[treatment==0])

#Permutation test
permutation.test <- function(treatment, outcome, n){
  distribution=c()
  result=0
  for(i in 1:n){
    distribution[i]=diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))
  }
  result=sum(abs(distribution) >= abs(original))/(n)
  return(list(result, distribution))
}

test1 <- permutation.test(treatment, outcome, 10000)
hist(test1[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=original, lwd=3, col="red")

test1[[1]]


#Compare to t-test
t.test(outcome~treatment)







hist(l$Helicolenus_dactylopterus[presAbs == 1, emodnetDTM])
hist(l$Mora_moro[presAbs == 1, emodnetDTM])



fit <- glm(presAbs ~ emodnetDTM,
           family = binomial(link = "logit"),
           data = l$Helicolenus_dactylopterus)


summary(fit)


