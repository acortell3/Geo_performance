


###### FUNCTIONS TO ROTATE GEOS (also for GeomeasuRe)
library(sf)

## x: is spatial object
## orientation: Which side do you want to orient the object to? Options are left, right, custom
## angle: What is the angle to rotate. Prefined options are 90, 180 and 270 degrees.
## The user can decide its custom transformation (option = "custom"), in which case
## it must define the value of theta. Default is 180
## just.info: In case you only want to see the original orientation of the object

## The arguments orientation and angle might seem contradictory. Orientation commands.
## If orientation == right or left, and the object is not oriented in the desired direction,
## the object is rotated 180º. If another value is desired for the rotation, then the
## orientation argument must be set to "custom" and the desired degrees must be put
## in the argument 'angle'. If the desired orientation is the original orientation of the 
## piece a warning will be issued (inconsequential). If this argument has any value other
## than right, left or custom, the function will throw an error.


## Prepare function
orient <- function(x, orientation = "left", angle = 180, just.info = FALSE){
  
  ## Prepare object
  fig <- st_coordinates(x)[,c(1,2)] # Extract coordinates
  bbox <- st_bbox(x) ## extract bbox
  bbox <- data.frame("X" = rep(bbox[c(1,3)],2), "Y" = c(rep(bbox[2],2), rep(bbox[4],2)))
  
  ## Individualise points bbox
  pNE <- bbox[which(bbox[,2] == max(bbox[,2]) & bbox[,1] == max(bbox[,1])),]
  pSE <- bbox[which(bbox[,2] == min(bbox[,2]) & bbox[,1] == max(bbox[,1])),]
  pNW <- bbox[which(bbox[,2] == max(bbox[,2]) & bbox[,1] == min(bbox[,1])),]
  pSW <- bbox[which(bbox[,2] == min(bbox[,2]) & bbox[,1] == min(bbox[,1])),]

  ## Compute distance between Western corner points
  ## Closest NW point in the figure to NW bbox point
  ## Euclidean distance between two points
  xsNW <- (fig[,1]-pNW[,1])^2
  ysNW <- (fig[,2]-pNW[,2])^2
  pgNW <- sqrt(xsNW+ysNW)
  fpNW <- fig[which(pgNW == min(pgNW))[1],] ## Select point
  
  ## Compute distance between Western corner points
  ## Closest SW point in the figure to SW bbox point
  ## Euclidean distance between two points
  xsSW <- (fig[,1]-pSW[,1])^2
  ysSW <- (fig[,2]-pSW[,2])^2
  pgSW <- sqrt(xsSW+ysSW)
  fpSW <- fig[which(pgSW == min(pgSW))[1],] ## Select point
  
  ## Compute distance between the two W corners of the figure
  xWdist <- (fpNW[1]-fpSW[1])^2
  yWdist <- (fpNW[2]-fpSW[2])^2
  Wdist <- sqrt(xWdist+yWdist)
  
  ## Compute distance between Eastern corner points
  ## Closest NE point in the figure to NE bbox point
  ## Euclidean distance between two points
  xsNE <- (fig[,1]-pNE[,1])^2
  ysNE <- (fig[,2]-pNE[,2])^2
  pgNE <- sqrt(xsNE+ysNE)
  fpNE <- fig[which(pgNE == min(pgNE))[1],] ## Select point
  
  ## Compute distance between Western corner points
  ## Closest SE point in the figure to SE bbox point
  ## Euclidean distance between two points
  xsSE <- (fig[,1]-pSE[,1])^2
  ysSE <- (fig[,2]-pSE[,2])^2
  pgSE <- sqrt(xsSE+ysSE)
  fpSE <- fig[which(pgSE == min(pgSE))[1],] ## Select point
  
  ## Compute distance between the two E corners of the figure
  xEdist <- (fpNE[1]-fpSE[1])^2
  yEdist <- (fpNE[2]-fpSE[2])^2
  Edist <- sqrt(xEdist+yEdist)
  
  o <- c()
  
  if (Edist > Wdist){
    o <- "left"
  } else {
    o <- "right"
  }
  
  if (just.info == TRUE){
    return(o)
  } else {
    if (orientation == "custom"){
      if (angle == 90){
        theta <- pi/2
      } else if (angle == 180){
        theta <- pi
      } else if (angle == 270){
        theta <- (3*pi)/2
      } else {
        theta <- angle
      }
    } else if (orientation == "right" | orientation == "left"){
      theta <- pi
    } else {
      stop("Orientation value not admitted. Please choose between right, left or custom")
    }
    
    if (orientation == o){
      warning("The shape is already oriented to the specified orientation")
      return(x)
    } else {

      ## Extract centroid
      centr <- st_coordinates(suppressWarnings(st_centroid(x)))
      centrx <- fig[,1]-centr[1]
      centry <- fig[,2]-centr[2]
      
      ## Perform rotation
      xi <- (centrx*cos(theta))-centry*sin(theta)
      yi <- (centrx*sin(theta))-centry*cos(theta)
      
      ## Assign rotated coordinates to df (putting back to centroid)
      rot_coords <- data.frame("X" = xi+centr[1],
                               "Y" = yi+centr[2])
      
      ## Convert back to spatial object
      res <- st_polygon(list(as.matrix(rot_coords)))
      res <- st_sfc(res, crs = st_crs(x))
      return(res)
      
    }
  }
}



## EXAMPLE
#geo <- st_read(paste0(dsn,"r.1.shp")) ## Load geometric spatial object
#fig <- st_coordinates(geo)[,c(1,2)] ## Convert to coordinate df
#bbox <- st_bbox(geo) ## extract bbox
#bbox <- data.frame("X" = rep(bbox[c(1,3)],2), "Y" = c(rep(bbox[2],2), rep(bbox[4],2)))
#geo_oriented <- orient(geo, orientation = "left")

#plot(geo_oriented)
#plot(geo, add = TRUE)

