### EVAPORATIVE DEMAND DROUGHT INDEX (EDDI) ###

YEARWEEKS <- 48 #Number of weeks per year
RANGELON <- 100 #Range of longitudes

#' Sequence of weeks over the period analyzed 
#'
#' @param week selected weeks
#' @param ntimes number of weeks available
#' @return sequence of weeks
#' @export

dates_week <- function(week, ntimes){
  return(seq(week, ntimes, by = YEARWEEKS))
}

#' Calculate eddi for a point/pixel (1-month time scale and weekly frequency)
#'
#' @param eto_ij sequence of ETo values for a point/pixel
#' @param scale selected time scale
#' @return sequence of EDDI values
#' @export
calculate_eddi <- function(eto_ij, scale = 1){
  eddi_ij <- as.numeric(array(NA, dim = c(length(eto_ij))))
  realscale = scale * YEARWEEKS/12

  eto_accum <- array(0, dim = length(eto_ij)) #ETo accumulated for each week according to the selected time scale
  eto_accum[1:(realscale-1)] <- NA
  j <- 1
  for (j in 1:realscale){
    eto_accum[realscale:length(eto_accum)] <- eto_accum[realscale:length(eto_accum)] + eto_ij[j:(length(eto_ij)-(realscale-j))]
  }

  week <- 1
  for (week in 1:YEARWEEKS){
    weeks <- dates_week(week, length(eto_ij))

    eto_weeks <- - eto_accum[weeks] #The sign is changed, thus the increases of the ETo are reflected in decreases in the index values 
    
    #Standardization of index values based on Log-logistic distibution (glo)
    eddi_ij[weeks] <- tryCatch({
      if(sum(!is.na(eto_weeks)) > 0){
        b <- samlmu(eto_weeks, nmom = 4, sort.data = TRUE, ratios = TRUE, trim = 0)
        pel_glo <- pelglo(b)
        cdf_glo <- cdfglo(eto_weeks, para = pel_glo)
        qnorm(cdf_glo)
      }else{
        NA
      }
    }, error = function(cond){
      print(paste("ERROR calculate_eddi"))
      return(NA)
    })
  }
  return(eddi_ij)
}

#' Transform an ETo NetCDF into an EDDI NetCDF (1-month time scale and weekly frequency)
#'
#' @param file_name set working directory of ETo NetCDF
#' @param scale selected time scale
#' @return None
#' @export

calculate_eddi_netcdf <- function(file_name, scale = 1) {
  nceto <- nc_open(file_name)
  lons <- ncvar_get(nceto, "lon") #Longitudes
  lats <- ncvar_get(nceto, "lat") #Latitudes
  time <- ncvar_get(nceto,"time") #Times(weeks)

  file_eddi <- gsub(".nc", paste0(scale, "_eddi.nc"), file_name)
  file.copy(file_name, file_eddi)
  nceddi <- nc_open(file_eddi, write = TRUE)

  lon_seqs <- seq(1, length(lons), by = RANGELON)

  lon_seq <- lon_seqs[1]
  for (lon_seq in lon_seqs){
    lon_length <- RANGELON
    if(lon_seq + lon_length > length(lons)){
      lon_length <- length(lons) - lon_seq
    }
    eto <- ncvar_get(nceto, nceto$var[[1]]$name, c(lon_seq, 1, 1), c(lon_length, -1, -1)) 
    eddi <- array(NA, dim = dim(eto))

    print(paste("lon", lon_seq, date()))

      i <- 1
    for (i in 1:lon_length){
      j <- 1
      for (j in 1:length(lats)){ 
        eto_ij <- eto[i, j, ]

        eddi_ij <- calculate_eddi(eto_ij, scale)
        eddi[i, j, ] <- round(eddi_ij, 4)
      }
    }
    #In some exceptions in which the probability may be zero an infinite value can be obtained. 
    #These Infinite values (Inf or -Inf) are replaced by +/-2.58 (corresponding with a return period of 1 in 200 years)
    eddi[is.infinite(eddi)&eddi<0] <- -2.58
    eddi[is.infinite(eddi)&eddi>0] <- 2.58
    ncvar_put(nceddi, nceddi$var[[1]]$name, eddi, c(lon_seq, 1, 1), c(lon_length, -1, -1))
  }
  nc_close(nceto)
  nc_close(nceddi)
}

