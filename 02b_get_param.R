get_nn <- function(i,yaiy,nearest_inds,weights) {
  
  # i is the pixel index
  
  # yaiy is a matrix of the values you want with the reference values
  # stacked on top of the target (NA) values that need to be filled in
  
  # nearest_inds is the vector of indices of the reference values that
  # are closest to the current pixel in sst and chla values
  
  # weights are how to weight each of the nearest indices, based on the
  # difference between days of the pixel and the reference value
  
  # Weights of the k nearest LUT records (nearest based on SST and CHLA),
  # where the weights were calculated above based on the difference between
  # the record doy (day of year) number and the satellite doy number.
  # as.numeric(nearest_inds[i,]) gives the k numeric indices of the
  # nearest neighbours in the LUT. Each have a weight assigned to them.
  w <- weights[as.numeric(nearest_inds[i,])]
  
  if (sum(w) <= 0.1) {#0.0001) {
    return(rep(NA,4))
  } else {
    # Get the matrix contain the k nearest neighbours (rows) and the
    # y variables to be found based on nearest neighbours (columns).
    nn_mat <- yaiy[as.numeric(nearest_inds[i,]),]
    # Find the mean of each column, weighting each of the k values (rows)
    # by their day number.
    nn <- apply(nn_mat, 2, weighted.mean, w=w, na.rm=T)
    return(nn)
  }
  
}

get_param <- function(logchla, sst, month, day, ref_tab, param_type, num_cl) {
    
    # Compute biomass profile or PI curve parameters for a set of records/pixels,
    # given input reference parameters.
    
    # chla and sst = vectors of logged chlorophyll-a and (not logged) sea surface
    #                temperatures that need parameters
    # ref_tab = table of reference chla/temperature pairings and their associated
    #           parameter values
    # month = month of the input chla and sst
    # param_type = either "bp" or "pi" (biomass profile or photosynthesis-irradiance curve)
    # num_cl = number of clusters to use when parallel processing
    
    library(yaImpute)
    library(parallel)
    num_cl <- min(detectCores()-1, num_cl)
    
    # How many nearest neighbours to use in the computation?
    k = 10
    
    
    
    # logchla <- logchla * 0.25
    # sst <- (sst - 13) * 0.06666667
    
    
    
    
    # FROM GEORGE'S SCRIPT: nn_pre.f90
    # check date: don't use fall data for spring
    # spring is May -- June (dayno=121--182),
    # fall is Sept. 15 -- Nov 15 (dayno=258--319)
    iday_s <- as.numeric(ref_tab[,"dn"]) # vector of days from reference table
    iday <- day # satellite image day
    bad_ind <- (121 <= iday & iday <= 182 & 258 <= iday_s & iday_s <= 319) | # satellite=spring and insitu=fall
               (258 <= iday & iday <= 319 & 121 <= iday_s & iday_s <= 182)   # satellite=fall and insitu=spring
    ref_tab <- ref_tab[!bad_ind,]
    
    
    # Weight of each possible match, based on the difference between the satellite
    # day and reference table days.
    weights = apply(cbind(abs(iday - iday_s[!bad_ind]),abs(iday + 365 - iday_s[!bad_ind])),1,min)
    weights = (1 - (weights / 183) )^4 # gaussian?
    
    
    
    
    # # only use chl in ref table that are between chl=0 and chl=chlmax as computed below?
    # if (zm < 0) {
    #   chlmax <- 32
    # } else {
    #   if (zm > 250) {
    #     zstar <- 250
    #   } else {
    #     zstar <- zm
    #   }
    #   hh <- h / (sigma ( sqrt(2*pi)))
    #   chlmax <- 32 * (B0 + hh * exp(-0.5 * (zm / sigma)^2)) / (B0 + hh * exp(-0.5 * ((zstar-zm) / sigma)^2))
    # }
    
    
    
    
    # Get chla and temperature names
    chlname <- colnames(ref_tab)[grepl("chl",colnames(ref_tab))]
    tempname <- colnames(ref_tab)[grepl("temp",colnames(ref_tab))]
    
    
    # Get temp and log(chla) data from parameter table.
    yaix <- ref_tab[,c(tempname,chlname)]
    yaix[,2] <- log(yaix[,2])
    colnames(yaix) <- c("temp","logchl")

    # Get temp and log(chla) data from satellite image.
    yaix_sat <- data.frame(cbind(sst,logchla),stringsAsFactors=F)
    colnames(yaix_sat) <- c("temp","logchl")

    # Add row for input temperature and log(chla).
    yaix <- rbind(yaix,yaix_sat)

    # Get parameters corresponding to temp and chla in table.
    # Note that you need to use log here to force some values to certain ranges so the results are realistic.
    if (param_type=="bp") {
        num_cols <- 4
        yaiy <- ref_tab[,c("h","sigma","zm","Bo")]
        yaiy[,c(1,2,4)] <- log(1+yaiy[,c(1,2,4)])
        colnames(yaiy) <- c("ln1ph","ln1psigma","zm","ln1pB0")
    } else if (param_type=="pi") {
        num_cols <- 2
        yaiy <- log(ref_tab[,c("pm","alpha")])
        colnames(yaiy) <- c("logpm","logalpha")
    }
    # Make sure rownames are in order, with no gaps (otherwise yai might not
    # retrieve the correct indices for imputation later).
    rownames(yaix) <- 1:nrow(yaix)
    rownames(yaiy) <- 1:nrow(yaiy)
    
    
    # Organize the data and get the k nearest neighbours.
    # This fills in blank "y" corresponding to x, calculated by finding the 10 nearest neighbours to that x.
    # NOTE THAT THIS REMOVES NA VALUES
    yaimethod <- "mahalanobis"#"euclidean"
    yai <- yai(x=yaix,y=yaiy,method=yaimethod,k=k)
    trg_ids <- yai$neiIdsTrgs
    
    # Create clusters for parallel processing.
    cl <- makeCluster(num_cl)
    # Load necessary variables and libraries into cluster.
    # NOTE: By default, clusterExport looks for variables in the global environment.
    #       To change this, you need the envir argument.
    # https://stackoverflow.com/questions/22739876/how-to-export-objects-to-parallel-clusters-within-a-function-in-r
    clusterExport(cl, c("get_nn","yaiy","weights","trg_ids"),envir=environment())
    param_wt <- parLapply(cl, 1:nrow(trg_ids), get_nn, yaiy=yaiy,
                          nearest_inds=trg_ids, weights=weights)
    # Stop parallel processing and return processing power to other operations.
    stopCluster(cl)
    
    param <- matrix(unlist(lapply(param_wt, function(x) {x[1:num_cols]})), ncol=num_cols, byrow=T)
    colnames(param) <- names(param_wt[[1]])
    
    # # For testing
    # param_wt <- list()
    # for (i in 1:nrow(trg_ids)) {
    #   param_wt[[i]] <- get_nn(i=i, yaiy=yaiy, nearest_inds=trg_ids, weights=weights)
    # }
    
    
    
    
    # Impute the mean value of the 10 nearest neighbours from "yai" to the actual value.
    #param <- impute(yai,method="mean",observed=F)
    #param <- impute(yai,method="dstWeighted",observed=F)
    
    
    
    
    # new_param_wt <- t(sapply(1:10, function(i) {apply(yaiy[as.numeric(trg_ids[i,]),],2,weighted.mean,w=weights[as.numeric(trg_ids[i,])],na.rm=T)}))
    # print(param[1:10,])
    # print(param_wt[1:10,])

    
    
    
    # # TEST FOR BP WITH ANN
    # testref <- ref_tab[,c(tempname,chlname)]
    # testref[,2] <- log(testref[,2])
    # colnames(testref) <- c("temp","logchl")
    # testref <- as.matrix(testref)
    # 
    # testtarget <- cbind(sst,logchla)
    # colnames(testtarget) <- c("temp","logchl")
    # #testtarget <- testtarget[1:100,]
    # 
    # test <- ann(ref=testref,target=testtarget,k=k)
    # 
    # # k nearest neighbours, then average it?
    # 
    # annmat <- test$knnIndexDist
    # nearest_ann <- annmat[,1]
    # 
    # best_param <- ref_tab[nearest_ann,c("h","sigma","zm","Bo")]
    # row.names(best_param) <- NULL
    # 
    # # an m x 2k matrix. Each row corresponds to a target point in target and columns
    # # 1:k hold the ref matrix row indices of the nearest neighbors, such that column
    # # 1 index holds the ref matrix row index for the first nearest neighbor and
    # # column k is the k-th nearest neighbor index. Columns k+1:2k hold the Euclidean
    # # distance from the target to each of the k nearest neighbors indexed in columns 1:k.
    
    
    
    
    # Some parameters are transformed before fitting, so here you must transform them back.
    if (param_type=="bp") {
        param[,c(1,2,4)] <- exp(param[,c(1,2,4)]) - 1
        colnames(param) <- c("h","sigma","zm","B0")
    } else if (param_type=="pi") {
        param <- exp(param)
        colnames(param) <- c("pmb","alphab")
    }

    # Get the "target" parameter rows, minus the "reference" rows.
    # (this includes NA values too, for reshaping purposes later (IF YOU USE THIS, YOU NEED NA_IND IN THE FUNCTION INPUT)
    # but if you're using binned ascii files and reshaping later, you shouldn't need this)
    # IF YOU'RE USING THE WEIGHTED CUSTOM FUNCTION...
    #best_param <- data.frame(matrix(nrow=length(logchla),ncol=4),stringsAsFactors=F)
    #best_param[!NA_ind,] <- param[((nrow(ref_tab)+1):nrow(param)),]
    best_param <- param
    colnames(best_param) <- colnames(param)
    # # IF YOU'RE USING THE IMPUTE FUNCTION...
    # best_param <- param[((nrow(ref_tab)+1):nrow(param)),]
    # colnames(best_param) <- colnames(param)
    # row.names(best_param) <- NULL
    
    if (param_type=="bp") {

        # Force B(0) = chla (i.e. surface chla given by biomass paramters = chla from satellite)
        # From George's scripts:
        #   The biomass profile will ultimately need to be scaled so the surface
        #   value B(0) agrees with the satellite value for each pixel:
        #       Bsat = approx B(0) = k * ((~B0) + ((~h)/((~sigma)*sqrt(2*pi)) * exp(-(~zm)^2/(2*(~sigma)^2))
        #   where k is a scale factor that must be determined for each pixel and
        #   ~B0, ~h, ~sigma, and ~zm are the imputed values derived from in situ data.
        #   The production calculations will then use:
        #       B0 = k*(~B0)
        #       h = k*(~h)
        #       sigma= ~sigma
        #       zm = ~zm
        h <- best_param[,1]
        sigma <- best_param[,2]
        zm <- best_param[,3]
        B0 <- best_param[,4]

        insitu_surface_chla <- B0 + (h/(sigma*sqrt(2*pi))) * exp(-0.5*((zm/sigma)^2))

        chla <- exp(logchla)
        newc <- chla/insitu_surface_chla

        best_param[,4] <- best_param[,4] * newc
        best_param[,1] <- best_param[,1] * newc
        B0 <- newc * B0
        h <- newc * h

    }
    
    return(best_param)
    
}
