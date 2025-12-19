# Written By:
#     ___         ___         ___                   ___                  _____        ___         ___               
#    /  /\       /__/\       /  /\      ___        /  /\                /  /::\      /  /\       /  /\        ___   
#   /  /:/       \  \:\     /  /::\    /  /\      /  /:/_              /  /:/\:\    /  /::\     /  /::\      /__/|  
#  /  /:/         \__\:\   /  /:/\:\  /  /:/     /  /:/ /\            /  /:/  \:\  /  /:/\:\   /  /:/\:\    |  |:|  
# /  /:/  ___ ___ /  /::\ /  /:/~/:/ /__/::\    /  /:/ /::\          /__/:/ \__\:|/  /:/  \:\ /  /:/~/:/    |  |:|  
#/__/:/  /  //__/\  /:/\:/__/:/ /:/__\__\/\:\__/__/:/ /:/\:\         \  \:\ /  /:/__/:/ \__\:/__/:/ /:/_____|__|:|  
#\  \:\ /  /:\  \:\/:/__\\  \:\/:::::/  \  \:\/\  \:\/:/~/:/          \  \:\  /:/\  \:\ /  /:\  \:\/:::::/__/::::\  
# \  \:\  /:/ \  \::/     \  \::/~~~~    \__\::/\  \::/ /:/            \  \:\/:/  \  \:\  /:/ \  \::/~~~~   ~\~~\:\ 
#  \  \:\/:/   \  \:\      \  \:\        /__/:/  \__\/ /:/              \  \::/    \  \:\/:/   \  \:\         \  \:\
#   \  \::/     \  \:\      \  \:\       \__\/     /__/:/                \__\/      \  \::/     \  \:\         \__\/
#    \__\/       \__\/       \__\/                 \__\/                             \__\/       \__\/              
#   
# V1 11/5/2025 - 11/14/2025
calculate_stream_depletions <- function(streams,
                                        streams_are_points = FALSE,
                                        stream_id_key = NULL,
                                        wells,
                                        wells_id_key = NULL,
                                        pumping,
                                        model_grid = NULL,
                                        grid_layer_key = 'lay',
                                        grid_stor_coef_key = 'Stor',
                                        grid_transmissivity_key = 'Tr',
                                        subwatersheds = NULL,
                                        influence_radius = NULL,
                                        proximity_criteria = 'whole domain',
                                        apportionment_criteria = 'inverse distance',
                                        geologic_apportionment = FALSE,
                                        analytical_model = 'glover',
                                        depletion_potential_criteria = 'fractional',
                                        sdf_averaging_criteria = 'fractional',
                                        custom_sdf_time = NULL,
                                        data_out_dir = getwd(),
                                        diag_out_dir = getwd(),
                                        suppress_loading_bar = FALSE,
                                        suppress_console_messages = FALSE,
                                        well_stor_coef_key = 'Stor',
                                        well_transmissivity_key = 'Tr',
                                        well_layer_key = 'lay',
                                        stream_transmissivity_key = NULL,
                                        leakance_key = NULL,
                                        lambda_key = NULL,
                                        prec = 80)
{

  ############################################################################################
  ######################################### HELPER FUNCTIONS #################################
  ############################################################################################
  # ==================================================================================================
  # Simple function to concatenate a loading bar
  # ==================================================================================================
  loading_bar <- function(optional_text = '', iter, total, width = 50) 
  {
    if(iter == 1){
      cat('\n')
    }
    # ------------------------------------------------------------------------------------------------
    # percent completion
    pct <- iter / total
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    # how much to fill the bar
    filled <- round(width * pct)
    bar <- paste0(rep("=", filled), collapse = "")
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    # how much is left
    space <- paste0(rep(" ", width - filled), collapse = "")
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    cat(sprintf("\r[%s%s] %3d%% %s",
                bar,
                space,
                round(100*pct),
                optional_text))
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    if (iter == total){
      cat("\n")
    }
    # ------------------------------------------------------------------------------------------------
  }
  # ------------------------------------------------------------------------------------------------
  
  
  
  #===========================================================================================
  # Finds between which indices in the fractional depletions equate to the target depletion
  #===========================================================================================
  calculate_custom_sdf_time <- function(average_fractional_depletions,
                                        target = custom_sdf_time)
  {
    #-------------------------------------------------------------------------------
    # find which index the target is between
    less_than_indices <- average_fractional_depletions < target
    greater_than_indices <- average_fractional_depletions > target
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # if all less than is FALSE then the number must be before the start, return -Inf
    # if all greater than is FALSE then the number must be after the end, return Inf
    if(all(less_than_indices == FALSE)){
      final_time <- -Inf
    } else if(all(greater_than_indices == FALSE)){
      final_time <- Inf
    } else {
      final_less_than_index <- tail(which(less_than_indices == TRUE),1)
      final_greater_than_index <- head(which(greater_than_indices == TRUE),1)
      #-------------------------------------------------------------------------------
      # if the last position less than the target is the final entry, return Inf
      # if the first position greater than the target is the first entry, return -Inf
      if(final_less_than_index == length(average_fractional_depletions)){
        final_time <- Inf
      } else if (final_greater_than_index == 1){
        final_time <- -Inf
      } else {
        
        dist_to_target <- target - average_fractional_depletions[final_less_than_index]
        dist_total <- average_fractional_depletions[final_greater_than_index] - 
          average_fractional_depletions[final_less_than_index]
        
        percent_of_total_distance <- dist_to_target/dist_total
        
        final_time <- final_less_than_index + percent_of_total_distance
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    return(final_time)
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  #===========================================================================================
  #
  #===========================================================================================
  calculate_depletion_potential <- function(depletion_potential_criteria = depletion_potential_criteria,
                                            depletions_potential_per_well_total = depletions_potential_per_well_total,
                                            distances = distances,
                                            fracs = fracs,
                                            pumping = pumping)
  {
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(depletion_potential_criteria) == 'Global'){
      average_fractional_depletions <- base::rowMeans(depletions_potential_per_well_total, na.rm = TRUE)
      
    } 
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(depletion_potential_criteria) == 'Fractional'){
      average_fractional_depletions <- sapply(c(1:ncol(pumping)), function(k){
        
        w <-  fracs/mean(fracs, na.rm = T)
        rm <- which(is.na(w) == TRUE &
                      is.nan(w) == FALSE)
        if(length(rm) > 0){
          w <- w[-c(rm)]
        } else {}
        
        if(all(is.nan(w)) == TRUE){
          w <- rep(1,length(w))
        }
        
        x <- depletions_potential_per_well_total[k,]
        x <- x[is.na(x) == FALSE]
        weighted_mean(x = x,
                      w = w, na.rm = T)
      })
      
    } 
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(depletion_potential_criteria) == 'Fractional+Pumping'){
      average_fractional_depletions <- sapply(c(1:ncol(pumping)), function(k){
        
        w <- (pumping[,k] * fracs)/mean(pumping[,k] * fracs, na.rm = T)
        rm <- which(is.na(w) == TRUE &
                      is.nan(w) == FALSE)
        if(length(rm) > 0){
          w <- w[-c(rm)]
        } else {}
        
        if(all(is.nan(w)) == TRUE){
          w <- rep(1,length(w))
        }
        
        x <- depletions_potential_per_well_total[k,]
        x <- x[is.na(x) == FALSE]
        weighted_mean(x = x,
                      w = w, na.rm = T)
      })
    } 
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(depletion_potential_criteria) == 'Distance'){
      
      average_fractional_depletions <- sapply(c(1:ncol(pumping)), function(k){
        w <- 1/unlist(distances)
        rm <- which(is.na(w) == TRUE &
                      is.nan(w) == FALSE)
        if(length(rm) > 0){
          w <- w[-c(rm)]
        } else {}
        
        if(all(is.nan(w)) == TRUE){
          w <- rep(1,length(w))
        }
        
        x <- depletions_potential_per_well_total[k,]
        x <- x[is.na(x) == FALSE]
        weighted_mean(x = x,
                      w = w, na.rm = T)
      })
      
    } 
    # ------------------------------------------------------------------------------------------------
    
    
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(depletion_potential_criteria) == 'Fractional+Distance'){
      
      average_fractional_depletions <- sapply(c(1:ncol(pumping)), function(k){
        
        w <- (fracs)/mean(fracs, na.rm = T)
        rm <- which(is.na(w) == TRUE &
                      is.nan(w) == FALSE)
        if(length(rm) > 0){
          w <- w[-c(rm)] * 1/(unlist(distances)*unlist(distances))
        } else {
          w <- w * 1/unlist(distances)
        }
        
        if(all(is.nan(w)) == TRUE){
          w <- rep(1,length(w))
        }
        
        x <- depletions_potential_per_well_total[k,]
        x <- x[is.na(x) == FALSE]
        weighted_mean(x = x,
                      w = w, na.rm = T)
      })
      
    }
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(depletion_potential_criteria) == 'All'){
      
      average_fractional_depletions <- sapply(c(1:ncol(pumping)), function(k){
        
        w <- (pumping[,k] * fracs)/mean(pumping[,k] * fracs, na.rm = T)
        rm <- which(is.na(w) == TRUE &
                      is.nan(w) == FALSE)
        if(length(rm) > 0){
          w <- w[-c(rm)] * 1/(unlist(distances)*unlist(distances))
        } else {
          w <- w * 1/unlist(distances)
        }
        
        if(all(is.nan(w)) == TRUE){
          w <- rep(1,length(w))
        }
        
        x <- depletions_potential_per_well_total[k,]
        x <- x[is.na(x) == FALSE]
        weighted_mean(x = x,
                      w = w, na.rm = T)
      })
      
    }
    # ------------------------------------------------------------------------------------------------
    
    return(average_fractional_depletions)
  }
  # ------------------------------------------------------------------------------------------------
  
  
  
  #===========================================================================================
  # 
  #===========================================================================================
  calculate_average_sdf <- function(sdf_averaging_criteria = sdf_averaging_criteria,
                                    sdf_vec = Jenk_SDF_per_well_total,
                                    fracs = fracs,
                                    pumping = pumping)
  {
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(sdf_averaging_criteria) == 'Global'){
      average_sdf <- mean(sdf_vec, na.rm = TRUE)
      
    } 
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(sdf_averaging_criteria) == 'Fractional'){
      w <-  fracs/mean(fracs, na.rm = T)
      rm <- which(is.na(w) == TRUE &
                    is.nan(w) == FALSE)
      if(length(rm) > 0){
        w <- w[-c(rm)]
      } else {}
      
      if(all(is.nan(w)) == TRUE){
        w <- rep(1,length(w))
      }
      
      x <- sdf_vec
      x <- x[is.na(x) == FALSE]
      
      average_sdf <- weighted_mean(x = x,
                                   w = w, na.rm = T)
      
    } 
    # ------------------------------------------------------------------------------------------------
    
    # ------------------------------------------------------------------------------------------------
    if(str_to_title(sdf_averaging_criteria) == 'Fractional+Pumping'){
      
      w <- (pumping[,k] * fracs)/mean(pumping[,k] * fracs, na.rm = T)
      rm <- which(is.na(w) == TRUE &
                    is.nan(w) == FALSE)
      if(length(rm) > 0){
        w <- w[-c(rm)]
      } else {}
      
      if(all(is.nan(w)) == TRUE){
        w <- rep(1,length(w))
      }
      
      x <- sdf_vec
      x <- x[is.na(x) == FALSE]
      
      
      average_sdf <- weighted_mean(x = x,
                                   w = w, na.rm = T)
    } 
    # ------------------------------------------------------------------------------------------------
    
    return(average_sdf)
  }
  # ------------------------------------------------------------------------------------------------
  
  
  
  
  #===========================================================================================
  # USAGE:
  # custom_sdf_time <- 0.28
  # custom_sdf_convergence_threshold <- 0.01
  # n_sdf_convergence_tries <- 1000
  # custom_SDF <- custom_sdf_gradient_descent(analytical_model = analytical_model,
  #                                           distance = distance,
  #                                           stor_coef = stor_coef,
  #                                           transmissivity = transmissivity,
  #                                           lambda = lambda,
  #                                           custom_sdf_time = custom_sdf_time)
  # by gradient descent approximates the inverse of the hunt and hantush model
  # within some user specified precision of the desired answer
  #===========================================================================================
  custom_sdf_gradient_descent <- function(analytical_model = analytical_model,
                                          custom_sdf_time = custom_sdf_time,
                                          distance = distance,
                                          stor_coef = stor_coef,
                                          transmissivity = transmissivity,
                                          lambda = NULL,
                                          leakance = NULL)
  {
    #-------------------------------------------------------------------------------
    # necessary variables for the while loop
    r <- 0
    test_time <- 1
    counter <- 0
    mod <- ((distance*distance)*stor_coef)/transmissivity
    history <- c()
    mod_subtract <- mod/10
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # if hunt model inverse is to be approximated
    if(str_to_title(analytical_model) == 'Hunt'){
      while((r < (custom_sdf_time + custom_sdf_convergence_threshold) &
             r > (custom_sdf_time - custom_sdf_convergence_threshold)) == FALSE){
        
        counter <- counter + 1
        
        #-------------------------------------------------------------------------------
        # assemble terms
        z <- Rmpfr::mpfr((sqrt((stor_coef * distance* distance)/
                                 (4*transmissivity*test_time))), prec = prec)
        t1 <- Rmpfr::erfc(z)
        
        
        t2_a <- Rmpfr::mpfr(((lambda*lambda*test_time)/(4*stor_coef*transmissivity)),  prec = prec)
        t2_b <- Rmpfr::mpfr(((lambda*distance)/(2*transmissivity)), prec = prec)
        t2 <- base::exp(t2_a + t2_b)
        
        
        t3_a <- Rmpfr::mpfr((sqrt((lambda*lambda*test_time)/(4*stor_coef*transmissivity))),  prec = prec)
        t3 <- Rmpfr::erfc(t3_a + z)
        r <- as.numeric(t1 - (t2*t3))
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # is it reasonable that an answer can be approached?
        # is transmissivity too low over the given distance that the answer can be approximated
        # given the precision
        if(is.nan(r) == FALSE){
          if(r > (custom_sdf_time + thresh)){
            test_time <- test_time-(mod*(1-erfc(r-custom_sdf_time)))
            history <- append(history, '+')
            if(test_time < 0){
              test_time <- 1
            }
          } else if(r < (custom_sdf_time + thresh)){
            test_time <- test_time+(mod*(1-erfc(custom_sdf_time-r)))
            history <- append(history, '-')
            if(test_time < 0){
              test_time <- 1
            }
          }
          
          #-------------------------------------------------------------------------------
          # if bouncing around the answer but taking too big of step sizes
          if(length(history) >= 10){
            if(paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('+','-'),5), collapse = '')|
               paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('-','+'),5), collapse = '')){
              mod <- mod - mod_subtract
              history <- c()
            }
            
            if(paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('-'),10), collapse = '')|
               paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('+'),10), collapse = '')){
              mod <- mod*2
              mod_subtract <- mod/10
              history <- c()
            }
          }
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # if max tries reached break
          if(counter == n_sdf_convergence_tries){
            test_time <- -9999
            r <- custom_sdf_time
          } else {}
          #-------------------------------------------------------------------------------
        } else {
          test_time <- -9999
          r <- custom_sdf_time
        }
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # if hantush model is to be approximated by gradient descent
    if(str_to_title(analytical_model) == 'Hantush'){
      while((r < (custom_sdf_time + custom_sdf_convergence_threshold) &
             r > (custom_sdf_time - custom_sdf_convergence_threshold)) == FALSE){
        
        counter <- counter + 1
        
        
        #-------------------------------------------------------------------------------
        # assemble terms
        z <- Rmpfr::mpfr((sqrt((stor_coef * distance* distance)/
                                 (4*transmissivity*test_time))),  prec = prec)
        t1 <- Rmpfr::erfc(z)
        
        
        t2_a <- Rmpfr::mpfr(((transmissivity*test_time)/(stor_coef*leakance*leakance)),  prec = prec)
        t2_b <- Rmpfr::mpfr((distance/leakance), prec = prec)
        t2 <- base::exp(t2_a + t2_b)
        
        
        t3_a <- Rmpfr::mpfr((sqrt((transmissivity*test_time)/(stor_coef*leakance*leakance))),  prec = prec)
        t3 <- Rmpfr::erfc(t3_a + z)
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # fill infinite indices with higher precision numbers
        r <- as.numeric(t1 - (t2*t3))
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        # is it reasonable that an answer can be approached?
        # is transmissivity too low over the given distance that the answer can be approximated
        # given the precision
        if(is.nan(r) == FALSE){
          if(r > (custom_sdf_time + thresh)){
            test_time <- test_time-(mod*(1-erfc(r-custom_sdf_time)))
            history <- append(history, '+')
            if(test_time < 0){
              test_time <- 1
            }
          } else if(r < (custom_sdf_time + thresh)){
            test_time <- test_time+(mod*(1-erfc(custom_sdf_time-r)))
            history <- append(history, '-')
            if(test_time < 0){
              test_time <- 1
            }
          }
          
          #-------------------------------------------------------------------------------
          # if bouncing around the answer but taking too big of step sizes
          if(length(history) >= 10){
            if(paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('+','-'),5), collapse = '')|
               paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('-','+'),5), collapse = '')){
              mod <- mod - mod_subtract
              history <- c()
            }
            
            if(paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('-'),10), collapse = '')|
               paste(history[(length(history)-9):length(history)], collapse = '') == paste(rep(c('+'),10), collapse = '')){
              mod <- mod*2
              mod_subtract <- mod/10
              history <- c()
            }
          }
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # if max tries reached break
          if(counter == n_sdf_convergence_tries){
            test_time <- -9999
            r <- custom_sdf_time
          } else {}
          #-------------------------------------------------------------------------------
        } else {
          test_time <- -9999
          r <- custom_sdf_time
        }
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    return(test_time)
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  #===========================================================================================
  # install required packages if not present
  #===========================================================================================
  require_package <- function(pkg)
  {
    if(require(pkg, character.only = TRUE) == FALSE){
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  # ------------------------------------------------------------------------------------------------

  #===========================================================================================
  # weighted mean function to avoid loading the stats library
  #===========================================================================================
  weighted_mean <- function(x, w, na.rm = FALSE)
  {
    # get where both x and weights are real values
    if(na.rm == T){
      valid_pos <- !(is.na(x) == TRUE | is.na(w) == TRUE)
      x <- x[valid_pos]
      w <- w[valid_pos]
    }
    wmean <- sum(x*w)/sum(w)
    return(wmean)
  }
  # ------------------------------------------------------------------------------------------------
  
  #===========================================================================================
  # must be the geometry itself and not the entire sf object
  #===========================================================================================
  Extract_SF_Linestring_Vertices <- function(geometry)
  {
    
    latitude <- geometry[1][[1]]
    latitude <- as.matrix(latitude)[,2]
    
    
    longitude <- geometry[1][[1]]
    longitude <- as.matrix(longitude)[,1]
    
    
    return(list(latitude,
                longitude))
  }
  #-------------------------------------------------------------------------------
  
  
  #===========================================================================================
  # complimentary error function
  # as x approaches 0 erfc returns 1. 
  # for example in glover model where Qs(t) = Qw*erfc(z)
  # where z = sqrt((d^2*S)/(4Tt))
  # as t approaches infinity the expression z approaches 0
  # and erfc returns 1
  # therefore at t = infinity the cumulative depletion will
  # be equal to the pumping rate
  # therefore when 99% of depletions have happened erfc (z) will return 0.99
  #===========================================================================================
  erfc <- function(x)
  {
    p <- 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
    return(p)
  }
  #-------------------------------------------------------------------------------
  
  
  #===========================================================================================
  # inverse complimentary error function
  #===========================================================================================
  erfcinv <- function(y) 
  {
    return(qnorm(1 - y / 2) / sqrt(2))
  }
  #-------------------------------------------------------------------------------
  
  #===========================================================================================
  # called to find stream points impacted by wells within the same watershed
  #===========================================================================================
  find_adjacent_stream_points <- function(wells,
                                          subwatersheds,
                                          stream_points_geometry)
  {
    #-------------------------------------------------------------------------------
    # for each well first check which subwatershed its within
    # then remove all other subwatersheds
    # in paired list check which stream points are also within that subwatershed
    impacted_points <- list()
    impacted_length <- list()
    for(i in 1:nrow(wells)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = '')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # check which subwatershed its within
      well_intersect_indices <- st_intersects(subwatersheds,
                                              wells[i, ])
      rm_empty_intersections <- which(lengths(well_intersect_indices) == 0)
      if(length(rm_empty_intersections) > 0){
        well_intersect_indices <- c(1:length(well_intersect_indices))[-c(rm_empty_intersections)]
      } else {
        well_intersect_indices <- c(1:length(well_intersect_indices))
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if for some reason outside domain then place NA
      # else find all stream points within that subwatershed
      if(length(well_intersect_indices) == 0){
        impacted_points[[i]] <- as.character(NA)
        impacted_length[[i]] <- as.character(NA)
      } else {
        #-------------------------------------------------------------------------------
        # check which subwatershed touches original subwatershed
        subwatershed_touches_indices <- st_touches(subwatersheds,
                                                   subwatersheds[well_intersect_indices, ])
        subwatershed_touches_indices <- which(lengths(subwatershed_touches_indices) == 0)
        subwatershed_touches_indices <- c(1:length(subwatershed_touches_indices))[-c(rm_empty_intersections)]
        
        if(length(subwatershed_touches_indices) > 0){
          well_intersect_indices <- append(well_intersect_indices,
                                           subwatershed_touches_indices)
        } else {}
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # find which streams are within adjacent watersheds
        strm_intersect_indices <- st_intersects(stream_points_geometry,
                                                st_geometry(subwatersheds[well_intersect_indices, ]))
        rm_empty_intersections <- which(lengths(strm_intersect_indices) == 0)
        if(length(rm_empty_intersections) > 0){
          strm_intersect_indices <- c(1:length(strm_intersect_indices))[-c(rm_empty_intersections)]
        } else {
          strm_intersect_indices <- c(1:length(strm_intersect_indices))
        }
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # if there are no streams within the subwatershed
        if(length(strm_intersect_indices) == 0){
          impacted_points[[i]] <- as.character(NA)
          impacted_length[[i]] <- as.character(NA)
        } else {
          impacted_points[[i]] <- as.character(strm_intersect_indices)

          #-------------------------------------------------------------------------------
          # finding impacted length by reach
          all_keys <- as.vector(unlist(st_drop_geometry(stream_points_geometry[strm_intersect_indices,stream_id_key])))
          
          coords <- st_coordinates(stream_points_geometry[strm_intersect_indices, ])
          coords <- cbind(coords, all_keys)
          coords <- as.data.frame(coords)
          colnames(coords) <- c('x','y',stream_id_key)
          
          split_coords <- split(coords, coords[ ,stream_id_key])
          
          lines_list <- lapply(split_coords, function(x) {
            st_linestring(as.matrix(x[, c("x", "y")]))
          })
          lines_geometry <- st_sf(id = names(split_coords),
                                  geometry = st_sfc(lines_list),
                                  crs = st_crs(wells))
          #-------------------------------------------------------------------------------
          
          impacted_length[[i]] <- as.character(round(sum(st_length(lines_geometry)),0))
        }
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # making sure theyre all the same dimension
    max_impacted_n <- max(lengths(impacted_points))
    impacted_points <- lapply(impacted_points,
                              function(x) append(x,
                                                 rep(NA,max_impacted_n - length(x))))
    impacted_points <- do.call(rbind, impacted_points)
    impacted_length <- do.call(rbind, impacted_length)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # formatting output
    w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
    wells$ImpLMet <- as.vector(impacted_length)
    out <- cbind(w_index, impacted_points)
    colnames(out) <- c('wellN',
                       paste0('PN',c(1:max_impacted_n)))
    out <- as.data.frame(out)
    return(list(out,
                wells))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------

  
  
  
  
  #===========================================================================================
  # called to find stream points impacted by wells within the same watershed
  #===========================================================================================
  find_adjacent_and_expanding_stream_points <- function(wells,
                                                        subwatersheds,
                                                        influence_radius,
                                                        stream_points_geometry)
  {
    #-------------------------------------------------------------------------------
    # for each well first check which subwatershed its within
    # then remove all other subwatersheds
    # in paired list check which stream points are also within that subwatershed
    impacted_points <- list()
    impacted_length <- list()
    for(i in 1:nrow(wells)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = '')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # check which subwatershed its within
      well_intersect_indices <- st_intersects(subwatersheds,
                                              wells[i, ])
      rm_empty_intersections <- which(lengths(well_intersect_indices) == 0)
      if(length(rm_empty_intersections) > 0){
        well_intersect_indices <- c(1:length(well_intersect_indices))[-c(rm_empty_intersections)]
      } else {
        well_intersect_indices <- c(1:length(well_intersect_indices))
      }
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # if for some reason outside domain then place NA
      # else find all stream points within that subwatershed
      if(length(well_intersect_indices) == 0){
        impacted_points[[i]] <- as.character(NA)
        impacted_length[[i]] <- as.character(NA)
      } else {
        #-------------------------------------------------------------------------------
        # check which subwatershed touches original subwatershed
        subwatershed_touches_indices <- st_touches(subwatersheds,
                                                   subwatersheds[well_intersect_indices, ])
        rm_empty_intersections <- which(lengths(subwatershed_touches_indices) == 0)
        if(length(rm_empty_intersections) > 0){
          subwatershed_touches_indices <- c(1:length(subwatershed_touches_indices))[-c(rm_empty_intersections)]
        } else {
          subwatershed_touches_indices <- c(1:length(subwatershed_touches_indices))
        }
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # append to the one well already within
        if(length(subwatershed_touches_indices) > 0){
          well_intersect_indices <- append(well_intersect_indices,
                                           subwatershed_touches_indices)
        } else {}
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # take st_union of buffer around well and all the watersheds identified above
        adjacent_expanding_geometry <- st_union(subwatersheds[well_intersect_indices, ],
                                                st_buffer(wells[i, ], influence_radius))
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # find which streams are within adjacent watersheds
        strm_intersect_indices <- st_intersects(stream_points_geometry,
                                                st_geometry(adjacent_expanding_geometry))
        rm_empty_intersections <- which(lengths(strm_intersect_indices) == 0)
        if(length(rm_empty_intersections) > 0){
          strm_intersect_indices <- c(1:length(strm_intersect_indices))[-c(rm_empty_intersections)]
        } else {
          strm_intersect_indices <- c(1:length(strm_intersect_indices))
        }
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # if there are no streams within the subwatershed
        if(length(strm_intersect_indices) == 0){
          impacted_points[[i]] <- as.character(NA)
          impacted_length[[i]] <- as.character(NA)
        } else {
          impacted_points[[i]] <- as.character(strm_intersect_indices)
          
          #-------------------------------------------------------------------------------
          # finding impacted length by reach
          all_keys <- as.vector(unlist(st_drop_geometry(stream_points_geometry[strm_intersect_indices,stream_id_key])))
          
          coords <- st_coordinates(stream_points_geometry[strm_intersect_indices, ])
          coords <- cbind(coords, all_keys)
          coords <- as.data.frame(coords)
          colnames(coords) <- c('x','y',stream_id_key)
          
          split_coords <- split(coords, coords[ ,stream_id_key])
          
          lines_list <- lapply(split_coords, function(x) {
            st_linestring(as.matrix(x[, c("x", "y")]))
          })
          lines_geometry <- st_sf(id = names(split_coords),
                                  geometry = st_sfc(lines_list),
                                  crs = st_crs(wells))
          #-------------------------------------------------------------------------------
          
          impacted_length[[i]] <- as.character(round(sum(st_length(lines_geometry)),0))
        }
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # making sure theyre all the same dimension
    max_impacted_n <- max(lengths(impacted_points))
    impacted_points <- lapply(impacted_points,
                              function(x) append(x,
                                                 rep(NA,max_impacted_n - length(x))))
    impacted_points <- do.call(rbind, impacted_points)
    impacted_length <- do.call(rbind, impacted_length)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # formatting output
    w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
    wells$ImpLMet <- as.vector(impacted_length)
    out <- cbind(w_index, impacted_points)
    colnames(out) <- c('wellN',
                       paste0('PN',c(1:max_impacted_n)))
    out <- as.data.frame(out)
    return(list(out,
                wells))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  #===========================================================================================
  # called to find stream points impacted by wells within an influence radius
  # calculated by Zipper et al (2019) https://doi.org/10.1029/2018WR024403
  # to be two times the maximum distance from any landpoint within the domain to its closest
  # stream segment
  # this base level assumption ensures 1 >= segments for each well
  #===========================================================================================
  find_local_stream_points <- function(wells,
                                       influence_radius,
                                       stream_points_geometry)
  {
    #-------------------------------------------------------------------------------
    # for each well first check which subwatershed its within
    # then remove all other subwatersheds
    # in paired list check which stream points are also within that subwatershed
    impacted_points <- list()
    impacted_length <- list()
    for(i in 1:nrow(wells)){
      
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = '')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # find all stream points within influence radius
      strm_intersect_indices <- st_intersects(stream_points_geometry,
                                              st_buffer(wells[i, ], influence_radius))
      rm_empty_intersections <- which(lengths(strm_intersect_indices) == 0)
      if(length(rm_empty_intersections) > 0){
        strm_intersect_indices <- c(1:length(strm_intersect_indices))[-c(rm_empty_intersections)]
      } else {
        strm_intersect_indices <- c(1:length(strm_intersect_indices))
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if there are no streams within the radius
      # (should not be possible with propper Zipper 2019 implementation)
      if(length(strm_intersect_indices) == 0){
        impacted_points[[i]] <- as.character(NA)
        impacted_length[[i]] <- as.character(NA)
      } else {
        impacted_points[[i]] <- as.character(strm_intersect_indices)
        
        #-------------------------------------------------------------------------------
        # finding impacted length by reach
        all_keys <- as.vector(unlist(st_drop_geometry(stream_points_geometry[strm_intersect_indices,stream_id_key])))
        
        coords <- st_coordinates(stream_points_geometry[strm_intersect_indices, ])
        coords <- cbind(coords, all_keys)
        coords <- as.data.frame(coords)
        colnames(coords) <- c('x','y',stream_id_key)
        
        split_coords <- split(coords, coords[ ,stream_id_key])
        
        lines_list <- lapply(split_coords, function(x) {
          st_linestring(as.matrix(x[, c("x", "y")]))
        })
        lines_geometry <- st_sf(id = names(split_coords),
                                geometry = st_sfc(lines_list),
                                crs = st_crs(wells))
        #-------------------------------------------------------------------------------
        
        impacted_length[[i]] <- as.character(round(sum(st_length(lines_geometry)),0))
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # making sure theyre all the same dimension
    max_impacted_n <- max(lengths(impacted_points))
    impacted_points <- lapply(impacted_points,
                              function(x) append(x,
                                                 rep(NA,max_impacted_n - length(x))))
    impacted_points <- do.call(rbind, impacted_points)
    impacted_length <- do.call(rbind, impacted_length)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # formatting output
    w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
    wells$ImpLMet <- as.vector(impacted_length)
    out <- cbind(w_index, impacted_points)
    colnames(out) <- c('wellN',
                       paste0('PN',c(1:max_impacted_n)))
    out <- as.data.frame(out)
    return(list(out,
                wells))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  #===========================================================================================
  # all wells impacted by all points
  #===========================================================================================
  find_whole_domain_points <- function(wells,
                                       stream_points_geometry)
  {
    #-------------------------------------------------------------------------------
    # for each well first check which subwatershed its within
    # then remove all other subwatersheds
    # in paired list check which stream points are also within that subwatershed
    impacted_points <- list()
    impacted_length <- list()
    for(i in 1:nrow(wells)){
      
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = '')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # finding impacted length by reach
      all_keys <- as.vector(unlist(st_drop_geometry(stream_points_geometry[ ,stream_id_key])))
      
      coords <- st_coordinates(stream_points_geometry)
      coords <- cbind(coords, all_keys)
      coords <- as.data.frame(coords)
      colnames(coords) <- c('x','y',stream_id_key)
      
      split_coords <- split(coords, coords[ ,stream_id_key])
      
      lines_list <- lapply(split_coords, function(x) {
        st_linestring(as.matrix(x[, c("x", "y")]))
      })
      lines_geometry <- st_sf(id = names(split_coords),
                              geometry = st_sfc(lines_list),
                              crs = st_crs(wells))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      impacted_points[[i]] <- as.character(c(1:nrow(stream_points_geometry)))
      impacted_length[[i]] <- as.character(round(sum(st_length(lines_geometry)),0))
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # making sure theyre all the same dimension
    max_impacted_n <- max(lengths(impacted_points))
    impacted_points <- lapply(impacted_points,
                              function(x) append(x,
                                                 rep(NA,max_impacted_n - length(x))))
    impacted_points <- do.call(rbind, impacted_points)
    impacted_length <- do.call(rbind, impacted_length)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # formatting output
    w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
    wells$ImpLMet <- as.vector(impacted_length)
    out <- cbind(w_index, impacted_points)
    colnames(out) <- c('wellN',
                       paste0('PN',c(1:max_impacted_n)))
    out <- as.data.frame(out)
    return(list(out,
                wells))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  #===========================================================================================
  # convert everything to the projection of the wells
  #===========================================================================================
  ensure_projections <- function(wells,
                                 geometry_list)
  {
    new_list <- list()
    for(i in 1:length(geometry_list)){
      new_list[[i]] <- st_transform(geometry_list[[i]],
                                    st_crs(wells))
    }
    return(new_list)
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  #===========================================================================================
  # For each well, find the closest point on the reaches that it impacts
  #===========================================================================================
  find_closest_points_per_segment <- function(wells,
                                              stream_points_geometry,
                                              stream_id_key){
    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s',
                              'For each well, finding the closest point on the reaches that it impacts'),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # blank data to be filled
    blank_matrix <- base::matrix(nrow = nrow(wells),
                                 ncol = nrow(streams))
    blank_matrix <- as.data.frame(blank_matrix)
    colnames(blank_matrix) <- c(1:ncol(blank_matrix))
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # get the closest point on each stream segment associated with each well
    stats <- list()
    for(i in 1:nrow(wells)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = 'Closest point by reach')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      imp_pts <- as.vector(unlist(impacted_points[i, ])) # get all impacted points for this well
      imp_pts <- as.numeric(imp_pts)[-c(1)] # get rid of well number
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if there are no impacted points there can not be a closest one
      # additionally, trying to index the stream points data frame by a NA, or numeric(0)
      # would cause an error
      if(all(is.na(imp_pts)) == TRUE){
        # pass
        stats[[i]] <- NA
      } else {
        rm_NA_index <- which(is.na(imp_pts) == TRUE)
        if(length(rm_NA_index) > 0){
          imp_pts <- imp_pts[-c(rm_NA_index)] # removing NA indexes
        }
        
        #-------------------------------------------------------------------------------
        # for each stream reach impacted by the well, what is the closest point on that segment
        stream_points_subset <- stream_points_geometry[imp_pts, ]
        stream_points_subset_inner <- split(stream_points_subset,
                                            st_drop_geometry(stream_points_subset[ ,stream_id_key]))
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        closest_points <- lapply(stream_points_subset_inner,
                                 function(x){
                                   rownames <- row.names(x)
                                   closest_index <- rownames[which.min(st_distance(wells[i, ], x))]
                                   as.numeric(closest_index)
                                 })
        stats[[i]] <- as.vector(unlist(closest_points))
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        blank_matrix[i, as.numeric(colnames(blank_matrix)) %in% as.numeric(names(closest_points))] <- 
          closest_points
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    max_closest_n <- max(lengths(stats))
    mean_closest_n <- round(mean(lengths(stats)),0)
    median_closest_n <- round(median(lengths(stats)),0)
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s %s',
                              'Max | Mean | Median number of reaches effected per well: ',
                              paste(max_closest_n,
                                    '|',
                                    mean_closest_n,
                                    '|',
                                    median_closest_n)),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    return(blank_matrix)
  }
  #-------------------------------------------------------------------------------
  
  
  
  #===========================================================================================
  # Apportion depletions by inverse distance method
  #===========================================================================================
  inverse_distance_apportionment <- function(power,
                                             wells,
                                             closest_points_per_segment,
                                             stream_points_geometry)
  {
    #-------------------------------------------------------------------------------
    # inverse distance equation in zipper
    # https://doi.org/10.1029/2018WR024403
    equation <- function(well,
                         power = power,
                         closest_points_per_well,
                         stream_points_geometry = stream_points_geometry)
    {
      #-------------------------------------------------------------------------------
      closest_points <- stream_points_geometry[closest_points_per_well, ]
      reaches <- st_drop_geometry(stream_points_geometry[closest_points_per_well,stream_id_key])
      reaches <- as.vector(unlist(reaches))
      
      dists <- as.vector(st_distance(well,
                                     st_geometry(closest_points)))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if geology is not to be considered simply apportion by inverse distance
      if(geologic_apportionment == FALSE){
        numerator <- 1/(dists**power)
        denominator <- sum(1/(dists**power))
      } 
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if apportioned by inverse distance and geology check whether model grid is used
      # or just information about the transmissivity of the aquifer at the stream and the
      # well
      if(geologic_apportionment == TRUE){

        #-------------------------------------------------------------------------------
        # if no model grid use information about aquifer at streams and wells
        if(is.null(model_grid) == TRUE){
          TR_prime <- as.vector(unlist(st_drop_geometry(stream_points_geometry[closest_points_per_well,
                                                                               stream_transmissivity_key])))
          
          TR <- as.vector(unlist(st_drop_geometry(well[,well_transmissivity_key])))
          
          numerator <- (1/(dists**power))*(TR_prime/TR)
          denominator <- sum((1/(dists**power))*(TR_prime/TR))
        } 
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # if model grid is present draw a line between the well and each stream segment
        # take the aquifer properties along that line
        if(is.null(model_grid) == FALSE){
          
          TR <- list()
          for(i in 1:length(closest_points_per_well)){
            #-------------------------------------------------------------------------------
            # make line between well and stream
            m <- matrix(c(st_coordinates(well)[,1],
                          st_coordinates(well)[,2],
                          st_coordinates(stream_points_geometry[closest_points_per_well[i], ])[,1],
                          st_coordinates(stream_points_geometry[closest_points_per_well[i], ])[,2]),
                        ncol = 2,
                        byrow = TRUE)
            line <- st_sf(st_sfc(st_linestring(m), crs = st_crs(wells)))
            st_geometry(line) <- 'geometry'
            #-------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------
            # selecting correct gridcells
            grid_layers <- as.vector(unlist(st_drop_geometry(model_grid[,grid_layer_key])))
            well_layers <- as.vector(unlist(st_drop_geometry(well[,well_layer_key])))
            gr <- model_grid[grid_layers %in% well_layers, ]
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # intersect line between well and stream with the grid
            int <- st_intersects(gr, line$geometry)
            grid_inds <- c(1:length(int))
            rm <- which(lengths(int) == 0)
            if(length(rm) > 0){
              grid_inds <- grid_inds[-c(rm)]
            } else {}
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            if(length(grid_inds) > 0){
              #-------------------------------------------------------------------------------
              # getting average properties from grid
              gr <- gr[grid_inds, ]
              transmissivity <- mean(as.vector(unlist(st_drop_geometry(gr[,grid_transmissivity_key]))), na.rm = T)
              TR[[i]] <- transmissivity
              #-------------------------------------------------------------------------------
            } else {
              TR[[i]] <- NA
            }
            #-------------------------------------------------------------------------------
          }
          #-------------------------------------------------------------------------------
          TR <- as.vector(unlist(TR))
          TR[is.nan(TR) == TRUE] <- 0
          TR[is.na(TR) == TRUE] <- 0
          numerator <- (1/(dists**power))*(TR/max(TR))
          denominator <- sum((1/(dists**power))*(TR/max(TR)))
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      fractions_of_depletions <- numerator/denominator
      
      return(list(fractions_of_depletions,
                  reaches))
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # use inverse distance weighting to assign fractions of depletions
    fractions_of_depletions <- as.data.frame(base::matrix(nrow = nrow(wells),
                                                          ncol = nrow(streams)))
    colnames(fractions_of_depletions) <- c(1:ncol(fractions_of_depletions))
    reaches <- as.data.frame(base::matrix(nrow = nrow(wells),
                                          ncol = nrow(streams)))
    colnames(reaches) <- c(1:ncol(reaches))
    for(i in 1:nrow(wells)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = 'Apportioning depletions')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      closest_points_per_well <- closest_points_per_segment[i,-c(1)]
      reference <- as.vector(unlist(closest_points_per_well))
      closest_points_per_well <- as.vector(unlist(closest_points_per_well))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # are there any depletions to apportion? if yes continue, if not append NA
      if(all(is.na(closest_points_per_well)) == FALSE){
        
        
        #-------------------------------------------------------------------------------
        # are there any NAs to remove to avoid indexing the stream points by NA
        rm_indices <- which(is.na(closest_points_per_well) == TRUE)
        if(length(rm_indices) > 0){
          #-------------------------------------------------------------------------------
          closest_points_per_well <- closest_points_per_well[-c(rm_indices)] %>%
            as.vector() %>% unlist()
          #-------------------------------------------------------------------------------
        } else {}
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # assign fractions
        output <- equation(well = wells[i, ],
                           power = power,
                           closest_points_per_well = closest_points_per_well,
                           stream_points_geometry = stream_points_geometry)
        fractions_of_depletions[i,
                                is.na(reference) == FALSE] <- output[[1]]
        reaches[i,
                is.na(reference) == FALSE] <- output[[2]]
        #-------------------------------------------------------------------------------
        
      } else {
        # pass
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
  
    
    #-------------------------------------------------------------------------------
    # writeout statistics
    max_dep <- max(fractions_of_depletions, na.rm = T)
    wm <- which(fractions_of_depletions == max_dep, arr.ind = TRUE)
    
    reach_max_dep <- reaches[wm]
    well_max_dep <- as.vector(unlist(st_drop_geometry(wells[wm[,1],wells_id_key])))
    #-------------------------------------------------------------------------------
    

    #-------------------------------------------------------------------------------
    # format writeout pt 2
    w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
    fractions_of_depletions <- cbind(w_index, fractions_of_depletions)
    fractions_of_depletions <- as.data.frame(fractions_of_depletions)
    colnames(fractions_of_depletions) <- c('wellN',
                                           paste0('RN',
                                                  1:(ncol(fractions_of_depletions)-1)))
    reaches <- cbind(w_index, reaches)
    reaches <- as.data.frame(reaches)
    colnames(reaches) <- c('wellN',
                            paste0('RN',
                                   1:(ncol(reaches)-1)))
    #-------------------------------------------------------------------------------


    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s %s',
                              'Max apportioned depletion fraction: ',
                              paste(round(max_dep,4),
                                    'for reach',
                                    paste(reach_max_dep, collapse = ','),
                                    'for well',
                                    paste(well_max_dep, collapse = ','))),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    return(list(fractions_of_depletions,
                reaches))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  #===========================================================================================
  # Apportion depletions by thiessen polygon method
  # method taken from Zipper (2018) https://doi.org/10.1029/2018WR022707
  #===========================================================================================
  thiessen_polygon_apportionment <- function(wells,
                                             closest_points_per_segment,
                                             stream_points_geometry)
  {
    #-------------------------------------------------------------------------------
    fractions_of_depletions <- as.data.frame(base::matrix(nrow = nrow(wells),
                                                          ncol = nrow(streams)))
    colnames(fractions_of_depletions) <- c(1:ncol(fractions_of_depletions))
    reaches <- as.data.frame(base::matrix(nrow = nrow(wells),
                                          ncol = nrow(streams)))
    colnames(reaches) <- c(1:ncol(reaches))
    for(i in 1:nrow(wells)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = 'Apportioning Depletions')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # find closest points on stream reaches for that well
      closest_points_subset <- closest_points_per_segment[i,-c(1)]
      rm <- which(is.na(closest_points_subset))
      if(length(rm) > 0){
        reference <- closest_points_subset
        closest_points_subset <- as.vector(unlist(closest_points_subset[-c(which(is.na(closest_points_subset)))]))
      } else {
        reference <- closest_points_subset
        closest_points_subset <- as.vector(unlist(closest_points_subset))
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if there are closest points continue
      if(length(closest_points_subset) > 0){
        #-------------------------------------------------------------------------------
        # generate thiessen polygons out of the closest points
        closest_stream_points <- stream_points_geometry[closest_points_subset, ]
        closest_stream_points$key <- paste(st_coordinates(closest_stream_points)[ ,1],
                                           st_coordinates(closest_stream_points)[ ,2],
                                           sep = '_')
        closest_stream_points <- closest_stream_points[ ,'key']
        closest_stream_points <- unique(closest_stream_points)
        closest_stream_points$BUFF <- rep(NA, nrow(closest_stream_points))
        keys <- closest_stream_points$key
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        # if there are enough entities to make Thiessen polygons out of
        if(length(closest_points_subset) > 1){
          #-------------------------------------------------------------------------------
          # generate thiessen polygons out of the closest points plus the well
          wells_tmp <- wells[i, ]
          wells_tmp$BUFF <- rep(NA,nrow(wells_tmp))
          wells_tmp <- wells_tmp[ ,'BUFF']
          closest_stream_points <- closest_stream_points[ ,'BUFF']
          closest_stream_points_plus_well <- rbind(wells_tmp, closest_stream_points)
          
          closest_voronoi_plus_well <- st_voronoi(st_union(st_geometry(closest_stream_points_plus_well)))
          closest_voronoi_plus_well <- st_collection_extract(closest_voronoi_plus_well,'POLYGON')
          closest_voronoi_plus_well <- st_sf(geometry = closest_voronoi_plus_well, crs = st_crs(closest_voronoi_plus_well))
          #-------------------------------------------------------------------------------
          
          
          #-------------------------------------------------------------------------------
          envelope <- st_bbox(closest_voronoi_plus_well)
          envelope <- st_as_sfc(envelope)
          # st_geometry(envelope) <- 'geometry'
          closest_voronoi <- st_voronoi(st_union(st_geometry(closest_stream_points)),
                                        envelope = envelope)
          closest_voronoi <- st_collection_extract(closest_voronoi,'POLYGON')
          closest_voronoi <- st_sf(geometry = closest_voronoi, crs = st_crs(closest_voronoi))
          closest_voronoi$key <- keys
          #-------------------------------------------------------------------------------
          
          
          #-------------------------------------------------------------------------------
          # which thiessen polygon does the well intersect
          wells_intersection <- st_intersects(closest_voronoi_plus_well, wells[i, ])
          rm <- which(lengths(wells_intersection) == 0)
          if(length(rm) > 0){
            wells_intersection <- c(1:length(wells_intersection))[-c(rm)]
          } else {
            wells_intersection <- c(1:length(wells_intersection))
          }
          wells_intersection <- closest_voronoi_plus_well[wells_intersection, ]
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get areas
          well_voronoi_area <- as.numeric(st_area(wells_intersection))
          closest_voronoi_intersection <- st_intersection(st_geometry(closest_voronoi),
                                                          st_geometry(wells_intersection))
          closest_voronoi_intersection <- st_sf(closest_voronoi_intersection)
          st_geometry(closest_voronoi_intersection) <- 'geometry'
          closest_voronoi_intersection$key <- closest_voronoi$key
          voronoi_intersected_areas <- as.numeric(st_area(closest_voronoi_intersection))
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # find which intersections get what depletion fraction
          if(geologic_apportionment == FALSE){
            apportioned_depletions <- voronoi_intersected_areas/well_voronoi_area
          } 
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # if geology is to be considered, then see whether model_grid argument passed
          if(geologic_apportionment == TRUE) {
            
            #-------------------------------------------------------------------------------
            # if model grid isnt passed, use information of aquifer at well and at stream location
            if(is.null(model_grid) == TRUE){
              TR <- as.vector(unlist(st_drop_geometry(wells[i,well_transmissivity_key])))
              TR_prime <- as.vector(unlist(st_drop_geometry(closest_stream_points[,stream_transmissivity_key])))
              apportioned_depletions <- (voronoi_intersected_areas*(TR_prime/TR))/sum((voronoi_intersected_areas*(TR_prime/TR)))
            } 
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # if model grid is passed then intersect it with the voronoi polygons and take the
            # properties of those gridcells
            if(is.null(model_grid) == FALSE){
              
              #-------------------------------------------------------------------------------
              # selecting correct gridcells
              grid_layers <- as.vector(unlist(st_drop_geometry(model_grid[,grid_layer_key])))
              well_layers <- as.vector(unlist(st_drop_geometry(wells[i,well_layer_key])))
              gr <- model_grid[grid_layers %in% well_layers, ]
              #-------------------------------------------------------------------------------
              
              #-------------------------------------------------------------------------------
              # if well takes water from stream
              if(length(voronoi_intersected_areas) > 0){
                #-------------------------------------------------------------------------------
                # if polygon doesnt intersect any gridcells then TR[[j]] will be set to mean(numeric(0))
                # which becomes NaN, which is then caught and set to 0, representing no flow
                TR <- list()
                for(j in 1:nrow(closest_voronoi_intersection)){
                  #-------------------------------------------------------------------------------
                  # making st_intersection and assigning TR without generating attribute warning
                  tr <- st_intersection(st_geometry(gr),
                                        closest_voronoi_intersection$geometry[j])
                  tr <- st_sf(tr)
                  st_geometry(tr) <- 'geometry'
                  a <- which(lengths(st_intersects(gr,tr)) != 0)
                  tr[,grid_transmissivity_key] <- as.vector(unlist(st_drop_geometry(gr[a,grid_transmissivity_key])))
                  #-------------------------------------------------------------------------------
                  
                  #-------------------------------------------------------------------------------
                  # ensure there are no overlapping geometries
                  # commented out for now while deciding what to do with multiperf wells
                  # where gridcells of separate layers are perfectly overlapping
                  # tr <- tr[!duplicated(tr$geometry), ]
                  tr <- as.vector(unlist(st_drop_geometry(tr[,grid_transmissivity_key])))
                  TR[[j]] <- mean(tr, na.rm = T)
                  #-------------------------------------------------------------------------------
                }
                TR <- as.vector(unlist(TR))
                TR[is.nan(TR) == TRUE] <- 0
                TR[is.na(TR) == TRUE] <- 0
                apportioned_depletions <- (voronoi_intersected_areas*(TR/max(TR)))/sum((voronoi_intersected_areas*(TR/max(TR))))
                #-------------------------------------------------------------------------------
              } else {
                apportioned_depletions <- rep(NA, nrow(closest_voronoi))
                closest_voronoi_intersection <- closest_voronoi
              }
              #-------------------------------------------------------------------------------
            }
            #-------------------------------------------------------------------------------
          }
          apportioned_depletions <- data.frame(dep = apportioned_depletions,
                                               key = closest_voronoi_intersection$key)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get which voronoi polygons these depletions should be assigned to
          surrounding_voronoi <- st_intersects(closest_voronoi, closest_voronoi_intersection)
          rm <- which(lengths(surrounding_voronoi) == 0)
          if(length(rm) > 0){
            surrounding_voronoi <- c(1:length(surrounding_voronoi))[-c(rm)]
          } else {
            surrounding_voronoi <- c(1:length(surrounding_voronoi))
          }
          surrounding_voronoi <- closest_voronoi[surrounding_voronoi, ]
          #-------------------------------------------------------------------------------
          
          
          #-------------------------------------------------------------------------------
          # re-creating original points
          closest_stream_points <- stream_points_geometry[closest_points_subset, ]
          closest_stream_points$key <- paste(st_coordinates(closest_stream_points)[ ,1],
                                             st_coordinates(closest_stream_points)[ ,2],
                                             sep = '_')
          closest_stream_points <- closest_stream_points[ ,'key']
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # in case more than one point got collapsed to a single voronoi polygon
          depletions_accounting_duplicates <- lapply(apportioned_depletions$key, function(x){
            length(which(closest_stream_points$key == x))
          })
          depletions_accounting_duplicates <- unlist(depletions_accounting_duplicates)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # dividing depletions among duplicates (if exists)
          apportioned_depletions$dep <- 
            apportioned_depletions$dep/depletions_accounting_duplicates
          apportioned_depletions$divisor <- depletions_accounting_duplicates
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # assigning depletions to points in the correct positions
          frac <- lapply(closest_stream_points$key, function(x){
            apportioned_depletions$dep[apportioned_depletions$key == x]
          })
          frac[lengths(frac) == 0] <- NA
          frac <- unlist(frac)
          closest_stream_points$frac <- frac
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # finding what reaches these depletions apply to
          reach_names <- as.numeric(row.names(closest_stream_points)[(is.na(closest_stream_points$frac) == FALSE)])
          
          reach_names <- st_drop_geometry(stream_points_geometry[reach_names, stream_id_key])
          reach_names <- as.vector(unlist(reach_names))
          closest_stream_points$reach_name <- rep(NA, nrow(closest_stream_points))
          closest_stream_points$reach_name[is.na(closest_stream_points$frac) == FALSE] <- reach_names
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # output
          fractions_of_depletions[i,
                                  as.numeric(colnames(fractions_of_depletions)) %in% as.numeric(reach_names)] <-
            frac[is.na(frac) == FALSE]
          reaches[i, as.numeric(colnames(reaches)) %in% as.numeric(reach_names)] <-
            as.vector(reach_names)
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        # if there are not enough entities to make polygons out of, all depletions go to the
        # only reach included
        if(length(closest_points_subset) == 1){
          #-------------------------------------------------------------------------------
          # finding what reaches these depletions apply to
          reach_names <- as.numeric(row.names(closest_stream_points))
          
          reach_names <- st_drop_geometry(stream_points_geometry[reach_names, stream_id_key])
          reach_names <- as.vector(unlist(reach_names))
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # output
          fractions_of_depletions[i,
                                  as.numeric(colnames(fractions_of_depletions)) %in% as.numeric(reach_names)] <- 1
          reaches[i, as.numeric(colnames(reaches)) %in% as.numeric(reach_names)] <-
            as.vector(reach_names)
          #-------------------------------------------------------------------------------
          
        }
        #-------------------------------------------------------------------------------
      } else {
        # pass
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # correcting for matrix full of NA values
    fractions_of_depletions[is.na(fractions_of_depletions) == TRUE] <- 0
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # writeout statistics
    max_dep <- max(fractions_of_depletions, na.rm = T)
    wm <- which(fractions_of_depletions == max_dep, arr.ind = TRUE)
    
    reach_max_dep <- reaches[wm]
    well_max_dep <- as.vector(unlist(st_drop_geometry(wells[wm[,1],wells_id_key])))
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # format writeout pt 2
    w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
    fractions_of_depletions <- cbind(w_index, fractions_of_depletions)
    fractions_of_depletions <- as.data.frame(fractions_of_depletions)
    colnames(fractions_of_depletions) <- c('wellN',
                                           paste0('RN',
                                                  1:(ncol(fractions_of_depletions)-1)))
    reaches <- cbind(w_index, reaches)
    reaches <- as.data.frame(reaches)
    colnames(reaches) <- c('wellN',
                           paste0('RN',
                                  1:(ncol(reaches)-1)))
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s %s',
                              'Max apportioned depletion fraction: ',
                              paste(round(max_dep,4),
                                    'for reach',
                                    paste(reach_max_dep, collapse = ','),
                                    'for well',
                                    paste(well_max_dep, collapse = ','))),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    return(list(fractions_of_depletions,
                reaches))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  #===========================================================================================
  # Apportion depletions by web or web squared method
  #===========================================================================================
  web_apportionment <- function(power,
                                wells,
                                closest_points_per_segment,
                                stream_points_geometry)
  {
    #-------------------------------------------------------------------------------
    # web equation in zipper (2018)
    # https://doi.org/10.1029/2018WR022707
    equation <- function(well,
                         power = power,
                         closest_points_per_well,
                         stream_points_geometry = stream_points_geometry)
    {
      #-------------------------------------------------------------------------------
      # finding the closest points and the reaches that those points represent
      closest_points <- stream_points_geometry[closest_points_per_well, ]
      reaches <- st_drop_geometry(stream_points_geometry[closest_points_per_well,stream_id_key])
      reaches <- as.vector(unlist(reaches))
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # get all the points per reach
      n_points_per_reach_proximity <- split(stream_points_geometry,
                                            st_drop_geometry(stream_points_geometry[ ,stream_id_key]))
      n_points_per_reach_proximity <- n_points_per_reach_proximity[reaches]
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # get sum per reach of the point distances to the well
      dists <- lapply(n_points_per_reach_proximity, function(x){
        s1 <- st_distance(well,
                          x)
        sum(as.vector(unlist(s1)))
      })
      dists <- as.vector(unlist(dists))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # repeat above but with all points
      n_points_per_reach_all <- split(stream_points_geometry,
                                      st_drop_geometry(stream_points_geometry[ ,stream_id_key]))
      n_points_per_reach_all <- n_points_per_reach_all[reaches]
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # repeat above but with all points
      dists_all <- lapply(n_points_per_reach_all, function(x){
        s1 <- st_distance(well,
                          x)
        sum(as.vector(unlist(s1)))
      })
      dists_all <- as.vector(unlist(dists_all))
      #-------------------------------------------------------------------------------
      
      
      
      
      #-------------------------------------------------------------------------------
      if(geologic_apportionment == FALSE){
        #-------------------------------------------------------------------------------
        numerator <- 1/(dists**power)
        denominator <- sum(1/(dists_all**power))
        fractions_of_depletions <- numerator/denominator
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      
      
      #-------------------------------------------------------------------------------
      if(geologic_apportionment == TRUE){
        #-------------------------------------------------------------------------------
        if(is.null(model_grid) == TRUE){
          #-------------------------------------------------------------------------------
          TR <- as.vector(unlist(st_drop_geometry(well[,well_transmissivity_key])))
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          TR_prime <- lapply(n_points_per_reach_all, function(x){
            
            mean(as.vector(unlist(st_drop_geometry(x[,stream_transmissivity_key]))))
          })
          TR_prime <- as.vector(unlist(TR_prime))
          #-------------------------------------------------------------------------------
          
          
          #-------------------------------------------------------------------------------
          numerator <- (1/(dists**power))*(TR_prime/TR)
          denominator <- sum((1/(dists_all**power))*(TR_prime/TR))
          fractions_of_depletions <- numerator/denominator
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        if(is.null(model_grid) == FALSE){

          #-------------------------------------------------------------------------------
          # get the coordinates of all points on the reach
          reach_coords <- lapply(1:length(n_points_per_reach_proximity),function(i){
            st_coordinates(n_points_per_reach_proximity[[i]])
          })
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # create linestrings from the well to each of these points
          reach_lines_from_well <- lapply(reach_coords, function(x){
            #-------------------------------------------------------------------------------
            # split the coordinates of reach[i] into a list
            # so split_coords[[1]] will be all the points of reach 1
            # with each x,y coordinate in its own list index for the below lapply
            split_coords <- split(x,
                                  1:nrow(x))
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # on each coordinate create a matrix of the well and the point
            # then a linestring out of the point pair
            # at the end, since this is done per reach, rbind the result and it will be a bunch of
            # lines from the well to each point on the reach
            line <- lapply(split_coords,function(x){
              m <- matrix(c(st_coordinates(well)[,1],
                            st_coordinates(well)[,2],
                            x[1],
                            x[2]),
                          ncol = 2,
                          byrow = TRUE)
              line <- st_sf(st_sfc(st_linestring(m), crs = st_crs(wells)))
              st_geometry(line) <- 'geometry'
              line
            })
            do.call(rbind, line)
            #-------------------------------------------------------------------------------
          })
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # selecting correct gridcells
          grid_layers <- as.vector(unlist(st_drop_geometry(model_grid[,grid_layer_key])))
          well_layers <- as.vector(unlist(st_drop_geometry(well[,well_layer_key])))
          gr <- model_grid[grid_layers %in% well_layers, ]
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # taking gridcells along lines connecting reaches to the well
          gridcells_along_lines <- lapply(reach_lines_from_well, function(x){
            step1 <- st_intersects(gr, x)
            inds <- c(1:length(step1))
            rm <- inds[lengths(step1) == 0]
            if(length(rm) > 0){
              inds <- inds[-c(rm)]
            } else {}
            if(length(inds) > 0){
              step1 <- gr[inds,]
            } else {
              step1 <- gr[1,]
              step1[,grid_transmissivity_key] <- NA
            }
          })
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          TR <- lapply(gridcells_along_lines, function(x){
            mean(as.vector(unlist(st_drop_geometry(x[,grid_transmissivity_key]))),
                 na.rm = T)
          })
          TR <- as.vector(unlist(TR))
          TR[is.nan(TR) == TRUE] <- 0
          TR[is.na(TR) == TRUE] <- 0
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          numerator <- (1/(dists**power))*(TR/max(TR))
          denominator <- sum((1/(dists_all**power))*(TR/max(TR)))
          fractions_of_depletions <- numerator/denominator
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      return(list(fractions_of_depletions,
                  reaches))
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # use inverse distance weighting to assign fractions of depletions
    fractions_of_depletions <- as.data.frame(base::matrix(nrow = nrow(wells),
                                                          ncol = nrow(streams)))
    colnames(fractions_of_depletions) <- c(1:ncol(fractions_of_depletions))
    reaches <- as.data.frame(base::matrix(nrow = nrow(wells),
                                          ncol = nrow(streams)))
    colnames(reaches) <- c(1:ncol(reaches))
    for(i in 1:nrow(wells)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = nrow(wells),
                    width = 50,
                    optional_text = 'Apportioning depletions')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      closest_points_per_well <- closest_points_per_segment[i,-c(1)]
      reference <- as.vector(unlist(closest_points_per_well))
      closest_points_per_well <- as.vector(unlist(closest_points_per_well))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # are there any depletions to apportion? if yes continue, if not append NA
      if(all(is.na(closest_points_per_well)) == FALSE){
        
        
        #-------------------------------------------------------------------------------
        # are there any NAs to remove to avoid indexing the stream points by NA
        rm_indices <- which(is.na(closest_points_per_well) == TRUE)
        if(length(rm_indices) > 0){
          #-------------------------------------------------------------------------------
          closest_points_per_well <- closest_points_per_well[-c(rm_indices)] %>%
            as.vector() %>% unlist()
          #-------------------------------------------------------------------------------
        } else {}
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        # assign fractions
        output <- equation(well = wells[i, ],
                           power = power,
                           closest_points_per_well = closest_points_per_well,
                           stream_points_geometry = stream_points_geometry)
        fractions_of_depletions[i,
                                is.na(reference) == FALSE] <- output[[1]]
        reaches[i,
                is.na(reference) == FALSE] <- output[[2]]
        #-------------------------------------------------------------------------------
        
      } else {
        # pass
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # writeout statistics
    max_dep <- max(fractions_of_depletions, na.rm = T)
    wm <- which(fractions_of_depletions == max_dep, arr.ind = TRUE)
    
    reach_max_dep <- reaches[wm]
    well_max_dep <- as.vector(unlist(st_drop_geometry(wells[wm[,1],wells_id_key])))
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # format writeout pt 2
    w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
    fractions_of_depletions <- cbind(w_index, fractions_of_depletions)
    fractions_of_depletions <- as.data.frame(fractions_of_depletions)
    colnames(fractions_of_depletions) <- c('wellN',
                                           paste0('RN',
                                                  1:(ncol(fractions_of_depletions)-1)))
    reaches <- cbind(w_index, reaches)
    reaches <- as.data.frame(reaches)
    colnames(reaches) <- c('wellN',
                           paste0('RN',
                                  1:(ncol(reaches)-1)))
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s %s',
                              'Max apportioned depletion fraction: ',
                              paste(round(max_dep,4),
                                    'for reach',
                                    paste(reach_max_dep, collapse = ','),
                                    'for well',
                                    paste(well_max_dep, collapse = ','))),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    return(list(fractions_of_depletions,
                reaches))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  #===========================================================================================
  # Uses glover model to calculate the depletions per reach
  #===========================================================================================
  glover_stream_depletion_calculations <- function(closest_points_per_segment = closest_points_per_segment,
                                                   stream_points_geometry = stream_points_geometry,
                                                   reach_impact_frac = reach_impact_frac,
                                                   wells = wells,
                                                   well_transmissivity_key = well_transmissivity_key,
                                                   well_stor_coef_key = well_stor_coef_key)
  {
    #===========================================================================================
    # Calculates stream depletions assuming a fully penetrating stream with no
    # clogging layer according to Glover and Balmer (1954) https://doi.org/10.1029/TR035i003p00468
    #===========================================================================================
    glover_stream_depletion_model <- function(stor_coef,
                                              transmissivity,
                                              distance,
                                              QW)
    {
      #-------------------------------------------------------------------------------
      # glover model equation in zipper (2019)
      # https://doi.org/10.1029/2018WR024403
      equation <- function(stor_coef,
                           transmissivity,
                           elapsed_time,
                           distance)
      {
        QA <- erfc(sqrt((stor_coef * distance**2)/
                          (4*transmissivity*elapsed_time)))
        
        return(QA)
      }
      #-------------------------------------------------------------------------------
      
      
      # EXPLANATION
      # The following matrices are an abstraction of the principle of linear superposition.
      # In these matrices, each column is a different pumping rate.
      # This method is necessary as analytical stream depletion functions do not return
      # the depletion at timestep t, but rather the cumulative depletion between 0 and t.
      # Therefore the depletion at timestep t is actually f(t) - f(t-1).
      
      # The matrix [timestep_mat] shows, as stated, each column as a pumping rate and
      # each row as the timesteps. This is the platonic ideal of if all pumping rates
      # started at timestep 1. In this case to get the cumulative depletion at step 1
      # we would just need to for each pumping rate evaluate f(1)*pump and sum them.
      
      # The matrices [starts_mat and stops_mat] represent for each pumping rate (column)
      # when they start and stop. For example column 1 starts has each row set to 0 (starts)
      # at time 0 in [start_mat]. These are less physical representations and more structures
      # that allow us to assemble a physical representation.
      
      # The same is true of [pumping_mat], each column is filled with its representative pumping rate
      # even if it is not active for that timestep
      
      # The matrix [starts_actual] assembles when each pumping rate actually starts,
      # and for how long it has been active. For example columns 1 and 2 may look like
      # 0 0
      # 1 0
      # 2 1
      # 3 2
      # ...
      # showing that at row 4 pumping rate 1 has been active for 3 timesteps, and pumping
      # rate 2 has been active for 2 timestep.
      
      # The matrix [stops_actual] assembles how much time we need to subtract from [starts_actual]
      # to get the impulse at that timestep only. For example columns 1 and 2 may look like
      # 0 0
      # 0 0
      # 1 0
      # 2 1
      # ...
      # so to get the depletion in timestep 4 for column 1 we can use [starts_actual and stops_actual] to evaluate
      # f(3) - f(2). Then for column 2 at timestep 4 we can evaluate f(2) - f(1). The sum of depletions at timestep
      # 4 will then be the sum of these evaluations.
      
      
      # FOR AN EXAMPLE RUN:
      # THIS WILL BE THE SAME AS A CONTINUOUS PUMPING RATE
      # timesteps = c(0,1,2,3,4)
      # start_pumping = c(0,1,2,3)
      # stop_pumping = c(1,2,3,4)
      # QW = c(10,10,10,10)
      #-------------------------------------------------------------------------------
      start_pumping <- c(0:(length(QW)-1))
      stop_pumping <- c(1:length(QW))
      timesteps <- c(0:length(QW)) 
      
      
      timestep_mat <- base::matrix(timesteps,
                                   nrow = length(timesteps),
                                   ncol = length(start_pumping))
      starts_mat <- base::matrix(start_pumping,
                                 nrow = length(timesteps),
                                 ncol = length(start_pumping),
                                 byrow = T)
      stops_mat <- base::matrix(stop_pumping,
                                nrow = length(timesteps),
                                ncol = length(stop_pumping),
                                byrow = T)
      pumping_mat <- base::matrix(QW,
                                  nrow = length(timesteps),
                                  ncol = length(QW),
                                  byrow = T)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # calculate time since each pumping interval starts/stops, bounded at 0
      starts_actual <- timestep_mat - starts_mat
      starts_actual[starts_actual < 0] <- 0
      
      stops_actual <- timestep_mat - stops_mat
      stops_actual[stops_actual < 0] <- 0
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # vectorize for calculations
      starts_actual_vec <- c(starts_actual)
      stops_actual_vec <- c(stops_actual)
      pumping_vec <- c(pumping_mat)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # vector of zeroes that will be filled with the function evaluations
      depletions_vec <- rep(0, length(starts_actual_vec))
      fractional_vec <- rep(0, length(starts_actual_vec))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # evaulate f(t) - f(t-1)
      depletions_vec[starts_actual_vec > 0] <-
        pumping_vec[starts_actual_vec > 0] *
        (equation(elapsed_time = starts_actual_vec[starts_actual_vec > 0],
                  distance = distance,
                  stor_coef = stor_coef,
                  transmissivity = transmissivity) -
           equation(elapsed_time = stops_actual_vec[starts_actual_vec > 0],
                    distance = distance,
                    stor_coef = stor_coef,
                    transmissivity = transmissivity))
      
      depletions_mat <- matrix(depletions_vec,
                               nrow = length(timesteps),
                               ncol = length(start_pumping))
      depletions_mat <- depletions_mat[-c(1), ]
      depletions <- base::rowSums(depletions_mat)
        
        
        
        
        

      fractional_vec[starts_actual_vec > 0] <-
        (equation(elapsed_time = starts_actual_vec[starts_actual_vec > 0],
                  distance = distance,
                  stor_coef = stor_coef,
                  transmissivity = transmissivity) -
           equation(elapsed_time = stops_actual_vec[starts_actual_vec > 0],
                    distance = distance,
                    stor_coef = stor_coef,
                    transmissivity = transmissivity))
        
      fractional_mat <- matrix(fractional_vec,
                               nrow = length(timesteps),
                               ncol = length(start_pumping))
      fractional_mat <- fractional_mat[-c(1), ]
      fractions <- base::rowSums(fractional_mat)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # memory allocation
      rm(list = c('depletions_vec','depletions_mat',
                  'fractional_mat','fractional_vec',
                  'starts_actual','stops_actual',
                  'timestep_mat','starts_mat',
                  'stops_mat','pumping_mat'))
      #-------------------------------------------------------------------------------


      #-------------------------------------------------------------------------------
      # sum and return
      Jenk_SDF <- (distance*distance*stor_coef)/transmissivity
      return(list(depletions,
                  fractions,
                  Jenk_SDF))
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    #-------------------------------------------------------------------------------
    # for each reach calculate sum of all depletions
    depletions_per_reach <- list()
    depletions_potential_per_reach <- list()
    pump_frac_per_reach <- list()
    jenk_sdf_per_reach <- list()
    custom_sdf_per_reach <- list()
    for(i in 1:ncol(closest_points_per_segment)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = ncol(closest_points_per_segment),
                    width = 50,
                    optional_text = '')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # what are the closest points to each well
      points <- as.vector(unlist(closest_points_per_segment[, i]))
      fracs <- as.vector(unlist(reach_impact_frac[ , i]))
      well_indices <- c(1:length(points))
      rm <- which(is.na(points))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if no wells are associated with the reach then pass
      if(all(is.na(points) == TRUE) == FALSE){
        #-------------------------------------------------------------------------------
        # remove any non-relevant wells from the for loop
        if(length(rm) > 0){
          well_indices <- well_indices[-c(rm)]
        } else {}
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        depletions_per_well <- list()
        depletions_potential_per_well <- list()
        pump_frac_per_well <- list()
        distances <- list()
        Jenk_SDF_per_well <- list()
        custom_SDF_per_well <- list()
        counter <- 0
        for(j in well_indices){
          #-------------------------------------------------------------------------------
          # increment list counter
          counter <- counter + 1
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get distances
          distance <- st_distance(wells[j, ],
                                  stream_points_geometry[points[j], ])
          distance <- as.numeric(distance)
          distances[[counter]] <- distance
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get aquifer properties depending on if its assigned on a per-well basis
          # or if its going to be taken as an average of the aquifer properties along a line
          # between the well and the nearest stream
          if(is.null(model_grid) == TRUE){
            transmissivity <- as.numeric(st_drop_geometry(wells[j,well_transmissivity_key]))
            stor_coef <- as.numeric(st_drop_geometry(wells[j, well_stor_coef_key]))
          } else {
            #-------------------------------------------------------------------------------
            # make line between well and stream
            m <- matrix(c(st_coordinates(wells[j, ])[,1],
                          st_coordinates(wells[j, ])[,2],
                          st_coordinates(stream_points_geometry[points[j], ])[,1],
                          st_coordinates(stream_points_geometry[points[j], ])[,2]),
                        ncol = 2,
                        byrow = TRUE)
            line <- st_sf(st_sfc(st_linestring(m), crs = st_crs(wells)))
            st_geometry(line) <- 'geometry'
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # selecting correct gridcells
            grid_layers <- as.vector(unlist(st_drop_geometry(model_grid[,grid_layer_key])))
            well_layers <- as.vector(unlist(st_drop_geometry(wells[,well_layer_key])))
            gr <- model_grid[grid_layers == well_layers[j], ]
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # intersect line between well and stream with the grid
            int <- st_intersects(gr, line$geometry)
            grid_inds <- c(1:length(int))
            rm <- which(lengths(int) == 0)
            if(length(rm) > 0){
              grid_inds <- grid_inds[-c(rm)]
            } else {}
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # getting average properties from grid
            gr <- gr[grid_inds, ]
            transmissivity <- mean(as.vector(unlist(st_drop_geometry(gr[,grid_transmissivity_key]))), na.rm = T)
            stor_coef <- mean(as.vector(unlist(st_drop_geometry(gr[,grid_stor_coef_key]))), na.rm = T)
            #-------------------------------------------------------------------------------
          }
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # calculate maximum stream depletion potential and
          # multiply by fraction of depletions of this well apportioned to this reach
          Q_out <- glover_stream_depletion_model(stor_coef = stor_coef,
                                                 transmissivity = transmissivity,
                                                 distance = distance,
                                                 QW = pumping[j, ])
          Q_final <- Q_out[[1]]*fracs[j]
          Q_fraction <- Q_out[[2]]
          Jenk_SDF <- Q_out[[3]]
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          depletions_per_well[[counter]] <- Q_final
          depletions_potential_per_well[[counter]] <- Q_fraction
          Jenk_SDF_per_well[[counter]] <- Jenk_SDF
          pump_frac_per_well[[counter]] <- pumping[j, ] * fracs[j]
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        depletions_total <- do.call(cbind, depletions_per_well)
        pump_frac_per_well_total <- do.call(cbind, pump_frac_per_well)
        depletions_potential_per_well_total <- do.call(cbind, depletions_potential_per_well)
        Jenk_SDF_per_well_total <- do.call(rbind, Jenk_SDF_per_well)
        #-------------------------------------------------------------------------------
        
        average_fractional_depletions <- calculate_depletion_potential(depletion_potential_criteria = depletion_potential_criteria,
                                                                       depletions_potential_per_well_total = depletions_potential_per_well_total,
                                                                       distances = distances,
                                                                       fracs = fracs,
                                                                       pumping = pumping)
        average_Jenk_SDF <- calculate_average_sdf(sdf_averaging_criteria = sdf_averaging_criteria,
                                                  sdf_vec = Jenk_SDF_per_well_total,
                                                  fracs = fracs,
                                                  pumping = pumping)
        
        jenk_sdf_per_reach[[i]] <- average_Jenk_SDF
        depletions_potential_per_reach[[i]] <- average_fractional_depletions
        depletions_per_reach[[i]] <- base::rowSums(depletions_total)
        pump_frac_per_reach[[i]] <- base::rowSums(pump_frac_per_well_total)
        
        
        if(is.null(custom_sdf_time) == FALSE){

          custom_sdf_per_reach[[i]] <- calculate_custom_sdf_time(average_fractional_depletions,
                                                                 target = custom_sdf_time)
        }
      } else{
        depletions_per_reach[[i]] <- rep(0, ncol(pumping)) # reach has no depletions
        pump_frac_per_reach[[i]] <- rep(0, ncol(pumping))
        depletions_potential_per_reach[[i]] <- rep(0, ncol(pumping))
        jenk_sdf_per_reach[[i]] <- NA
        
        
        if(is.null(custom_sdf_time) == FALSE){
          custom_sdf_per_reach[[i]] <- NA
        }
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    

    
    #-------------------------------------------------------------------------------
    # stats
    start_of_depletions <- lapply(depletions_per_reach, function(x){
      rle(x)$lengths[1]
    })
    start_of_depletions <- unlist(start_of_depletions)
    rm <- which(start_of_depletions == ncol(pumping))
    if(length(rm) > 0){
      start_of_depletions[-c(rm)] # if never started remove
    } else {}
    
    
    
    mean_start_of_depletions <- mean(start_of_depletions, na.rm = TRUE)
    median_start_of_depletions <- median(start_of_depletions, na.rm = TRUE)
    
    final_depletions <- lapply(depletions_per_reach, function(x){
      tail(x, 1)
    })
    n_timesteps <- ncol(pumping)
    mean_final_depletions <- mean(unlist(final_depletions), na.rm = TRUE)
    median_final_depletions <- median(unlist(final_depletions), na.rm = TRUE)
    which_max_final_depletions <- which.max(unlist(final_depletions))
    max_final_depletions <- max(unlist(final_depletions), na.rm = TRUE)
    #-------------------------------------------------------------------------------

    
    
    #-------------------------------------------------------------------------------
    # write status to log
    u <- paste0('(',units,'^3','):')

    
    writeLines(text = sprintf('%s %s',
                              'Mean | Median start of stream depletions (timestep): ',
                              paste(round(mean_start_of_depletions,2),
                                    ' | ',
                                    median_start_of_depletions)),
               con = log_file)
    
    
    writeLines(text = sprintf('%s %s',
                              paste('Mean | Median final depletions',u),
                              paste(round(mean_final_depletions,4),
                                    '|',
                                    round(median_final_depletions,4),
                                    'at timestep (t_final)',
                                    n_timesteps)),
               con = log_file)
    
    writeLines(text = sprintf('%s %s',
                              paste('Max final depletions',u),
                              paste(round(max_final_depletions,4),
                                    'at timestep (t_final)',
                                    n_timesteps,
                                    'for reach',
                                    which_max_final_depletions)),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    

    #-------------------------------------------------------------------------------
    # output
    depletions_potential_per_reach <- do.call(rbind, depletions_potential_per_reach)
    depletions_per_reach <- do.call(rbind, depletions_per_reach)
    pump_frac_per_reach <- do.call(rbind, pump_frac_per_reach)
    jenk_sdf_per_reach <- do.call(rbind, jenk_sdf_per_reach)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # output and stats 2
    writeLines(text = sprintf('%s %s',
                              paste('Mean | Median Jenkins SDF'),
                              paste(round(mean(jenk_sdf_per_reach, na.rm = TRUE),4),
                                    '|',
                                    round(median(jenk_sdf_per_reach, na.rm = TRUE),4))),
               con = log_file)
    if(is.null(custom_sdf_time) == TRUE){
      custom_sdf_per_reach <- NULL
    } else{
      custom_sdf_per_reach <- do.call(rbind, custom_sdf_per_reach)
      writeLines(text = sprintf('%s %s',
                                paste('Mean | Median Custom SDF',paste0('(',custom_sdf_time,')')),
                                paste(round(mean(custom_sdf_per_reach, na.rm = TRUE),4),
                                      '|',
                                      round(median(custom_sdf_per_reach, na.rm = TRUE),4))),
                 con = log_file)
    }
    #-------------------------------------------------------------------------------
    return(list(depletions_per_reach,
                depletions_potential_per_reach,
                pump_frac_per_reach,
                jenk_sdf_per_reach,
                custom_sdf_per_reach))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  
  #===========================================================================================
  # Uses hunt model to calculate the depletions per reach
  #===========================================================================================
  hunt_stream_depletion_calculations <- function(closest_points_per_segment = closest_points_per_segment,
                                                 stream_points_geometry = stream_points_geometry,
                                                 reach_impact_frac = reach_impact_frac,
                                                 wells = wells,
                                                 well_transmissivity_key = well_transmissivity_key,
                                                 well_stor_coef_key = well_stor_coef_key,
                                                 lambda_key = lambda_key)
  {
    #===========================================================================================
    # Calculates stream depletions assuming a fully penetrating stream with no
    # clogging layer according to Hunt (1999) https://doi.org/10.1111/j.1745-6584.1999.tb00962.x
    # identical to Hantush (1965) model in special case that lambda = 2*(T/L)
    #===========================================================================================
    hunt_stream_depletion_model <- function(stor_coef,
                                            transmissivity,
                                            distance,
                                            lambda,
                                            QW)
    {
      #-------------------------------------------------------------------------------
      # hunt model equation in zipper (2019)
      # https://doi.org/10.1029/2018WR024403
      # for wide rivers, low aquifer transmissivity, and long simulations
      # exponential terms can easily be greater than exp(1e5) so uses
      # Rmpfr approach of streamdepletR to store higher precision numbers
      # this takes more memory and time so this is only done if infinite
      # numbers are produced to speed processing
      equation <- function(stor_coef,
                           transmissivity,
                           elapsed_time,
                           distance,
                           lambda)
      {
        #-------------------------------------------------------------------------------
        t2_a <- ((lambda*lambda*elapsed_time)/(4*stor_coef*transmissivity))
        t2_b <- ((lambda*distance)/(2*transmissivity))
        t2 <- base::exp(t2_a + t2_b)
        
        infinite_indices <- which(is.infinite(t2))
        all_indices <- c(1:length(elapsed_time))
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        # check whether exponentials contain any infinities
        if(length(infinite_indices) == 0){
          #-------------------------------------------------------------------------------
          # assemble terms
          z <- (sqrt((stor_coef * distance* distance)/
                       (4*transmissivity*elapsed_time)))
          t1 <- erfc(z)
          
          t3_a <- (sqrt((lambda*lambda*elapsed_time)/(4*stor_coef*transmissivity)))
          t3 <- erfc(t3_a + z)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # give final answer
          QA <- as.numeric(t1 - (t2*t3))
          #-------------------------------------------------------------------------------
        } else {
          
          #-------------------------------------------------------------------------------
          # blank fill
          QA <- rep(NA, length(elapsed_time))
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # assemble terms
          z <- Rmpfr::mpfr((sqrt((stor_coef * distance* distance)/
                                   (4*transmissivity*elapsed_time[infinite_indices]))), prec = prec)
          t1 <- Rmpfr::erfc(z)
          
          
          t2_a <- Rmpfr::mpfr(((lambda*lambda*elapsed_time[infinite_indices])/(4*stor_coef*transmissivity)),  prec = prec)
          t2_b <- Rmpfr::mpfr(((lambda*distance)/(2*transmissivity)), prec = prec)
          t2 <- base::exp(t2_a + t2_b)
          
          
          t3_a <- Rmpfr::mpfr((sqrt((lambda*lambda*elapsed_time[infinite_indices])/(4*stor_coef*transmissivity))),  prec = prec)
          t3 <- Rmpfr::erfc(t3_a + z)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # fill infinite indices with higher precision numbers
          QA[infinite_indices] <- as.numeric(t1 - (t2*t3))
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # are there non-infinity generating terms in base packages that should be calculated
          if(length(infinite_indices) != length(elapsed_time)){
            #-------------------------------------------------------------------------------
            # assemble terms
            z <- (sqrt((stor_coef * distance* distance)/
                         (4*transmissivity*elapsed_time[all_indices[-c(infinite_indices)]])))
            t1 <- erfc(z)
            
            t2_a <- ((lambda*lambda*elapsed_time[all_indices[-c(infinite_indices)]])/(4*stor_coef*transmissivity))
            t2_b <- ((lambda*distance)/(2*transmissivity))
            t2 <- base::exp(t2_a + t2_b)
            
            
            
            t3_a <- (sqrt((lambda*lambda*elapsed_time[all_indices[-c(infinite_indices)]])/(4*stor_coef*transmissivity)))
            t3 <- erfc(t3_a + z)
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # give final answer
            QA[all_indices[-c(infinite_indices)]] <- as.numeric(t1 - (t2*t3))
            #-------------------------------------------------------------------------------
            
          }
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        rm(list = c('t1','t2','t2_a','t2_b','t3','t3_a','z'))
        return(QA)
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      
      # EXPLANATION
      # check comments on glover model
      #-------------------------------------------------------------------------------
      start_pumping <- c(0:(length(QW)-1))
      stop_pumping <- c(1:length(QW))
      timesteps <- c(0:length(QW)) 
      
      
      timestep_mat <- base::matrix(timesteps,
                                   nrow = length(timesteps),
                                   ncol = length(start_pumping))
      starts_mat <- base::matrix(start_pumping,
                                 nrow = length(timesteps),
                                 ncol = length(start_pumping),
                                 byrow = T)
      stops_mat <- base::matrix(stop_pumping,
                                nrow = length(timesteps),
                                ncol = length(stop_pumping),
                                byrow = T)
      pumping_mat <- base::matrix(QW,
                                  nrow = length(timesteps),
                                  ncol = length(QW),
                                  byrow = T)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # calculate time since each pumping interval starts/stops, bounded at 0
      starts_actual <- timestep_mat - starts_mat
      starts_actual[starts_actual < 0] <- 0
      
      stops_actual <- timestep_mat - stops_mat
      stops_actual[stops_actual < 0] <- 0
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # vectorize for calculations
      starts_actual_vec <- c(starts_actual)
      stops_actual_vec <- c(stops_actual)
      pumping_vec <- c(pumping_mat)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # vector of zeroes that will be filled with the function evaluations
      depletions_vec <- rep(0, length(starts_actual_vec))
      fractional_vec <- rep(0, length(starts_actual_vec))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # evaulate f(t) - f(t-1)
      depletions_vec[starts_actual_vec > 0] <-
        pumping_vec[starts_actual_vec > 0] *
        (equation(elapsed_time = starts_actual_vec[starts_actual_vec > 0],
                  distance = distance,
                  stor_coef = stor_coef,
                  transmissivity = transmissivity,
                  lambda = lambda) -
           equation(elapsed_time = stops_actual_vec[starts_actual_vec > 0],
                    distance = distance,
                    stor_coef = stor_coef,
                    transmissivity = transmissivity,
                    lambda = lambda))
      
      depletions_mat <- matrix(depletions_vec,
                               nrow = length(timesteps),
                               ncol = length(start_pumping))
      depletions_mat <- depletions_mat[-c(1), ]
      depletions <- base::rowSums(depletions_mat)
      
      
      
      
      
      
      fractional_vec[starts_actual_vec > 0] <-
        (equation(elapsed_time = starts_actual_vec[starts_actual_vec > 0],
                  distance = distance,
                  stor_coef = stor_coef,
                  transmissivity = transmissivity,
                  lambda = lambda) -
           equation(elapsed_time = stops_actual_vec[starts_actual_vec > 0],
                    distance = distance,
                    stor_coef = stor_coef,
                    transmissivity = transmissivity,
                    lambda = lambda))
      
      fractional_mat <- matrix(fractional_vec,
                               nrow = length(timesteps),
                               ncol = length(start_pumping))
      fractional_mat <- fractional_mat[-c(1), ]
      fractions <- base::rowSums(fractional_mat)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # memory allocation
      rm(list = c('depletions_vec','depletions_mat',
                  'fractional_mat','fractional_vec',
                  'starts_actual','stops_actual',
                  'timestep_mat','starts_mat',
                  'stops_mat','pumping_mat'))
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # sum and return
      Jenk_SDF <- (distance*distance*stor_coef)/transmissivity
      return(list(depletions,
                  fractions,
                  Jenk_SDF))
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    #-------------------------------------------------------------------------------
    # for each reach calculate sum of all depletions
    depletions_per_reach <- list()
    pump_frac_per_reach <- list()
    depletions_potential_per_reach <- list()
    jenk_sdf_per_reach <- list()
    custom_sdf_per_reach <- list()
    for(i in 1:ncol(closest_points_per_segment)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = ncol(closest_points_per_segment),
                    width = 50,
                    optional_text = '')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # what are the closest points to each well
      points <- as.vector(unlist(closest_points_per_segment[, i]))
      fracs <- as.vector(unlist(reach_impact_frac[ , i]))
      RN <- st_drop_geometry(stream_points_geometry[points, stream_id_key])
      RN <- as.vector(unlist(RN))
      well_indices <- c(1:length(points))
      rm <- which(is.na(points))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if no wells are associated with the reach then pass
      if(all(is.na(points) == TRUE) == FALSE){
        #-------------------------------------------------------------------------------
        # remove any non-relevant wells from the for loop
        if(length(rm) > 0){
          well_indices <- well_indices[-c(rm)]
        } else {}
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        depletions_per_well <- list()
        pump_frac_per_well <- list()
        depletions_potential_per_well <- list()
        Jenk_SDF_per_well <- list()
        custom_SDF_per_well <- list()
        counter <- 0
        for(j in well_indices){
          #-------------------------------------------------------------------------------
          # increment list counter
          counter <- counter + 1
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get distances
          distance <- st_distance(wells[j, ],
                                  stream_points_geometry[points[j], ])
          distance <- as.numeric(distance)
          distances[[counter]] <- distance
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get aquifer properties depending on if its assigned on a per-well basis
          # or if its going to be taken as an average of the aquifer properties along a line
          # between the well and the nearest stream
          if(is.null(model_grid) == TRUE){
            transmissivity <- as.numeric(st_drop_geometry(wells[j,well_transmissivity_key]))
            stor_coef <- as.numeric(st_drop_geometry(wells[j, well_stor_coef_key]))
          } else {
            #-------------------------------------------------------------------------------
            # make line between well and stream
            m <- matrix(c(st_coordinates(wells[j, ])[,1],
                          st_coordinates(wells[j, ])[,2],
                          st_coordinates(stream_points_geometry[points[j], ])[,1],
                          st_coordinates(stream_points_geometry[points[j], ])[,2]),
                        ncol = 2,
                        byrow = TRUE)
            line <- st_sf(st_sfc(st_linestring(m), crs = st_crs(wells)))
            st_geometry(line) <- 'geometry'
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # selecting correct gridcells
            grid_layers <- as.vector(unlist(st_drop_geometry(model_grid[,grid_layer_key])))
            well_layers <- as.vector(unlist(st_drop_geometry(wells[,well_layer_key])))
            gr <- model_grid[grid_layers == well_layers[j], ]
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # intersect line between well and stream with the grid
            int <- st_intersects(gr, line$geometry)
            grid_inds <- c(1:length(int))
            rm <- which(lengths(int) == 0)
            if(length(rm) > 0){
              grid_inds <- grid_inds[-c(rm)]
            } else {}
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # getting average properties from grid
            gr <- gr[grid_inds, ]
            transmissivity <- mean(as.vector(unlist(st_drop_geometry(gr[,grid_transmissivity_key]))), na.rm = T)
            stor_coef <- mean(as.vector(unlist(st_drop_geometry(gr[,grid_stor_coef_key]))), na.rm = T)
            #-------------------------------------------------------------------------------
          }
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # take weighted mean of lambda along the considered reach
          # logic is that closer points will control more
          reaches <- as.vector(unlist(st_drop_geometry(stream_points_geometry[,stream_id_key])))
          stream_inds <- reaches == RN[j]
          all_distances <- st_distance(wells[j, ],
                                       stream_points_geometry[stream_inds, ])
          all_distances <- c(1:length(all_distances))[order(as.numeric(all_distances))]
          
          lambda <- as.vector(unlist(st_drop_geometry(stream_points_geometry[stream_inds, lambda_key])))

          lambda <- weighted_mean(x = lambda, w = all_distances, na.rm = T)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # calculate maximum stream depletion potential and
          # multiply by fraction of depletions of this well apportioned to this reach
          Q_out <- hunt_stream_depletion_model(stor_coef = stor_coef,
                                               transmissivity = transmissivity,
                                               distance = distance,
                                               QW = pumping[j, ],
                                               lambda = lambda)
          Q_final <- Q_out[[1]]*fracs[j]
          Q_fraction <- Q_out[[2]]
          Jenk_SDF <- Q_out[[3]]
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # gradient descent example for future use
          # if(is.null(custom_sdf_time) == FALSE){
          #   custom_SDF <- custom_sdf_gradient_descent(analytical_model = analytical_model,
          #                                             distance = distance,
          #                                             stor_coef = stor_coef,
          #                                             transmissivity = transmissivity,
          #                                             lambda = lambda,
          #                                             custom_sdf_time = custom_sdf_time)
          #   custom_SDF_per_well[[counter]] <- custom_SDF
          # }
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          depletions_per_well[[counter]] <- Q_final
          depletions_potential_per_well[[counter]] <- Q_fraction
          Jenk_SDF_per_well[[counter]] <- Jenk_SDF
          pump_frac_per_well[[counter]] <- pumping[j, ] * fracs[j]
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        depletions_total <- do.call(cbind, depletions_per_well)
        pump_frac_per_well_total <- do.call(cbind, pump_frac_per_well)
        depletions_potential_per_well_total <- do.call(cbind, depletions_potential_per_well)
        Jenk_SDF_per_well_total <- do.call(rbind, Jenk_SDF_per_well)
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        average_fractional_depletions <- calculate_depletion_potential(depletion_potential_criteria = depletion_potential_criteria,
                                                                       depletions_potential_per_well_total = depletions_potential_per_well_total,
                                                                       distances = distances,
                                                                       fracs = fracs,
                                                                       pumping = pumping)
        average_Jenk_SDF <- calculate_average_sdf(sdf_averaging_criteria = sdf_averaging_criteria,
                                                  sdf_vec = Jenk_SDF_per_well_total,
                                                  fracs = fracs,
                                                  pumping = pumping)
        
        jenk_sdf_per_reach[[i]] <- average_Jenk_SDF
        depletions_potential_per_reach[[i]] <- average_fractional_depletions
        depletions_per_reach[[i]] <- base::rowSums(depletions_total)
        pump_frac_per_reach[[i]] <- base::rowSums(pump_frac_per_well_total)
        
        if(is.null(custom_sdf_time) == FALSE){
          
          custom_sdf_per_reach[[i]] <- calculate_custom_sdf_time(average_fractional_depletions,
                                                                 target = custom_sdf_time)
        }
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        # gradient descent example for future use
        # if(is.null(custom_sdf_time) == FALSE){
        #   custom_SDF_per_well_total <- do.call(rbind, custom_SDF_per_well)
        #   custom_SDF_per_well_total[custom_SDF_per_well_total == -9999] <- NA
        #   if(all(is.na(custom_SDF_per_well_total)) == FALSE){
        #     average_custom_sdf <- calculate_average_sdf(sdf_averaging_criteria = sdf_averaging_criteria,
        #                                                 sdf_vec = custom_SDF_per_well_total,
        #                                                 fracs = fracs,
        #                                                 pumping = pumping)
        #     custom_sdf_per_reach[[i]] <- average_custom_sdf
        #   } else {
        #     custom_sdf_per_reach[[i]] <- -9999
        #   }
        # }
        #-------------------------------------------------------------------------------
      } else{
        depletions_per_reach[[i]] <- rep(0, ncol(pumping)) # reach has no depletions
        pump_frac_per_reach[[i]] <- rep(0, ncol(pumping))
        depletions_potential_per_reach[[i]] <- rep(0, ncol(pumping))
        jenk_sdf_per_reach[[i]] <- NA
        
        if(is.null(custom_sdf_time) == FALSE){
          custom_sdf_per_reach[[i]] <- NA
        }
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    
    
    #-------------------------------------------------------------------------------
    # stats
    start_of_depletions <- lapply(depletions_per_reach, function(x){
      rle(x)$lengths[1]
    })
    start_of_depletions <- unlist(start_of_depletions)
    rm <- which(start_of_depletions == ncol(pumping))
    if(length(rm) > 0){
      start_of_depletions[-c(rm)] # if never started remove
    } else {}
    
    
    
    mean_start_of_depletions <- mean(start_of_depletions, na.rm = TRUE)
    median_start_of_depletions <- median(start_of_depletions, na.rm = TRUE)
    
    final_depletions <- lapply(depletions_per_reach, function(x){
      tail(x, 1)
    })
    n_timesteps <- ncol(pumping)
    mean_final_depletions <- mean(unlist(final_depletions), na.rm = TRUE)
    median_final_depletions <- median(unlist(final_depletions), na.rm = TRUE)
    which_max_final_depletions <- which.max(unlist(final_depletions))
    max_final_depletions <- max(unlist(final_depletions), na.rm = TRUE)
    #-------------------------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------------------------
    # write status to log
    u <- paste0('(',units,'^3','):')
    
    
    writeLines(text = sprintf('%s %s',
                              'Mean | Median start of stream depletions (timestep): ',
                              paste(round(mean_start_of_depletions,2),
                                    ' | ',
                                    median_start_of_depletions)),
               con = log_file)
    
    
    writeLines(text = sprintf('%s %s',
                              paste('Mean | Median final depletions',u),
                              paste(round(mean_final_depletions,4),
                                    '|',
                                    round(median_final_depletions,4),
                                    'at timestep (t_final)',
                                    n_timesteps)),
               con = log_file)
    
    writeLines(text = sprintf('%s %s',
                              paste('Max final depletions',u),
                              paste(round(max_final_depletions,4),
                                    'at timestep (t_final)',
                                    n_timesteps,
                                    'for reach',
                                    which_max_final_depletions)),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    
    
    
    
    
    
    #-------------------------------------------------------------------------------
    # output
    depletions_potential_per_reach <- do.call(rbind, depletions_potential_per_reach)
    depletions_per_reach <- do.call(rbind, depletions_per_reach)
    pump_frac_per_reach <- do.call(rbind, pump_frac_per_reach)
    jenk_sdf_per_reach <- do.call(rbind, jenk_sdf_per_reach)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # output and stats 2
    writeLines(text = sprintf('%s %s',
                              paste('Mean | Median Jenkins SDF'),
                              paste(round(mean(jenk_sdf_per_reach, na.rm = TRUE),4),
                                    '|',
                                    round(median(jenk_sdf_per_reach, na.rm = TRUE),4))),
               con = log_file)
    if(is.null(custom_sdf_time) == TRUE){
      custom_sdf_per_reach <- NULL
    } else{
      custom_sdf_per_reach <- do.call(rbind, custom_sdf_per_reach)
      writeLines(text = sprintf('%s %s',
                                paste('Mean | Median Custom SDF',paste0('(',custom_sdf_time,')')),
                                paste(round(mean(custom_sdf_per_reach, na.rm = TRUE),4),
                                      '|',
                                      round(median(custom_sdf_per_reach, na.rm = TRUE),4))),
                 con = log_file)
    }
    #-------------------------------------------------------------------------------
    return(list(depletions_per_reach,
                depletions_potential_per_reach,
                pump_frac_per_reach,
                jenk_sdf_per_reach,
                custom_sdf_per_reach))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  #===========================================================================================
  # Uses hunt model to calculate the depletions per reach
  #===========================================================================================
  hantush_stream_depletion_calculations <- function(closest_points_per_segment = closest_points_per_segment,
                                                    stream_points_geometry = stream_points_geometry,
                                                    reach_impact_frac = reach_impact_frac,
                                                    wells = wells,
                                                    well_transmissivity_key = well_transmissivity_key,
                                                    well_stor_coef_key = well_stor_coef_key,
                                                    leakance_key = leakance_key)
  {
    #===========================================================================================
    # Calculates stream depletions assuming a fully penetrating stream with no
    # clogging layer according to Hantush (1965) https://doi.org/10.1029/JZ070i012p02829
    # identical to Hunt (1999) model in special case that lambda = 2*(T/L)
    # by observation of author for leakance values <= 10 this is just a slower glover solution
    # value of leakance < 1 not acceptable
    #===========================================================================================
    hantush_stream_depletion_model <- function(stor_coef,
                                               transmissivity,
                                               distance,
                                               leakance,
                                               QW)
    {
      #-------------------------------------------------------------------------------
      # hantush model equation in reeves (2008)
      # https://doi.org/10.3133/ofr20081166
      # for wide rivers, low aquifer transmissivity, and long simulations
      # exponential terms can easily be greater than exp(1e5) so uses
      # Rmpfr approach of streamdepletR to store higher precision numbers
      # this takes more memory and time so this is only done if infinite
      # numbers are produced to speed processing
      equation <- function(stor_coef,
                           transmissivity,
                           elapsed_time,
                           distance,
                           leakance)
      {
        
        #-------------------------------------------------------------------------------
        t2_a <- ((transmissivity*elapsed_time)/(stor_coef*leakance*leakance))
        t2_b <- (distance/leakance)
        t2 <- base::exp(t2_a + t2_b)
        
        infinite_indices <- which(is.infinite(t2))
        all_indices <- c(1:length(elapsed_time))
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        # check whether exponentials contain any infinities
        if(length(infinite_indices) == 0){
          #-------------------------------------------------------------------------------
          # assemble terms
          z <- (sqrt((stor_coef * distance* distance)/
                       (4*transmissivity*elapsed_time)))
          t1 <- erfc(z)
          
          t3_a <- (sqrt((transmissivity*elapsed_time)/(stor_coef*leakance*leakance)))
          t3 <- erfc(t3_a + z)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # give final answer
          QA <- as.numeric(t1 - (t2*t3))
          #-------------------------------------------------------------------------------
        } else {
          
          #-------------------------------------------------------------------------------
          # blank fill
          QA <- rep(NA, length(elapsed_time))
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # assemble terms
          z <- Rmpfr::mpfr((sqrt((stor_coef * distance* distance)/
                                   (4*transmissivity*elapsed_time[infinite_indices]))),  prec = prec)
          t1 <- Rmpfr::erfc(z)
          
          
          t2_a <- Rmpfr::mpfr(((transmissivity*elapsed_time[infinite_indices])/(stor_coef*leakance*leakance)),  prec = prec)
          t2_b <- Rmpfr::mpfr((distance/leakance), prec = prec)
          t2 <- base::exp(t2_a + t2_b)
          
          
          t3_a <- Rmpfr::mpfr((sqrt((transmissivity*elapsed_time[infinite_indices])/(stor_coef*leakance*leakance))),  prec = prec)
          t3 <- Rmpfr::erfc(t3_a + z)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # fill infinite indices with higher precision numbers
          QA[infinite_indices] <- as.numeric(t1 - (t2*t3))
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # are there non-infinity generating terms in base packages that should be calculated
          if(length(infinite_indices) != length(elapsed_time)){
            #-------------------------------------------------------------------------------
            # assemble terms
            z <- (sqrt((stor_coef * distance* distance)/
                         (4*transmissivity*elapsed_time[all_indices[-c(infinite_indices)]])))
            t1 <- erfc(z)
            
            t2_a <- ((transmissivity*elapsed_time[all_indices[-c(infinite_indices)]])/(stor_coef*leakance*leakance))
            t2_b <- (distance/leakance)
            t2 <- base::exp(t2_a + t2_b)
            
            t3_a <- (sqrt((transmissivity*elapsed_time[all_indices[-c(infinite_indices)]])/(stor_coef*leakance*leakance)))
            t3 <- erfc(t3_a + z)
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # give final answer
            QA[all_indices[-c(infinite_indices)]] <- as.numeric(t1 - (t2*t3))
            #-------------------------------------------------------------------------------
            
          }
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------
        rm(list = c('t1','t2','t2_a','t2_b','t3','t3_a','z'))
        return(QA)
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      
      # EXPLANATION
      # check comments on glover model
      #-------------------------------------------------------------------------------
      start_pumping <- c(0:(length(QW)-1))
      stop_pumping <- c(1:length(QW))
      timesteps <- c(0:length(QW)) 
      
      
      timestep_mat <- base::matrix(timesteps,
                                   nrow = length(timesteps),
                                   ncol = length(start_pumping))
      starts_mat <- base::matrix(start_pumping,
                                 nrow = length(timesteps),
                                 ncol = length(start_pumping),
                                 byrow = T)
      stops_mat <- base::matrix(stop_pumping,
                                nrow = length(timesteps),
                                ncol = length(stop_pumping),
                                byrow = T)
      pumping_mat <- base::matrix(QW,
                                  nrow = length(timesteps),
                                  ncol = length(QW),
                                  byrow = T)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # calculate time since each pumping interval starts/stops, bounded at 0
      starts_actual <- timestep_mat - starts_mat
      starts_actual[starts_actual < 0] <- 0
      
      stops_actual <- timestep_mat - stops_mat
      stops_actual[stops_actual < 0] <- 0
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # vectorize for calculations
      starts_actual_vec <- c(starts_actual)
      stops_actual_vec <- c(stops_actual)
      pumping_vec <- c(pumping_mat)
      #-------------------------------------------------------------------------------
      
      
      
      #-------------------------------------------------------------------------------
      # vector of zeroes that will be filled with the function evaluations
      depletions_vec <- rep(0, length(starts_actual_vec))
      fractional_vec <- rep(0, length(starts_actual_vec))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # evaulate f(t) - f(t-1)
      depletions_vec[starts_actual_vec > 0] <-
        pumping_vec[starts_actual_vec > 0] *
        (equation(elapsed_time = starts_actual_vec[starts_actual_vec > 0],
                  distance = distance,
                  stor_coef = stor_coef,
                  transmissivity = transmissivity,
                  leakance = leakance) -
           equation(elapsed_time = stops_actual_vec[starts_actual_vec > 0],
                    distance = distance,
                    stor_coef = stor_coef,
                    transmissivity = transmissivity,
                    leakance = leakance))
      
      depletions_mat <- matrix(depletions_vec,
                               nrow = length(timesteps),
                               ncol = length(start_pumping))
      depletions_mat <- depletions_mat[-c(1), ]
      depletions <- base::rowSums(depletions_mat)
      
      
      
      
      
      
      fractional_vec[starts_actual_vec > 0] <-
        (equation(elapsed_time = starts_actual_vec[starts_actual_vec > 0],
                  distance = distance,
                  stor_coef = stor_coef,
                  transmissivity = transmissivity,
                  leakance = leakance) -
           equation(elapsed_time = stops_actual_vec[starts_actual_vec > 0],
                    distance = distance,
                    stor_coef = stor_coef,
                    transmissivity = transmissivity,
                    leakance = leakance))
      
      fractional_mat <- matrix(fractional_vec,
                               nrow = length(timesteps),
                               ncol = length(start_pumping))
      fractional_mat <- fractional_mat[-c(1), ]
      fractions <- base::rowSums(fractional_mat)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # memory allocation
      rm(list = c('depletions_vec','depletions_mat',
                  'fractional_mat','fractional_vec',
                  'starts_actual','stops_actual',
                  'timestep_mat','starts_mat',
                  'stops_mat','pumping_mat'))
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # sum and return
      Jenk_SDF <- (distance*distance*stor_coef)/transmissivity
      return(list(depletions,
                  fractions,
                  Jenk_SDF))
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    #-------------------------------------------------------------------------------
    # for each reach calculate sum of all depletions
    depletions_per_reach <- list()
    pump_frac_per_reach <- list()
    depletions_potential_per_reach <- list()
    jenk_sdf_per_reach <- list()
    custom_sdf_per_reach <- list()
    for(i in 1:ncol(closest_points_per_segment)){
      #-------------------------------------------------------------------------------
      if(suppress_loading_bar == FALSE){
        #-------------------------------------------------------------------------------
        # user message
        loading_bar(iter = i,
                    total = ncol(closest_points_per_segment),
                    width = 50,
                    optional_text = '')
        #-------------------------------------------------------------------------------
      }
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      # what are the closest points to each well
      points <- as.vector(unlist(closest_points_per_segment[, i]))
      fracs <- as.vector(unlist(reach_impact_frac[ , i]))
      RN <- st_drop_geometry(stream_points_geometry[points, stream_id_key])
      RN <- as.vector(unlist(RN))
      well_indices <- c(1:length(points))
      rm <- which(is.na(points))
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # if no wells are associated with the reach then pass
      if(all(is.na(points) == TRUE) == FALSE){
        #-------------------------------------------------------------------------------
        # remove any non-relevant wells from the for loop
        if(length(rm) > 0){
          well_indices <- well_indices[-c(rm)]
        } else {}
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        depletions_per_well <- list()
        pump_frac_per_well <- list()
        depletions_potential_per_well <- list()
        Jenk_SDF_per_well <- list()
        custom_SDF_per_well <- list()
        counter <- 0
        for(j in well_indices){
          #-------------------------------------------------------------------------------
          # increment list counter
          counter <- counter + 1
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get distances
          distance <- st_distance(wells[j, ],
                                  stream_points_geometry[points[j], ])
          distance <- as.numeric(distance)
          distances[[counter]] <- distance
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # get aquifer properties depending on if its assigned on a per-well basis
          # or if its going to be taken as an average of the aquifer properties along a line
          # between the well and the nearest stream
          if(is.null(model_grid) == TRUE){
            transmissivity <- as.numeric(st_drop_geometry(wells[j,well_transmissivity_key]))
            stor_coef <- as.numeric(st_drop_geometry(wells[j, well_stor_coef_key]))
          } else {
            #-------------------------------------------------------------------------------
            # make line between well and stream
            m <- matrix(c(st_coordinates(wells[j, ])[,1],
                          st_coordinates(wells[j, ])[,2],
                          st_coordinates(stream_points_geometry[points[j], ])[,1],
                          st_coordinates(stream_points_geometry[points[j], ])[,2]),
                        ncol = 2,
                        byrow = TRUE)
            line <- st_sf(st_sfc(st_linestring(m), crs = st_crs(wells)))
            st_geometry(line) <- 'geometry'
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # selecting correct gridcells
            grid_layers <- as.vector(unlist(st_drop_geometry(model_grid[,grid_layer_key])))
            well_layers <- as.vector(unlist(st_drop_geometry(wells[,well_layer_key])))
            gr <- model_grid[grid_layers == well_layers[j], ]
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # intersect line between well and stream with the grid
            int <- st_intersects(gr, line$geometry)
            grid_inds <- c(1:length(int))
            rm <- which(lengths(int) == 0)
            if(length(rm) > 0){
              grid_inds <- grid_inds[-c(rm)]
            } else {}
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            # getting average properties from grid
            gr <- gr[grid_inds, ]
            transmissivity <- mean(as.vector(unlist(st_drop_geometry(gr[,grid_transmissivity_key]))), na.rm = T)
            stor_coef <- mean(as.vector(unlist(st_drop_geometry(gr[,grid_stor_coef_key]))), na.rm = T)
            #-------------------------------------------------------------------------------
          }
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # take weighted mean of leakance along the reach, closer points are weighted more
          # logic is that closer points will control more of the flow
          reaches <- as.vector(unlist(st_drop_geometry(stream_points_geometry[,stream_id_key])))
          stream_inds <- reaches == RN[j]
          all_distances <- st_distance(wells[j, ],
                                       stream_points_geometry[stream_inds, ])
          all_distances <- c(1:length(all_distances))[order(as.numeric(all_distances))]
          
          leakance <- as.vector(unlist(st_drop_geometry(stream_points_geometry[stream_inds, leakance_key])))
          
          leakance <- weighted_mean(x = leakance, w = all_distances, na.rm = T)
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          # calculate maximum stream depletion potential and
          # multiply by fraction of depletions of this well apportioned to this reach
          Q_out <- hantush_stream_depletion_model(stor_coef = stor_coef,
                                                  transmissivity = transmissivity,
                                                  distance = distance,
                                                  QW = pumping[j, ],
                                                  leakance = leakance)
          Q_final <- Q_out[[1]]*fracs[j]
          Q_fraction <- Q_out[[2]]
          Jenk_SDF <- Q_out[[3]]
          #-------------------------------------------------------------------------------
          
          
          #-------------------------------------------------------------------------------
          # gradient descent example for future use
          # if(is.null(custom_sdf_time) == FALSE){
          #   custom_SDF <- custom_sdf_gradient_descent(analytical_model = analytical_model,
          #                                             distance = distance,
          #                                             stor_coef = stor_coef,
          #                                             transmissivity = transmissivity,
          #                                             leakance = leakance,
          #                                             custom_sdf_time = custom_sdf_time)
          #   custom_SDF_per_well[[counter]] <- custom_SDF
          # }
          #-------------------------------------------------------------------------------
          
          #-------------------------------------------------------------------------------
          depletions_per_well[[counter]] <- Q_final
          depletions_potential_per_well[[counter]] <- Q_fraction
          pump_frac_per_well[[counter]] <- pumping[j, ] * fracs[j]
          Jenk_SDF_per_well[[counter]] <- Jenk_SDF
          #-------------------------------------------------------------------------------
        }
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        depletions_total <- do.call(cbind, depletions_per_well)
        pump_frac_per_well_total <- do.call(cbind, pump_frac_per_well)
        depletions_potential_per_well_total <- do.call(cbind, depletions_potential_per_well)
        Jenk_SDF_per_well_total <- do.call(rbind, Jenk_SDF_per_well)
        #-------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------
        average_fractional_depletions <- calculate_depletion_potential(depletion_potential_criteria = depletion_potential_criteria,
                                                                       depletions_potential_per_well_total = depletions_potential_per_well_total,
                                                                       distances = distances,
                                                                       fracs = fracs,
                                                                       pumping = pumping)
        average_Jenk_SDF <- calculate_average_sdf(sdf_averaging_criteria = sdf_averaging_criteria,
                                                  sdf_vec = Jenk_SDF_per_well_total,
                                                  fracs = fracs,
                                                  pumping = pumping)
        
        jenk_sdf_per_reach[[i]] <- average_Jenk_SDF
        depletions_potential_per_reach[[i]] <- average_fractional_depletions
        depletions_per_reach[[i]] <- base::rowSums(depletions_total)
        pump_frac_per_reach[[i]] <- base::rowSums(pump_frac_per_well_total)
        if(is.null(custom_sdf_time) == FALSE){
          
          custom_sdf_per_reach[[i]] <- calculate_custom_sdf_time(average_fractional_depletions,
                                                                 target = custom_sdf_time)
        }
        #-------------------------------------------------------------------------------

        
        #-------------------------------------------------------------------------------
        # if(is.null(custom_sdf_time) == FALSE){
        #   custom_SDF_per_well_total <- do.call(rbind, custom_SDF_per_well)
        #   custom_SDF_per_well_total[custom_SDF_per_well_total == -9999] <- NA
        #   if(all(is.na(custom_SDF_per_well_total)) == FALSE){
        #     average_custom_sdf <- calculate_average_sdf(sdf_averaging_criteria = sdf_averaging_criteria,
        #                                                 sdf_vec = custom_SDF_per_well_total,
        #                                                 fracs = fracs,
        #                                                 pumping = pumping)
        #     custom_sdf_per_reach[[i]] <- average_custom_sdf
        #   } else {
        #     custom_sdf_per_reach[[i]] <- -9999
        #   }
        # }
        #-------------------------------------------------------------------------------
      } else{
        depletions_per_reach[[i]] <- rep(0, ncol(pumping)) # reach has no depletions
        pump_frac_per_reach[[i]] <- rep(0, ncol(pumping))
        depletions_potential_per_reach[[i]] <- rep(0, ncol(pumping))
        jenk_sdf_per_reach[[i]] <- NA
        
        if(is.null(custom_sdf_time) == FALSE){
          custom_sdf_time[[i]] <- NA
        }
      }
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    
    
    #-------------------------------------------------------------------------------
    # stats
    start_of_depletions <- lapply(depletions_per_reach, function(x){
      rle(x)$lengths[1]
    })
    start_of_depletions <- unlist(start_of_depletions)
    rm <- which(start_of_depletions == ncol(pumping))
    if(length(rm) > 0){
      start_of_depletions[-c(rm)] # if never started remove
    } else {}
    
    
    
    mean_start_of_depletions <- mean(start_of_depletions, na.rm = TRUE)
    median_start_of_depletions <- median(start_of_depletions, na.rm = TRUE)
    
    final_depletions <- lapply(depletions_per_reach, function(x){
      tail(x, 1)
    })
    n_timesteps <- ncol(pumping)
    mean_final_depletions <- mean(unlist(final_depletions), na.rm = TRUE)
    median_final_depletions <- median(unlist(final_depletions), na.rm = TRUE)
    which_max_final_depletions <- which.max(unlist(final_depletions))
    max_final_depletions <- max(unlist(final_depletions), na.rm = TRUE)
    #-------------------------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------------------------
    # write status to log
    u <- paste0('(',units,'^3','):')
    
    
    writeLines(text = sprintf('%s %s',
                              'Mean | Median start of stream depletions (timestep): ',
                              paste(round(mean_start_of_depletions,2),
                                    ' | ',
                                    median_start_of_depletions)),
               con = log_file)
    
    
    writeLines(text = sprintf('%s %s',
                              paste('Mean | Median final depletions',u),
                              paste(round(mean_final_depletions,4),
                                    '|',
                                    round(median_final_depletions,4),
                                    'at timestep (t_final)',
                                    n_timesteps)),
               con = log_file)
    
    writeLines(text = sprintf('%s %s',
                              paste('Max final depletions',u),
                              paste(round(max_final_depletions,4),
                                    'at timestep (t_final)',
                                    n_timesteps,
                                    'for reach',
                                    which_max_final_depletions)),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    
    
    
    
    
    
    #-------------------------------------------------------------------------------
    # output
    depletions_potential_per_reach <- do.call(rbind, depletions_potential_per_reach)
    depletions_per_reach <- do.call(rbind, depletions_per_reach)
    pump_frac_per_reach <- do.call(rbind, pump_frac_per_reach)
    jenk_sdf_per_reach <- do.call(rbind, jenk_sdf_per_reach)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # output and stats 2
    writeLines(text = sprintf('%s %s',
                              paste('Mean | Median Jenkins SDF'),
                              paste(round(mean(jenk_sdf_per_reach, na.rm = TRUE),4),
                                    '|',
                                    round(median(jenk_sdf_per_reach, na.rm = TRUE),4))),
               con = log_file)
    if(is.null(custom_sdf_time) == TRUE){
      custom_sdf_per_reach <- NULL
    } else{
      custom_sdf_per_reach <- do.call(rbind, custom_sdf_per_reach)
      writeLines(text = sprintf('%s %s',
                                paste('Mean | Median Custom SDF',paste0('(',custom_sdf_time,')')),
                                paste(round(mean(custom_sdf_per_reach, na.rm = TRUE),4),
                                      '|',
                                      round(median(custom_sdf_per_reach, na.rm = TRUE),4))),
                 con = log_file)
    }
    #-------------------------------------------------------------------------------
    return(list(depletions_per_reach,
                depletions_potential_per_reach,
                pump_frac_per_reach,
                jenk_sdf_per_reach,
                custom_sdf_per_reach))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ############################################################################################
  ######################################### MAIN FUNCTIONS ###################################
  ############################################################################################
  
  #===========================================================================================
  # Takes desired proximity criteria and decides which stream segments are
  # affected by each well
  #===========================================================================================
  find_impacted_stream_segments <- function(streams,
                                            wells,
                                            subwatersheds,
                                            proximity_criteria)
  {
    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s %s',
                              '####',
                              'Calculating Impacted Segments For Each Well'),
               con = log_file)
    
    writeLines(text = sprintf('%s %s',
                              'Well segment proximity criteria: ',
                               str_to_title(proximity_criteria)),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    
    
    
    #-------------------------------------------------------------------------------
    # if streams do not already come pre-processed as points of interest,
    # process them
    if(streams_are_points == FALSE){
      #-------------------------------------------------------------------------------
      # write status to log
      writeLines(text = sprintf('%s',
                                'Streams not preprocessed to points ... converting to points'),
                 con = log_file)
      #-------------------------------------------------------------------------------
      

      #-------------------------------------------------------------------------------
      # getting points from stream linestrings
      stream_points_list <- list()
      average_length <- list()
      id_list <- list()
      lambda_list <- list()
      leakance_list <- list()
      stream_transmissivity_list <- list()
      for(i in 1:nrow(streams)){
        coords <- Extract_SF_Linestring_Vertices(streams$geometry[i])
        stream_points_list[[i]] <- cbind(coords[[2]],
                                         coords[[1]])

        id_list[[i]] <- rep(st_drop_geometry(streams[i,stream_id_key]),
                            length(coords[[1]]))
        
        if(is.null(lambda_key) == FALSE){
          if(lambda_key %in% colnames(streams)){
            lambda_list[[i]] <- rep(st_drop_geometry(streams[i,lambda_key]),
                                    length(coords[[1]]))
          }
        }

        if(is.null(leakance_key) == FALSE){
          if(leakance_key %in% colnames(streams)){
            leakance_list[[i]] <- rep(st_drop_geometry(streams[i,leakance_key]),
                                      length(coords[[1]]))
          }
        }
        
        if(is.null(stream_transmissivity_key) == FALSE){
          if(stream_transmissivity_key %in% colnames(streams)){
            stream_transmissivity_list[[i]] <- rep(st_drop_geometry(streams[i,stream_transmissivity_key]),
                                                   length(coords[[1]]))
          }
        }
        
        average_length[[i]] <- length(coords[[1]])
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # turning into data frame for later export and diag
      stream_points_df <- do.call(rbind, stream_points_list)
      stream_points_df <- as.data.frame(stream_points_df)
      colnames(stream_points_df) <- c('x','y')
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # turning into geometry object for calculation of distances and intersections
      stream_points_geometry <- st_as_sf(stream_points_df,
                                         coords = c('x','y'),
                                         na.fail = FALSE,
                                         crs = st_crs(streams))
      cnams <- colnames(stream_points_geometry)
      stream_points_geometry <- cbind(unlist(id_list), stream_points_geometry)
      colnames(stream_points_geometry) <- c(stream_id_key, cnams)
      
      if(is.null(lambda_key) == FALSE){
        if(lambda_key %in% colnames(streams)){
          cnams <- colnames(stream_points_geometry)
          stream_points_geometry <- cbind(unlist(lambda_list), stream_points_geometry)
          colnames(stream_points_geometry) <- c(lambda_key, cnams)
        }
      }
      
      if(is.null(leakance_key) == FALSE){
        if(leakance_key %in% colnames(streams)){
          cnams <- colnames(stream_points_geometry)
          stream_points_geometry <- cbind(unlist(leakance_list), stream_points_geometry)
          colnames(stream_points_geometry) <- c(leakance_key, cnams)
        }
      }
      
      if(is.null(stream_transmissivity_key) == FALSE){
        if(stream_transmissivity_key %in% colnames(streams)){
          cnams <- colnames(stream_points_geometry)
          stream_points_geometry <- cbind(unlist(stream_transmissivity_list), stream_points_geometry)
          colnames(stream_points_geometry) <- c(stream_transmissivity_key, cnams)
        }
      }
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # save memory
      rm(stream_points_list)
      rm(id_list)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # write status to log
      average_length <- do.call(rbind, average_length)
      writeLines(text = sprintf('%s %s',
                                'Number of unique stream segments: ',
                                nrow(streams)),
                 con = log_file)
      writeLines(text = sprintf('%s %s',
                                'Mean | Median number of points per segment: ',
                                paste(round(mean(average_length, na.rm = TRUE),0),
                                      '|',
                                      round(median(average_length, na.rm = TRUE),0))),
                 con = log_file)
      rm(average_length) # save memory
      #-------------------------------------------------------------------------------
    } else {
      stream_points_geometry <- streams
      #-------------------------------------------------------------------------------
      # write status to log file
      writeLines(text = sprintf('%s',
                                'Streams came already preprocessed to points ... passing'),
                 con = log_file)
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------------------------
    # calculate which wells impact which stream points
    if(str_to_title(proximity_criteria) == 'Adjacent'){
      
      #-------------------------------------------------------------------------------
      # ensure everything is in the same projection
      proj_output <- ensure_projections(wells = wells,
                                        geometry_list = list(subwatersheds,
                                                             stream_points_geometry))
      subwatersheds <- proj_output[[1]]
      stream_points_geometry <- proj_output[[2]]
      rm(proj_output)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # find adjacent stream points
      writeout <- find_adjacent_stream_points(wells = wells,
                                              subwatersheds = subwatersheds,
                                              stream_points_geometry = stream_points_geometry)
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # calculate which wells impact which stream points
    if(str_to_title(proximity_criteria) == 'Adjacent+Expanding'){
      
      #-------------------------------------------------------------------------------
      # ensure everything is in the same projection
      proj_output <- ensure_projections(wells = wells,
                                        geometry_list = list(subwatersheds,
                                                             stream_points_geometry))
      subwatersheds <- proj_output[[1]]
      stream_points_geometry <- proj_output[[2]]
      rm(proj_output)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # find adjacent stream points
      writeout <- find_adjacent_and_expanding_stream_points(wells = wells,
                                                            subwatersheds = subwatersheds,
                                                            influence_radius = influence_radius,
                                                            stream_points_geometry = stream_points_geometry)
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------------------------
    # calculate which wells impact which stream points
    if(str_to_title(proximity_criteria) == 'Local Area' |
       str_to_title(proximity_criteria) == 'Expanding'){
      
      #-------------------------------------------------------------------------------
      # ensure everything is in the same projection
      proj_output <- ensure_projections(wells = wells,
                                        geometry_list = list(stream_points_geometry))
      stream_points_geometry <- proj_output[[1]]
      rm(proj_output)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # find adjacent stream points
      writeout <- find_local_stream_points(wells = wells,
                                           influence_radius = influence_radius,
                                           stream_points_geometry = stream_points_geometry)
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # calculate which wells impact which stream points
    if(str_to_title(proximity_criteria) == 'Whole Domain'){
      #-------------------------------------------------------------------------------
      # ensure everything is in the same projection
      proj_output <- ensure_projections(wells = wells,
                                        geometry_list = list(stream_points_geometry))
      stream_points_geometry <- proj_output[[1]]
      rm(proj_output)
      #-------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------
      # associate all stream points with all wells
      writeout <- find_whole_domain_points(wells = wells,
                                           stream_points_geometry = stream_points_geometry)
      #-------------------------------------------------------------------------------
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # log message
    m1 <- mean(as.numeric(writeout[[2]]$ImpLMet), na.rm = T)
    m2 <- median(as.numeric(writeout[[2]]$ImpLMet), na.rm = T)
    writeLines(text = sprintf('%s %s',
                              paste0('Mean | Median impacted segment length in ',units,': '),
                              paste0(as.character(round(m1,0)),
                                     '|',
                                     as.character(round(m2,0)))),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Found impacted segments without error'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Moving to next step ...'),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # return to higher level
    stream_points_geometry$Index <- c(1:nrow(stream_points_geometry))
    return(list(writeout[[1]],
                writeout[[2]],
                stream_points_geometry))
    #-------------------------------------------------------------------------------
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  #===========================================================================================
  # Apportions depletions based on given criteria
  #===========================================================================================
  calculate_depletion_apportionments <- function(wells,
                                                 apportionment_criteria,
                                                 stream_points_geometry,
                                                 stream_id_key){
    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s %s',
                              '####',
                              'Apportioning Depletions For Each Well'),
               con = log_file)
    
    writeLines(text = sprintf('%s %s',
                              'Well apportionment criteria: ',
                              str_to_title(apportionment_criteria)),
               con = log_file)
    
    if(geologic_apportionment == TRUE){
      writeLines(text = sprintf('%s',
                                'Geologic appoortionment set to TRUE'),
                 con = log_file)
    }
    #-------------------------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------------------------
    # find what points are important for each well
    if(str_to_title(apportionment_criteria) %in% c('Inverse Distance',
                                                   'Inverse Distance Squared',
                                                   'Thiessen Polygon',
                                                   'Web',
                                                   'Web Squared')){
      
      w_index <- as.vector(unlist(st_drop_geometry(wells[ ,wells_id_key])))
      closest_points_per_segment <- find_closest_points_per_segment(wells = wells,
                                                                    stream_points_geometry = stream_points_geometry,
                                                                    stream_id_key = stream_id_key)
      closest_points_per_segment <- cbind(w_index,
                                          closest_points_per_segment)
      closest_points_per_segment <- as.data.frame(closest_points_per_segment)
      colnames(closest_points_per_segment) <- c('wellN',
                                                paste0('PN',
                                                       1:(ncol(closest_points_per_segment)-1)))
      
    } else {
      closest_points_per_segment <- NULL # otherwise if some future method implemented where this is not important
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # Inverse distance apportionment
    if(str_to_title(apportionment_criteria) %in% c('Inverse Distance',
                                                   'Inverse Distance Squared')){
      if(str_to_title(apportionment_criteria) == 'Inverse Distance Squared'){
        power <- 2
      } else {
        power <- 1
      }

      output <- inverse_distance_apportionment(power = power,
                                               wells = wells,
                                               closest_points_per_segment = closest_points_per_segment,
                                               stream_points_geometry = stream_points_geometry)

    }
    #-------------------------------------------------------------------------------


    #-------------------------------------------------------------------------------
    # Inverse distance apportionment
    if(str_to_title(apportionment_criteria) %in% c('Thiessen Polygon')){

      output <- thiessen_polygon_apportionment(wells = wells,
                                               closest_points_per_segment = closest_points_per_segment,
                                               stream_points_geometry = stream_points_geometry)

    }
    #-------------------------------------------------------------------------------



    #-------------------------------------------------------------------------------
    # Inverse distance apportionment
    if(str_to_title(apportionment_criteria) %in% c('Web',
                                                   'Web Squared')){
      if(str_to_title(apportionment_criteria) == 'Web Squared'){
        power <- 2
      } else {
        power <- 1
      }

      output <- web_apportionment(power = power,
                                  wells = wells,
                                  closest_points_per_segment = closest_points_per_segment,
                                  stream_points_geometry = stream_points_geometry)

    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # log message
    writeLines(text = sprintf('%s',
                              'Apportioned depletions without error'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Moving to next step ...'),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    return(list(output[[1]],
                output[[2]],
                closest_points_per_segment))
  }
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  #===========================================================================================
  # Apportions depletions based on given criteria
  #===========================================================================================
  calculate_stream_depletion_per_reach <- function(closest_points_per_segment = closest_points_per_segment,
                                                   stream_points_geometry = stream_points_geometry,
                                                   wells = wells,
                                                   well_transmissivity_key = well_transmissivity_key,
                                                   well_stor_coef_key = well_stor_coef_key,
                                                   lambda_key = lambda_key,
                                                   leakance_key = leakance_key,
                                                   analytical_model = analytical_model){
    #-------------------------------------------------------------------------------
    # write status to log
    writeLines(text = sprintf('%s %s',
                              '####',
                              'Calculating Depletions For Each Well'),
               con = log_file)
    
    writeLines(text = sprintf('%s %s',
                              'Using analytical model: ',
                              str_to_title(analytical_model)),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    
    
    #-------------------------------------------------------------------------------
    # find what points are important for each well
    if(str_to_title(analytical_model) %in% c('Glover')){
      output <- glover_stream_depletion_calculations(closest_points_per_segment = closest_points_per_segment,
                                                     reach_impact_frac = reach_impact_frac,
                                                     stream_points_geometry = stream_points_geometry,
                                                     wells = wells,
                                                     well_transmissivity_key = well_transmissivity_key,
                                                     well_stor_coef_key = well_stor_coef_key)

    } else if (str_to_title(analytical_model) %in% c('Hunt')){
      output <- hunt_stream_depletion_calculations(closest_points_per_segment = closest_points_per_segment,
                                                   stream_points_geometry = stream_points_geometry,
                                                   reach_impact_frac = reach_impact_frac,
                                                   wells = wells,
                                                   well_transmissivity_key = well_transmissivity_key,
                                                   well_stor_coef_key = well_stor_coef_key,
                                                   lambda_key = lambda_key)

    } else if (str_to_title(analytical_model) %in% c('Hantush')){
      output <- hantush_stream_depletion_calculations(closest_points_per_segment = closest_points_per_segment,
                                                      stream_points_geometry = stream_points_geometry,
                                                      reach_impact_frac = reach_impact_frac,
                                                      wells = wells,
                                                      well_transmissivity_key = well_transmissivity_key,
                                                      well_stor_coef_key = well_stor_coef_key,
                                                      leakance_key = leakance_key)

    } else {}
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # format output
    depletions_per_reach <- cbind(as.vector(unlist(st_drop_geometry(streams[,stream_id_key]))),
                                  output[[1]])
    depletions_per_reach <- as.data.frame(depletions_per_reach)
    colnames(depletions_per_reach) <- c('RN', paste0('T',1:(ncol(depletions_per_reach)-1)))
    
    
    
    depletions_potential_per_reach <- cbind(as.vector(unlist(st_drop_geometry(streams[,stream_id_key]))),
                                            output[[2]])
    depletions_potential_per_reach <- as.data.frame(depletions_potential_per_reach)
    colnames(depletions_potential_per_reach) <- c('RN', paste0('T',1:(ncol(depletions_potential_per_reach)-1)))
    
    
    
    pump_frac_per_reach <- cbind(as.vector(unlist(st_drop_geometry(streams[,stream_id_key]))),
                                  output[[3]])
    pump_frac_per_reach <- as.data.frame(pump_frac_per_reach)
    colnames(pump_frac_per_reach) <- c('RN', paste0('T',1:(ncol(pump_frac_per_reach)-1)))
    
    
    
    jenk_sdf_per_reach <- cbind(as.vector(unlist(st_drop_geometry(streams[,stream_id_key]))),
                                 output[[4]])
    jenk_sdf_per_reach <- as.data.frame(jenk_sdf_per_reach)
    colnames(jenk_sdf_per_reach) <- c('RN', 'T')
    
    
    if(is.null(custom_sdf_time) == FALSE){
      custom_sdf_per_reach <- cbind(as.vector(unlist(st_drop_geometry(streams[,stream_id_key]))),
                                    output[[5]])
      custom_sdf_per_reach <- as.data.frame(custom_sdf_per_reach)
      colnames(custom_sdf_per_reach) <- c('RN', 'T')
    } else {
      custom_sdf_per_reach <- output[[5]]
    }
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # log message
    writeLines(text = sprintf('%s',
                              'Calculated depletions without error'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Finished ...'),
               con = log_file)
    #-------------------------------------------------------------------------------
    
    return(list(depletions_per_reach,
                depletions_potential_per_reach,
                pump_frac_per_reach,
                jenk_sdf_per_reach,
                custom_sdf_per_reach))
  }
  #-------------------------------------------------------------------------------
  
  
  


  
  
  
  
  
  
  ############################################################################################
  ######################################### RUN FUNCTIONS ####################################
  ############################################################################################
  required_packages <- c('sf','sp','raster','terra','lubridate','stringr','Rmpfr')
  for(i in 1:length(required_packages)){
    require_package(required_packages[i])
  }

  units <- units(st_distance(wells[1, ], wells[1, ]))$numerator
  if(units == 'm'){
    units <- 'meters'
  } else if (units == 'B0'){
    units <- 'degrees'
  } else if (units == 'US_survey_foot'){
    units <- 'feet'
  }
  
  #-------------------------------------------------------------------------------
  # open log file to write program execution
  log_file <- file(file.path(diag_out_dir,'log.txt'), 'w')
  writeLines(text = sprintf('%s %s',
                            '####',
                            'Opening program to calculate stream depletion'),
             con = log_file)
  writeLines(text = '',
             con = log_file)
  #-------------------------------------------------------------------------------
  
  
  
  ############################################################################################
  # errors
  #-------------------------------------------------------------------------------
  # before proceeding to well calculations do I have the information I need?
  if(str_to_title(proximity_criteria) %in% c('Adjacent','Adjacent+Expanding') &
     is.null(subwatersheds) == TRUE){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              'Proximity criteria required subwatersheds but none supplied'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'proximity criteria required subwatersheds (Adjacent | Adjacent+Expanding)\n',
                'but none supplied\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  } else if (str_to_title(proximity_criteria) %in% c('Local Area',
                                                     'Expanding',
                                                     'Adjacent+Expanding') &
             is.null(influence_radius) == TRUE){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              'Proximity criteria required influence radius but none supplied'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'If using Local Area critiera please calculate influence radius as per:'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Zipper et al. (2019) https://doi.org/10.1029/2018WR024403'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'proximity criteria required influence radius (Local Area | Expanding | Adjacent+Expanding)\n',
                'but none supplied\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  } else {}
  
  if(streams_are_points == TRUE &
     is.null(stream_id_key) == TRUE){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              'Identifying column for streams required to calculate impacted length but none supplied'),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'Identifying column for streams required to calculate impacted length',
                'but none supplied\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  }
  
  if(str_to_title(analytical_model) == 'Glover' &
     is.null(well_stor_coef_key) == TRUE  |
     is.null(well_transmissivity_key) == TRUE){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              paste0('Identifying column for storage coefficient or transmissivity in ',
                                     str_to_title(analytical_model),
                                     ' model required but not present in well set')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'Identifying column for storage coefficient or transmissivity in ',
                str_to_title(analytical_model),
                ' model required but not present in well set\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  }
  
  
  if(!well_stor_coef_key %in% colnames(wells) |
     !well_transmissivity_key %in% colnames(wells)){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              paste0('Identifying column for storage coefficient or transmissivity in ',
                                     'all models',
                                     ' required but not present in well set')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'Identifying column for storage coefficient or transmissivity in ',
                'all models',
                ' required but not present in well set\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  }
  
  
  if(is.null(model_grid) == FALSE){
    if(is.null(well_layer_key) == TRUE){
      #-------------------------------------------------------------------------------
      writeLines(text = sprintf('%s',
                                paste0('Identifying column for layer required ',
                                       'when model grid passed as argument',
                                       ' but not present in well set')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                'Exiting program ...'),
                 con = log_file)
      close(log_file)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                  'Identifying column for layer required ',
                  'when model grid passed as argument',
                  ' but not present in well set\n',
                  'exiting program ...'))
      #-------------------------------------------------------------------------------
    } else if(!well_layer_key %in% colnames(wells)){
      #-------------------------------------------------------------------------------
      writeLines(text = sprintf('%s',
                                paste0('Identifying column for layer required ',
                                       'when model grid passed as argument',
                                       ' but not present in well set')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                'Exiting program ...'),
                 con = log_file)
      close(log_file)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                  'Identifying column for layer required ',
                  'when model grid passed as argument',
                  ' but not present in well set\n',
                  'exiting program ...'))
      #-------------------------------------------------------------------------------
    }
  }
  

  
  if(str_to_title(analytical_model) == 'Hunt' &
     is.null(lambda_key) == TRUE){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              paste0('Identifying column for lambda in ',
                                     str_to_title(analytical_model),
                                     ' model required but not present in streams set')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'Identifying column for lambda in ',
                str_to_title(analytical_model),
                ' model required but not present in streams set\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  }
  
  
  
  
  if(str_to_title(analytical_model) == 'Hunt' &
     is.null(lambda_key) == FALSE){
    if(!lambda_key %in% colnames(streams) == TRUE){
      #-------------------------------------------------------------------------------
      writeLines(text = sprintf('%s',
                                paste0('Identifying column for lambda in',
                                       str_to_title(analytical_model),
                                       ' model required but not present in streams set')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                'Exiting program ...'),
                 con = log_file)
      close(log_file)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                  'Identifying column for lambda in ',
                  str_to_title(analytical_model),
                  ' model required but not present in streams set\n',
                  'exiting program ...'))
      #-------------------------------------------------------------------------------
    }
  }
    
  
  if(str_to_title(analytical_model) == 'Hantush' &
     is.null(leakance_key) == TRUE){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              paste0('Identifying column for leakance in ',
                                     str_to_title(analytical_model),
                                     ' model required but not present in streams set')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'Identifying column for leakance in ',
                str_to_title(analytical_model),
                ' model required but not present in streams set\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  }
  
  
  if(str_to_title(analytical_model) == 'Hantush' &
     is.null(leakance_key) == FALSE){
    if(!leakance_key %in% colnames(streams) == TRUE){
      #-------------------------------------------------------------------------------
      writeLines(text = sprintf('%s',
                                paste0('Identifying column for leakance in',
                                       str_to_title(analytical_model),
                                       ' model required but not present in streams set')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                'Exiting program ...'),
                 con = log_file)
      close(log_file)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                  'Identifying column for leakance in ',
                  str_to_title(analytical_model),
                  ' model required but not present in streams set\n',
                  'exiting program ...'))
      #-------------------------------------------------------------------------------
    }
  }
  
  
  if(is.null(model_grid) == FALSE){
    proj_output <- ensure_projections(wells = wells,
                                      geometry_list = list(model_grid))
    model_grid <- proj_output[[1]]
    
    if(!grid_layer_key %in% colnames(model_grid) |
       !grid_stor_coef_key %in% colnames(model_grid) |
       !grid_transmissivity_key %in% colnames(model_grid)){
      #-------------------------------------------------------------------------------
      writeLines(text = sprintf('%s',
                                paste0('Identifiers passed for storage coefficient, transmissivity, or layer')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                paste0('not found in model grid column names.')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                'Exiting program ...'),
                 con = log_file)
      close(log_file)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                  'Identifiers paseed for storage coefficient,transmissivity, or layer\n',
                  'not found in model grid column names.\n',
                  'exiting program ...'))
      #-------------------------------------------------------------------------------
    }
  }
  
  
  if(str_to_title(apportionment_criteria) == 'Thiessen Polygon' &
     nrow(streams) == 1){
    #-------------------------------------------------------------------------------
    writeLines(text = sprintf('%s',
                              paste0('Apportionment criteria of ',
                                     apportionment_criteria,
                                     ' selected but only one stream present')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              paste0(apportionment_criteria,'s',
                                     ' cannot be made out of only one entity')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'Exiting program ...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                'Apportionment criteria of ',
                apportionment_criteria,
                ' selected but only one stream present\n',
                apportionment_criteria,'s',
                ' cannot be made out of only one entity\n',
                'please select another apportionment method\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  }
    
    
  if(geologic_apportionment == TRUE &
     is.null(stream_transmissivity_key) == FALSE){
    if(!stream_transmissivity_key %in% colnames(streams) == TRUE){
      #-------------------------------------------------------------------------------
      writeLines(text = sprintf('%s',
                                paste0('Geologic apportionment selected but stream transmissivity key of \'',
                                       stream_transmissivity_key,
                                       '\' not present in stream set')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                'Exiting program ...'),
                 con = log_file)
      close(log_file)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                  'Geologic apportionment selected but stream transmissivity key of \'',
                  stream_transmissivity_key,
                  '\' not present in stream set\n',
                  'exiting program ...'))
      #-------------------------------------------------------------------------------
    }
  }
    
    
  if(geologic_apportionment == TRUE &
     is.null(stream_transmissivity_key) == TRUE){
    if(is.null(model_grid) == TRUE){
      #-------------------------------------------------------------------------------
      writeLines(text = sprintf('%s',
                                paste0('Geologic apportionment selected but neither stream transmissivity key')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                paste0('nor model grid given.')),
                 con = log_file)
      writeLines(text = sprintf('%s',
                                'Exiting program ...'),
                 con = log_file)
      close(log_file)
      #-------------------------------------------------------------------------------
      
      
      #-------------------------------------------------------------------------------
      stop(paste0('\ncalculate_stream_depletions.R encountered Error:    \n',
                  'Geologic apportionment selected but stream transmissivity key not given\n',
                  'exiting program ...'))
      #-------------------------------------------------------------------------------
    }
  }  
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  ############################################################################################
  # run stream impacted segments
  
  #-------------------------------------------------------------------------------
  # capture any error output and write to log file
  tryCatch(expr = {
    #-------------------------------------------------------------------------------
    # user message
    if(suppress_console_messages == FALSE){
      cat('Finding which stream reaches are impacted by wells: Step (1/3)')
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # find impacted points by proximity criteria
    if(is.null(wells_id_key) == TRUE){
      wells_id_key <- 'ID'
      wells$ID <- c(1:nrow(wells))
    } else {}
    
    if(is.null(stream_id_key) == TRUE){
      stream_id_key <- 'ID'
      streams$ID <- c(1:nrow(streams))
    } else {}
    
    output <- find_impacted_stream_segments(streams,
                                            wells,
                                            subwatersheds,
                                            proximity_criteria)
    impacted_points <- output[[1]]
    wells <- output[[2]]
    stream_points_geometry <- output[[3]]
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    # writeout
    write.csv(impacted_points,
              file.path(data_out_dir,
                        'impacted_stream_points_by_well.csv'),
              row.names = FALSE)

    st_write(wells,
              file.path(data_out_dir,
                        'wells_with_impacted_length.shp'),
              append = FALSE,
              quiet = TRUE)

    st_write(stream_points_geometry,
             file.path(data_out_dir,
                       'stream_points.shp'),
             append = FALSE,
             quiet = TRUE)
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # save space
    rm(output)
    writeLines(text = sprintf('%s',
                              ''),
               con = log_file)
    #-------------------------------------------------------------------------------
  }, error = function(e){
    #-------------------------------------------------------------------------------
    # write error to log file
    status <- 'find_impacted_stream_segments'
    writeLines(text = sprintf('%s %s',
                              'ENCOUNTERED ERROR: ',
                              class(e)[1]),
               con = log_file)
    writeLines(text = sprintf('%s %s',
                              'ON COMMAND: ',
                              paste0(capture.output(e$call),collapse = ' ')),
               con = log_file)
    writeLines(text = sprintf('%s %s',
                              'FOR REASON: ',
                              paste0(e$message, collapse = ' ')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'exiting program...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # write error to console
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    ', class(e)[1],'\n',
                'during:    ', status,'\n',
                'on command:    ', paste0(capture.output(e$call),collapse = ' '),'\n',
                'for reason:    ', paste0(e$message, collapse = ' '),'\n',
                'for more information see the log.txt file output\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  })
  #-------------------------------------------------------------------------------
  
  
  
  
  
  
  ############################################################################################
  # run depletion apportionments
  
  #-------------------------------------------------------------------------------
  # capture any error output and write to log file
  tryCatch(expr = {
    #-------------------------------------------------------------------------------
    # user message
    if(suppress_console_messages == FALSE){
      cat('\nApportioning depletions per well: Step (2/3)')
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # find impacted points by proximity criteria
    output <- calculate_depletion_apportionments(stream_points_geometry = stream_points_geometry,
                                                 stream_id_key = stream_id_key,
                                                 wells = wells,
                                                 apportionment_criteria = apportionment_criteria)
    reach_impact_frac <- output[[1]]
    impacted_reaches <- output[[2]]
    closest_points_per_segment <- output[[3]]
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # writeout
    write.csv(reach_impact_frac,
              file.path(data_out_dir,
                        'reach_depletion_fraction_by_well.csv'),
              row.names = FALSE)

    write.csv(impacted_reaches,
              file.path(data_out_dir,
                        'impacted_reaches_by_well.csv'),
              row.names = FALSE)
    
    if(str_to_title(apportionment_criteria) %in% c('Inverse Distance',
                                                   'Inverse Distance Squared',
                                                   'Thiessen Polygon',
                                                   'Web',
                                                   'Web Squared')){
      write.csv(closest_points_per_segment,
                file.path(data_out_dir,
                          'closest_points_per_reach_per_well.csv'),
                row.names = FALSE)
    }
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    # subtracting a column from a two column data frame turns it into a vector
    # when this behavior is not wanted, this is a workaround for this special case
    # of only one well being present
    if(ncol(reach_impact_frac) < 3){
      cnams <- colnames(reach_impact_frac)
      m <- matrix(data = as.vector(unlist(reach_impact_frac[,-c(1)])),
                  nrow = nrow(reach_impact_frac),
                  ncol = ncol(reach_impact_frac) - 1,
                  byrow = TRUE)
      reach_impact_frac <- as.data.frame(m)
      colnames(reach_impact_frac) <- cnams[-c(1)]
      
      
      cnams <- colnames(impacted_reaches)
      m <- matrix(data = as.vector(unlist(impacted_reaches[,-c(1)])),
                  nrow = nrow(impacted_reaches),
                  ncol = ncol(impacted_reaches) - 1,
                  byrow = TRUE)
      impacted_reaches <- as.data.frame(m)
      colnames(impacted_reaches) <- cnams[-c(1)]
      
      
      cnams <- colnames(closest_points_per_segment)
      m <- matrix(data = as.vector(unlist(closest_points_per_segment[,-c(1)])),
                  nrow = nrow(closest_points_per_segment),
                  ncol = ncol(closest_points_per_segment) - 1,
                  byrow = TRUE)
      closest_points_per_segment <- as.data.frame(m)
      colnames(closest_points_per_segment) <- cnams[-c(1)]
    } else {
      reach_impact_frac <- reach_impact_frac[,-c(1)]
      impacted_reaches <- impacted_reaches[,-c(1)]
      closest_points_per_segment <- closest_points_per_segment[,-c(1)]
    }
    #-------------------------------------------------------------------------------
    

    #-------------------------------------------------------------------------------
    # save space
    rm(output)
    writeLines(text = sprintf('%s',
                              ''),
               con = log_file)
    #-------------------------------------------------------------------------------
  }, error = function(e){
    #-------------------------------------------------------------------------------
    # write error to log file
    status <- 'calculate_depletion_apportionments'
    writeLines(text = sprintf('%s %s',
                              'ENCOUNTERED ERROR: ',
                              class(e)[1]),
               con = log_file)
    writeLines(text = sprintf('%s %s',
                              'ON COMMAND: ',
                              paste0(capture.output(e$call),collapse = ' ')),
               con = log_file)
    writeLines(text = sprintf('%s %s',
                              'FOR REASON: ',
                              paste0(e$message, collapse = ' ')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'exiting program...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # write error to console
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    ', class(e)[1],'\n',
                'during:    ', status,'\n',
                'on command:    ', paste0(capture.output(e$call),collapse = ' '),'\n',
                'for reason:    ', paste0(e$message, collapse = ' '),'\n',
                'for more information see the log.txt file output\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  })
  #-------------------------------------------------------------------------------
  
  
  

  
  
  
  
  
  ############################################################################################
  # run calculate depletions
  
  #-------------------------------------------------------------------------------
  # capture any error output and write to log file
  tryCatch(expr = {
    #-------------------------------------------------------------------------------
    # user message
    if(suppress_console_messages == FALSE){
      cat('\nCalculating depletions per reach: Step (3/3)')
    }
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    # calculate streamflow depletions based on either
    # a volumetric or fractional approach
    # fractional gives number between 0 and 1, where 1 is depletion is equal to
    # the pumping rate
    output <- calculate_stream_depletion_per_reach(closest_points_per_segment = closest_points_per_segment,
                                                   stream_points_geometry = stream_points_geometry,
                                                   wells = wells,
                                                   well_transmissivity_key = well_transmissivity_key,
                                                   well_stor_coef_key = well_stor_coef_key,
                                                   lambda_key = lambda_key,
                                                   leakance_key = leakance_key,
                                                   analytical_model = analytical_model)
    depletions_by_reach <- output[[1]]
    depletions_potential_by_reach <- output[[2]]
    pump_frac_by_reach <- output[[3]]
    jenk_sdf_by_reach <- output[[4]]
    custom_sdf_by_reach <- output[[5]]
    #-------------------------------------------------------------------------------

    
    #-------------------------------------------------------------------------------
    write.csv(depletions_by_reach,
              file.path(data_out_dir,
                        paste0('volumetric_depletions_by_reach.csv')),
              row.names = FALSE)
    write.csv(depletions_potential_by_reach,
              file.path(data_out_dir,
                        paste0('fractional_depletions_by_reach.csv')),
              row.names = FALSE)
    write.csv(pump_frac_by_reach,
              file.path(data_out_dir,
                        paste0('pump_frac_by_reach.csv')),
              row.names = FALSE)
    write.csv(jenk_sdf_by_reach,
              file.path(data_out_dir,
                        paste0('jenk_sdf_by_reach.csv')),
              row.names = FALSE)
    
    if(is.null(custom_sdf_time) == FALSE){
      write.csv(custom_sdf_by_reach,
                file.path(data_out_dir,
                          paste0('custom_sdf_by_reach.csv')),
                row.names = FALSE)
    }
    #-------------------------------------------------------------------------------

    
    #-------------------------------------------------------------------------------
    # save space
    rm(output)
    #-------------------------------------------------------------------------------
  }, error = function(e){
    #-------------------------------------------------------------------------------
    # write error to log file
    status <- 'calculate_depletions'
    writeLines(text = sprintf('%s %s',
                              'ENCOUNTERED ERROR: ',
                              class(e)[1]),
               con = log_file)
    writeLines(text = sprintf('%s %s',
                              'ON COMMAND: ',
                              paste0(capture.output(e$call),collapse = ' ')),
               con = log_file)
    writeLines(text = sprintf('%s %s',
                              'FOR REASON: ',
                              paste0(e$message, collapse = ' ')),
               con = log_file)
    writeLines(text = sprintf('%s',
                              'exiting program...'),
               con = log_file)
    close(log_file)
    #-------------------------------------------------------------------------------
    
    
    #-------------------------------------------------------------------------------
    # write error to console
    stop(paste0('\ncalculate_stream_depletions.R encountered Error:    ', class(e)[1],'\n',
                'during:    ', status,'\n',
                'on command:    ', paste0(capture.output(e$call),collapse = ' '),'\n',
                'for reason:    ', paste0(e$message, collapse = ' '),'\n',
                'for more information see the log.txt file output\n',
                'exiting program ...'))
    #-------------------------------------------------------------------------------
  })
  #-------------------------------------------------------------------------------
  
  
  
  #-------------------------------------------------------------------------------
  # close log file
  close(log_file)
  #-------------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------