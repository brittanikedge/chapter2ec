run_se <- function(data){
  
  #--------------------------
  # Spatial Error Model estimation
  #--------------------------
  cent <- getSpPPolygonsLabptSlots(as(data,"Spatial"))
  dnb <- dnearneigh(cent,0,15)
  dsts <- nbdists(dnb,cent)
  idw <- lapply(dsts, function(x) 1/x)
  wq <- nb2listw(dnb,glist=idw,style = "W",zero.policy = TRUE)
  se_res<-errorsarlm(yield~input_rate+I(input_rate*input_rate), 
                     data = data, listw=wq, zero.policy=TRUE, na.action=na.omit)
  
  return(se_res)
}

run_se_ec <- function(data){
  # data = analysis_res_g$data[[1]]
  # predict_results <- predict(se_res, data = data, listw = wq, zero.policy = TRUE)
  if(c("ec_0_2") %in% names(data)){
    ec_col_number <- match(c("ec_0_2"), names(data))
    if(is.na(ec_col_number) == FALSE){
      data <- rename(data, "ecs" = names(data)[ec_col_number])
    }
  }
  
  if(sum(ifelse(is.na(data$n) == TRUE, 1, 0)) > nrow(data)/2){
    #--------------------------
    # Spatial Error Model estimation
    #--------------------------
    cent <- getSpPPolygonsLabptSlots(as(data, "Spatial"))
    dnb <- dnearneigh(cent, 0, 15, longlat = FALSE, row.names = data$obs_id)
    dsts <- nbdists(dnb, cent)
    idw <- lapply(dsts, function(x) 1/x)
    wq <- nb2listw(dnb, glist = idw, style = "W", zero.policy = TRUE)
    se_res <- spatialreg::errorsarlm(yield  ~ s + I(s*s) + ecs + I(s*ecs), 
                                     data = data, listw=wq, zero.policy=TRUE, na.action=na.omit, tol.solve = 1e-14)
  }else{
    #--------------------------
    # Spatial Error Model estimation
    #--------------------------
    cent <- getSpPPolygonsLabptSlots(as(data, "Spatial"))
    dnb <- dnearneigh(cent, 0, 15, longlat = FALSE, row.names = data$obs_id)
    dsts <- nbdists(dnb, cent)
    idw <- lapply(dsts, function(x) 1/x)
    wq <- nb2listw(dnb, glist = idw, style = "W", zero.policy = TRUE)
    se_res <- spatialreg::errorsarlm(yield  ~ s + n + I(s*s) + I(n*n) + ecs + I(n*ecs) + I(s*ecs) + I(s*n), 
                                     data = data, listw=wq, zero.policy=TRUE, na.action=na.omit, tol.solve = 1e-14)
    
    
  }
  
  return(se_res)
}

run_scam_gam <- function(data, field_vars){
  
  results <- NULL
  
  results <- tryCatch(
    {
      withTimeout(
        {
          formula <- paste0(
            "yield ~ s(input_rate, bs = \"micv\") + s(X, k = 5) + s(Y, k = 5) + te(X, Y, k = c(5, 5))",
            ifelse(
              length(field_vars) > 0,
              paste0(" + ", paste0(field_vars, collapse = " + ")),
              ""
            )
          ) %>% formula()
          
          scam(formula, data = data)
          
        },
        timeout = 20, # 20 seconds,
        onTimeout = "silent"
      )
    },
    error = function(cond){
      return(NULL)
    }
  )
  
  
  if (is.null(results)) {
    
    formula <- paste0(
      "yield ~ s(input_rate, k = 3) + s(X, k = 4) + s(Y, k = 4) + te(X, Y, k = c(5, 5))",
      ifelse(
        length(field_vars) > 0,
        paste0(" + ", paste0(field_vars, collapse = " + ")),
        ""
      )
    ) %>% formula()
    
    results <- gam(formula, data = data)
    
  }
  
  return(results) 
  
}

eval_with_timeout <- function(expr, envir = parent.frame(), timeout, on_timeout = c("error", "warning", "silent")) {
  # substitute expression so it is not executed as soon it is used
  expr <- substitute(expr)
  
  # match on_timeout
  on_timeout <- match.arg(on_timeout)
  
  # execute expr in separate fork
  myfork <- parallel::mcparallel({
    eval(expr, envir = envir)
  }, silent = FALSE)
  
  # wait max n seconds for a result.
  myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout)
  # kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL)
  tools::pskill(-1 * myfork$pid, tools::SIGKILL)
  
  # clean up:
  parallel::mccollect(myfork, wait = FALSE)
  
  # timeout?
  if (is.null(myresult)) {
    if (on_timeout == "error") {
      stop("reached elapsed time limit")
    } else if (on_timeout == "warning") {
      warning("reached elapsed time limit")
    } else if (on_timeout == "silent") {
      myresult <- NA
    }
  }
  
  # move this to distinguish between timeout and NULL returns
  myresult <- myresult[[1]]
  
  if ("try-error" %in% class(myresult)) {
    stop(attr(myresult, "condition"))
  }
  
  # send the buffered response
  return(myresult)
}

run_scam_gam_ec <- function(data, field_vars){
  results <- NULL
  
  if(sum(ifelse(is.na(data$n) == FALSE, 1, 0)) > nrow(data)/2){
    results <- tryCatch(
      {
        eval_with_timeout(
          {
            formula <- paste0(
              "yield ~ s(n, bs = \"micv\") + s(s, bs = \"cv\") + I(n*ecs) + I(s*ecs) + I(s*n) + s(ecs, k = 3) + s(X, k = 5) + s(Y, k = 5) + te(X, Y, k = c(5, 5))",
              ifelse(
                length(field_vars) > 0,
                paste0(" + ", paste0(field_vars, collapse = " + ")),
                ""
              )
            ) %>% formula()
            
            scam(formula, data = data)
            
          },
          timeout = 60, # 20 seconds,
          on_timeout = "silent"
        )
      },
      error = function(cond){
        return(NULL)
      }
    )
    
    
    if (is.null(results) | is.na(results)) {
      
      formula <- paste0(
        "yield ~ s(s, k = 3) + s(n, k = 3) + s(ecs, k = 3) + I(n*ecs) + I(s*ecs) + I(s*n) + s(X, k = 4) + s(Y, k = 4) + te(X, Y, k = c(5, 5))",
        ifelse(
          length(field_vars) > 0,
          paste0(" + ", paste0(field_vars, collapse = " + ")),
          ""
        )
      ) %>% formula()
      
      results <- gam(formula, data = data)
      
    }}else{
      results <- tryCatch(
        {
          eval_with_timeout(
            {
              formula <- paste0(
                "yield ~ s(s, bs = \"cv\") + I(s*ecs) + s(ecs, k = 3) + s(X, k = 5) + s(Y, k = 5) + te(X, Y, k = c(5, 5))",
                ifelse(
                  length(field_vars) > 0,
                  paste0(" + ", paste0(field_vars, collapse = " + ")),
                  ""
                )
              ) %>% formula()
              
              scam(formula, data = data)
              
            },
            timeout = 60,
            on_timeout = "silent"
          )
        },
        error = function(cond){
          return(NULL)
        }
      )
      
      
      if (is.null(results) | is.na(results)) {
        
        formula <- paste0(
          "yield ~ s(s, k = 3) + s(ecs, k = 3) + I(s*ecs) + s(X, k = 4) + s(Y, k = 4) + te(X, Y, k = c(5, 5))",
          ifelse(
            length(field_vars) > 0,
            paste0(" + ", paste0(field_vars, collapse = " + ")),
            ""
          )
        ) %>% formula()
        
        results <- gam(formula, data = data)
        
      }
      
    }
  
  return(results)
  
}

run_scam_gam_ec_notimeout <- function(data, field_vars){
  results <- NULL

  if(sum(ifelse(is.na(data$n) == FALSE, 1, 0)) > nrow(data)/2){
    formula <- paste0(
              "yield ~ s(n, bs = \"micv\") + s(s, bs = \"cv\") + I(n*ecs) + I(s*ecs) + I(s*n) + s(ecs, k = 3) + s(X, k = 5) + s(Y, k = 5) + te(X, Y, k = c(5, 5))",
              ifelse(
                length(field_vars) > 0,
                paste0(" + ", paste0(field_vars, collapse = " + ")),
                ""
              )
            ) %>% formula()

           results <-  scam(formula, data = data)
}else{
     formula <- paste0(
                "yield ~ s(s, bs = \"cv\") + I(s*ecs) + s(ecs, k = 3) + s(X, k = 5) + s(Y, k = 5) + te(X, Y, k = c(5, 5))",
                ifelse(
                  length(field_vars) > 0,
                  paste0(" + ", paste0(field_vars, collapse = " + ")),
                  ""
                )
              ) %>% formula()

      results <- scam(formula, data = data)

    }

  return(results)

}

analysis_make_report <- function(ffy, rerun = TRUE){
  
  analysis_temp_rmd <- read_rmd(
    "Code/chapter2ec/analysis_template.Rmd"
  )
  
  analysis <- read_rmd(
    "Code/chapter2ec/1_analysis.Rmd"
  )
  
  weather <- read_rmd(here("Code/chapter2ec/0_1_weather_soil.Rmd"))
  
  analysis_rmd_y <- c(analysis_temp_rmd, weather, analysis)%>%
    gsub("field-year-here", ffy, .) %>% 
    gsub("title-here", "Analysis Report", .)
  
  
  #/*=================================================*/
  #' # Write out the rmd and render
  #/*=================================================*/
  analysis_report_rmd_file_name <- here(
    "Data", 
    ffy, 
    "analysis_report_exp.Rmd"
  )
  
  analysis_report_r_file_name <- here(
    "Data", 
    ffy, 
    "for_analysis_debug.R"
  )
  
  writeLines(analysis_rmd_y, con = analysis_report_rmd_file_name)
  
  purl(analysis_report_rmd_file_name, output = analysis_report_r_file_name)
  
}

data_summary <- function(ffy){
  data_exists <- here("Data", ffy, "analysis_data.shp") %>% 
    file.exists()
  
  data_temp <- data.table()
  
  if (data_exists){
    data_temp <- here("Data", ffy, "analysis_data.shp") %>% 
      read_sf()%>% 
      filter(!is.na(yild_vl))%>%
      data.table()%>%
      dplyr::mutate(farm = strsplit(ffy, "_")[[1]][1], 
                    field = strsplit(ffy, "_")[[1]][2], 
                    year = strsplit(ffy, "_")[[1]][3])%>%
      dplyr::select('farm', 'field', 'year', 'yild_vl', 's_rate', 'n_rate')%>%
      mutate_if(is.numeric, round, digits = 1)
  } 
  
  return(data_temp)
}

make_var_name_consistent <- function(data, dictionary) {
  col_list <- dictionary[, column]
  
  # col <- col_list[1]
  
  for (col in col_list) {
    temp_names_ls <- dictionary[column == col, names][[1]]
    
    matches <- temp_names_ls %in% names(data)
    
    if (any(matches)) { # if there is a match
      
      data <- setnames(data, temp_names_ls[matches][1], col)
    } else {
      data <- mutate(data, !!col := NA)
    }
  }
  
  return(data)
}

read_data <- function(ffy, var_ls){
  data_temp <- here("Data", ffy, "analysis_data.shp") %>% 
    read_sf() 
  dictionary <- fromJSON(
    file.path(
      here("Data"),
      "variable_name_dictionary.json"
    ),
    flatten = TRUE
  ) %>% 
    data.table()
  select_set <- as.vector(intersect(names(data_temp),var_ls))
  data_temp <- data_temp %>%         
    dplyr::select(select_set) %>%
    filter(!is.na(yild_vl)) %>%
    cbind(., st_coordinates(st_centroid(.))) %>%
    mutate(polygon_area := st_area(.)) %>%
    make_var_name_consistent(., dictionary[type == "final", ])
  return(data_temp)
}

assign_gc_rate <- function(data, input_type, gc_type, gc_rate) {
  
  data_temp <- tryCatch(
    {
      if (gc_type == "Rx") {
        #--------------------------
        # Read Rx data
        #--------------------------
        Rx <- st_read(gc_rate) %>% 
          st_set_crs(4326) %>% 
          st_transform(st_crs(data)) %>%
          # st_make_valid() %>%
          setnames(names(.), tolower(names(.)))
        
        dict_input <- dictionary[type == paste0("Rx-", tolower(input_type)), ]
        col_list <- dict_input[, column]
        
        Rx <- make_var_name_consistent(
          Rx, 
          dict_input 
        )
        
        #/*----------------------------------*/
        #' ## Unit conversion
        #/*----------------------------------*/
        if (input_type == "N") {
          Rx <- mutate(Rx, 
                       tgti = convert_N_unit(
                         input_data_n$form, 
                         input_data_n$unit, 
                         tgti, 
                         field_data$reporting_unit
                       ) 
                       # + n_base_rate # add base N rate
          )
        } else if (input_type == "S") {
          #--- seed rate conversion ---#
          if (any(Rx$tgti > 10000)){
            #--- convert to K ---#
            Rx <- mutate(Rx, tgti = tgti / 1000)
          }
        }
        
        #=== map ===#
        # tm_shape(Rx) +
        #   tm_fill(col = "tgti")
        
        #--------------------------
        # Identify grower-chosen rate by observation
        #--------------------------
        obs_tgti <- st_intersection(data, Rx) %>% 
          mutate(area = as.numeric(st_area(.))) %>% 
          data.table() %>% 
          .[, .SD[area == max(area)], by = obs_id] %>% 
          .[, num_obs_per_zone := .N, tgti] %>% 
          .[, analyze := FALSE] %>% 
          .[num_obs_per_zone >= 200, analyze := TRUE] %>% 
          .[, .(obs_id, tgti, analyze)] 
        
        data <- left_join(data, obs_tgti, by = "obs_id") %>% 
          rename(gc_rate = tgti)
      }
      
    },
    error = function(cond){
      data$gc_rate <- mean(Rx$tgti)
      return(data)
    }
  )
  
  if (gc_type == "uniform") {
    
    data$gc_rate <- gc_rate 
    
  } else {
    data <- data_temp
  }
  
  return(data)
}

get_whole_pi_test <- function(data, gam_res, crop_price, input_price) {
  
  test_data <- data.table(data) 
  
  whole_profits_test <- rbind(
    #=== opt (V) vs gc ===#
    get_dif_stat(
      test_data, 
      "input_rate", 
      "opt_input", 
      "gc_rate",
      gam_res,
      crop_price,
      input_price = input_price 
    ) %>% 
      .[, type := "VR \n vs \n GC"] %>% 
      .[, type_short := "ov vs g"],
    
    #=== site-specific vs optimal uniform without zones ===#
    get_dif_stat(
      test_data, 
      "input_rate", 
      "opt_input", 
      "opt_input_u",
      gam_res,
      crop_price,
      input_price = input_price 
    ) %>% 
      .[, type := "VR \n vs \n UR without zones"] %>% 
      .[, type_short := "ov vs ounoz"],
    
    #=== site-specific vs optimal uniform with zones ===#
    get_dif_stat(
      test_data, 
      "input_rate", 
      "opt_input", 
      "opt_input_u_zone",
      gam_res,
      crop_price,
      input_price = input_price 
    ) %>% 
      .[, type := "VR \n vs \n UR with zones"] %>% 
      .[, type_short := "ov vs ouwz"],
    
    #=== optimal uniform without zones vs optimal uniform with zones ===#
    get_dif_stat(
      test_data, 
      "input_rate", 
      "opt_input_u_zone", 
      "opt_input_u_full",
      gam_res,
      crop_price,
      input_price = input_price 
    ) %>% 
      .[, type := "UR with zones \n vs \n UR without zones"] %>% 
      .[, type_short := "ounoz vs ouwz"],
    
    #=== opt (u) with zone vs gc ===#
    get_dif_stat(
      test_data, 
      "input_rate", 
      "opt_input_u_zone", 
      "gc_rate",
      gam_res,
      crop_price,
      input_price = input_price 
    ) %>% 
      .[, type := "UR with zone \n vs \n GC"] %>% 
      .[, type_short := "ou vs g"]
  )
  
  
  return(whole_profits_test)
  
}

make_data_for_eval <- function(data, field_vars, est) {
  
  # data <- analysis_res_gcr$data[[1]]
  # est <- analysis_res_gcr$gam_res[[1]]
  
  data_dt <- data.table(data)
  
  var_names_ls <- c(field_vars, "ecs")
  
  if("zone_txt" %in% var_names_ls){
    data_for_eval <- data_dt[, ..var_names_ls] %>% 
      .[, lapply(.SD, mean), by = zone_txt]
  }else{
    data_for_eval <- data_dt[, ..var_names_ls] %>% 
      .[, lapply(.SD, mean, na.rm = TRUE)]
  }
  
  return(data_for_eval)
  
}


predict_se <- function(rate, est, profit_data){
  subset <- profit_data %>%
    .[treat == rate] %>%
    .[, yield_hat :=  as.matrix(predict(est, data = ., listw = wq, zero.policy = TRUE))[,2]] %>%
    data.frame() %>%
    dplyr::select("yield_hat") %>%
    as.vector()
  return(subset)
}

get_profit_data <- function(data_sf, est) {
  data_dt <- data.table(data_sf)
  
  if(class(est) == "gam" | class(est) == "scam"){
    if(sum(ifelse(is.na(data_dt$n) == TRUE, 1, 0)) > nrow(data_dt)/2){
      s_seq <- seq(
        quantile(data_dt$s, prob = 0.025), 
        quantile(data_dt$s, prob = 0.975), 
        by = .25
      )
      
      test_rates <- expand.grid(s = s_seq) %>%
        mutate(treat = s) %>%
        .[, "treat"] %>%
        as.vector()
      
      profit_data <- data_dt %>% 
        .[rep(1:nrow(.), length(test_rates)), ] %>% 
        .[, treat := rep(test_rates, each = nrow(.)/length(test_rates))] %>% 
        .[, s := as.numeric(treat)] %>%
        .[, yield_hat := predict(est, newdata = .)] %>% 
        .[, profit_hat := crop_price*yield_hat - s_price*s ] 
    }else{
      n_seq <- seq(
        quantile(data_dt$n, prob = 0.025), 
        quantile(data_dt$n, prob = 0.975), 
        by = 2
      )
      
      s_seq <- seq(
        quantile(data_dt$s, prob = 0.025), 
        quantile(data_dt$s, prob = 0.975), 
        by = .25
      )
      
      test_rates <- expand.grid(n = n_seq, 
                                s = s_seq) %>%
        mutate(treat = paste0(s, "_", n)) %>%
        .[, "treat"] %>%
        as.vector()
      
      profit_data <- data_dt %>% 
        .[rep(1:nrow(.), length(test_rates)), ] %>% 
        .[, treat := rep(test_rates, each = nrow(.)/length(test_rates))] %>% 
        .[, c("s", "n") := tstrsplit(treat, "_", fixed=TRUE)] %>%
        .[, s := as.numeric(s)] %>%
        .[, n :=as.numeric(n)] %>% 
        .[, yield_hat := predict(est, newdata = .)] %>% 
        .[, profit_hat := crop_price*yield_hat - n_price*n - s_price*s]
    }
  }else{
    if(sum(ifelse(is.na(data_dt$n) == TRUE, 1, 0)) > nrow(data_dt)/2){
      s_seq <- seq(
        quantile(data_dt$s, prob = 0.025), 
        quantile(data_dt$s, prob = 0.975), 
        by = .25
      )
      
      test_rates <- expand.grid(s = s_seq) %>%
        mutate(treat = s) %>%
        .[, "treat"] %>%
        as.vector()
      
      profit_data <- data_dt %>% 
        .[rep(1:nrow(.), length(test_rates)), ] %>% 
        .[, treat := rep(test_rates, each = nrow(.)/length(test_rates))] %>% 
        .[, s := as.numeric(treat)] %>%
        .[, yield_hat := est$coefficients[1] + est$coefficients[2]*s + est$coefficients[3]*s*s + est$coefficients[4]*ecs + est$coefficients[5]*ecs*s] %>%
        .[, profit_hat := crop_price*yield_hat - s_price*s]
      
      # yield_hat <- lapply(test_rates, predict_se,
      #                     est = est,
      #                     profit_data = profit_data) %>% do.call("rbind", .)
      
      # cbind(profit_data, yield_hat) %>%
      #   .[, profit_hat := crop_price*yield_hat - s_price*s]
    }else{
      n_seq <- seq(
        quantile(data_dt$n, prob = 0.025), 
        quantile(data_dt$n, prob = 0.975), 
        by = 1
      )
      
      s_seq <- seq(
        quantile(data_dt$s, prob = 0.025), 
        quantile(data_dt$s, prob = 0.975), 
        by = .25
      )
      
      test_rates <- expand.grid(n = n_seq, 
                                s = s_seq) %>%
        mutate(treat = paste0(s, "_", n)) %>%
        .[, "treat"] %>%
        as.vector()
      
      profit_data <- data_dt %>% 
        .[rep(1:nrow(.), length(test_rates)), ] %>% 
        .[, treat := rep(test_rates, each = nrow(.)/length(test_rates))] %>% 
        .[, c("s", "n") := tstrsplit(treat, "_", fixed=TRUE)] %>%
        .[, s := as.numeric(s)] %>%
        .[, n :=as.numeric(n)] %>%
        .[, yield_hat :=  est$coefficients[1] + est$coefficients[2]*s + est$coefficients[3]*n + est$coefficients[4]*s*s + est$coefficients[5]*n*n + est$coefficients[6]*ecs + est$coefficients[7]*n*ecs + est$coefficients[8]*s*ecs + est$coefficients[9]*s*n] %>%
        .[, profit_hat := crop_price*yield_hat - n_price*n - s_price*s]
      
      # yield_hat <- lapply(test_rates, predict_se,
      #                     est = est,
      #                     profit_data = profit_data) %>% do.call("rbind", .)
      # 
      # profit_data <- cbind(profit_data, yield_hat) %>%
      #   .[, profit_hat := crop_price*yield_hat - n_price*n - s_price*s]
    }
    
  }
  
  return(profit_data)
  
}

get_vr_rates <- function(profit_data){
  vr_rates <- setDT(profit_data)[ , .SD[which.max(profit_hat)], by = obs_id]
  return(vr_rates)
}

get_ur_rate <- function(profit_data){
  ur_rate <- profit_data %>%
  .[, mean_profit := mean(profit_hat, na.rm = TRUE), by = treat] %>%
  .[ , .SD[which.max(mean_profit)]]
  return(ur_rate)
}

get_profit_diff <- function(profit_data, ur_rate, vr_rates){
  ur_profit <- profit_data %>%
    .[treat == ur_rate$treat ,] %>%
    .[, profit_area := profit_hat*0.000247105*as.numeric(polygon_area)] %>%
    .[, sum(profit_area, na.rm = TRUE)]
  
  ur_acres <- profit_data %>%
    .[treat == ur_rate$treat ,] %>%
    data.frame() %>%
    mutate(non_na_area = case_when(
      is.na(profit_hat) == FALSE ~ 0.000247105*as.numeric(polygon_area),
      is.na(profit_hat) == TRUE ~ 0)) %>%
    select(non_na_area) %>%
    data.table() %>%
    .[, sum(non_na_area, na.rm = TRUE)] 
  
  vr_profit <- vr_rates %>%
    .[, profit_area := profit_hat*0.000247105*as.numeric(polygon_area)] %>%
    .[, sum(profit_area)]
  
  vr_acres <- vr_rates %>%
    .[, sum(0.000247105*as.numeric(polygon_area))]
  
  
  profit_diff <- (vr_profit/vr_acres - ur_profit/ur_acres)
  return(profit_diff)
}

predict_yield <- function(input_rate_seq, eval_data, est){
  
  cov_matrix <- as.matrix(est$resvar[-c(1,2), -c(1,2)])
  
  pred_data <- expand.grid(input_rate = input_rate_seq,
                           ecs = eval_data) %>% 
    data.frame() %>% 
    rowwise() %>%
    mutate(yield_hat = est$coefficients[1] + est$coefficients[2]*input_rate + est$coefficients[3]*input_rate^2 + est$coefficients[4]*ecs + + est$coefficients[5]*ecs*input_rate) %>%
    mutate(yield_hat_se = sqrt(t(as.matrix(c(1, input_rate, input_rate*input_rate, 1, ecs*input_rate)))%*%cov_matrix%*%(as.matrix(c(1, input_rate, input_rate*input_rate, 1, ecs*input_rate))))) %>%
    data.table()
}

convert_N_unit <- function(form, unit, rate, reporting_unit, conversion_type = "to_n_equiv") {
  conv_table <-
    fromJSON(
      here("Data","nitrogen_conversion.json"),
      flatten = TRUE
    ) %>%
    data.table() %>%
    .[, conv_factor := as.numeric(conv_factor)] %>%
    .[, form_unit := paste(type, unit, sep = "_")] %>%
    as.data.frame()
  
  if (form == "N_equiv") {
    conv_factor_n <- 1
  } else {
    conv_factor_n <- which(conv_table[, "form_unit"] %in% paste(form, unit, sep = "_")) %>%
      conv_table[., "conv_factor"]
  }
  
  if (reporting_unit == "metric") {
    conv_factor_n <- conv_factor_n * conv_unit(1, "lbs", "kg") * conv_unit(1, "hectare", "acre")
  }
  
  if (conversion_type == "to_n_equiv") {
    converted_rate <- (conv_factor_n) * rate
  } else {
    converted_rate <- (1 / conv_factor_n) * rate
  }
  
  return(as.numeric(converted_rate))
}

get_dif_stat <- function(data, test_var, opt_var, gc_var, gam_res, crop_price, input_price){
  
  if ("scam" %in% class(gam_res)) {
    gam_coef <- gam_res$coefficients.t
    gam_V <- gam_res$Ve.t
  } else {
    gam_coef <- gam_res$coefficients
    gam_V <- gam_res$Ve
  }
  
  base_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(gc_var)]
  
  comp_data <- data.table::copy(data) %>% 
    .[, (test_var) := get(opt_var)]
  
  #/*----------------------------------*/
  #' ## Profit (gc)
  #/*----------------------------------*/
  Xmat_base <- predict(gam_res, newdata = base_data, type = "lpmatrix") 
  # predict(gam_res, newdata = base_data) %>% mean
  
  #--- vector of 1s for summation divided by the number of observations for averaging ---#
  ones <- matrix(1 / dim(Xmat_base)[1], 1, dim(Xmat_base)[1])
  
  #--- average yield ---#
  yhat_base <- ones %*% Xmat_base %*% gam_coef
  
  #--- point estimate of profit differential ---#
  pi_gc <- crop_price * yhat_base - (input_price * ones %*% base_data$input_rate)  
  
  big_mat_base <- ones %*% Xmat_base
  
  #--- se of the profit differential  ---# 
  pi_gc_se <- crop_price * sqrt(big_mat_base %*% gam_V %*% t(big_mat_base))
  
  #/*----------------------------------*/
  #' ## Profit (optimal)
  #/*----------------------------------*/
  Xmat_comp <- predict(gam_res, newdata = comp_data, type = "lpmatrix") 
  
  #--- vector of 1s for summation divided by the number of observations for averaging ---#
  ones <- matrix(1 / dim(Xmat_comp)[1], 1, dim(Xmat_comp)[1])
  
  #--- average yield ---#
  yhat_comp <- ones %*% Xmat_comp %*% gam_coef
  
  #--- point estimate of profit differential ---#
  pi_opt <- crop_price * yhat_comp - (input_price * ones %*% comp_data$input_rate)  
  
  big_mat_comp <- ones %*% Xmat_comp
  
  #--- se of the profit differential  ---# 
  pi_opt_se <- crop_price * sqrt(big_mat_comp %*% gam_V %*% t(big_mat_comp))
  
  #/*----------------------------------*/
  #' ## Profit differential
  #/*----------------------------------*/
  #--- difference in X mat ---#
  X_dif_mat <- Xmat_comp - Xmat_base
  
  #--- vector of 1s for summation divided by the number of observations for averaging ---#
  ones <- matrix(1 / dim(X_dif_mat)[1], 1, dim(X_dif_mat)[1])
  
  #--- X_dif_mat summed ---#
  big_mat_dif <- ones %*% X_dif_mat
  
  #--- point estimate of profit differential ---#
  pi_dif <- ones %*% ((crop_price * X_dif_mat %*% gam_coef) - input_price * (comp_data$input_rate - base_data$input_rate))  
  
  #--- se of the profit differential  ---# 
  pi_dif_se <- crop_price * sqrt(big_mat_dif %*% gam_V %*% t(big_mat_dif))
  
  #--- t-stat ---#
  t_stat <- (pi_dif/pi_dif_se) %>% round(digits = 2)  
  
  return_data <- data.table(
    yhat_est_gc = yhat_base[1, 1],
    point_est_gc = pi_gc[1, 1],
    point_est_gc_se = pi_gc_se[1, 1],
    yhat_est_opt = yhat_comp[1, 1],
    point_est_opt = pi_opt[1, 1],
    point_est_opt_se = pi_opt_se[1, 1],
    point_est_dif = pi_dif[1, 1],
    point_est_dif_se = pi_dif_se[1, 1],
    t = t_stat[1, 1]
  )
  
  return(return_data)
  
}

read_rmd <- function(file_name) {
  file_name <-
    here(file_name)
  
  rmd_file <- readLines(file_name)
  
  return(rmd_file)
}

find_field_vars <- function(data_sf) {
  
  #/*----------------------------------*/
  #' ## pick field vars
  #/*----------------------------------*/
  #=== keep only the ones that are available ===#
  field_var_ls <- c(
    #=== topography ===#
    # "twi", 
    "tpi", "elevation", "slope", 
    #=== ssurgo ===#
    "clay", "sand", "water_storage",
    #=== ec ===#
    "om"
  ) %>% 
    .[. %in% names(data_sf)]
  
  #=== find variables to keep ===#
  keep_vars <- data_sf[, field_var_ls] %>% 
    st_drop_geometry() %>% 
    #=== if missing for more than 10%, then drop ===#
    data.table() %>% 
    .[, lapply(
      .SD, 
      function(x) {
        (sum(is.na(x))/nrow(data_sf)) < 1
      }
    )] %>% 
    as.matrix() %>% 
    as.vector()
  
  field_var_ls <- field_var_ls[keep_vars]
  
  return(field_var_ls)
}

get_county <- function(ffy){
  counties <- st_read(here("Data/US_County_Boundaries.shp")) %>%
    st_make_valid()
  data <- st_read(here("Data", ffy, "analysis_data.shp")) %>%
    st_transform(4326) %>%
    st_make_valid()
  
  county <- st_intersection(data[1,], counties) %>%
    as.data.frame() %>%
    dplyr::select("COUNTY", "STATE")
    
  
  return(county)
}

