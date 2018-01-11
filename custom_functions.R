
#### Then we determine the small area of the policy holder
sample_cluster_data <- function(iter_id, dt_smallarea, dt_clusterprod, cluster_id) {

    ### First we sample from the appropriate small areas to get some data
    N <- length(cluster_id[cluster_id == iter_id])
    use_dt <- dt_smallarea[cluster_id %in% iter_id]

    ### We sample according to population
    pop_cumsum <- cumsum(use_dt$total2011)

    sample_val <- sample(1:pop_cumsum[nrow(use_dt)], N, replace = TRUE)

    idx <- sapply(1:N, function(iter) {
            length(pop_cumsum[sample_val[iter] <= pop_cumsum])
    })

    sa_id <- use_dt[idx]$sa_id


    ### Now we use the product data to calculate the policy type
    use_dt <- dt_clusterprod[cluster_id %in% iter_id]
    prod_type_values <- use_dt$variable
    prod_type_prop   <- use_dt$value

    prod_type <- sample(prod_type_values, N, prob = prod_type_prop, replace = TRUE)

    dt_clusterdata <- data.table(cluster_id = iter_id
                                 ,sa_id     = sa_id
                                 ,prod_type = prod_type
                                 ,idx       = idx)

    return(dt_clusterdata)
}


### This function generates policy data based on the cluster type
sample_product_data <- function(iter_type, dt_productdata, premape_shape, premape_scale) {
    use_dt <- dt_productdata[prod_type == iter_type]

    N <- dt_policy[prod_type == iter_type, .N]

    policy_start_age <- round(rnorm(N, use_dt$age_mean, use_dt$age_sd), 0)
    policy_start_age <- pmax(use_dt$age_min, policy_start_age)
    policy_start_age <- pmin(use_dt$age_max, policy_start_age)

    prem_ape <- rgamma(N, shape = premape_shape, scale = premape_scale) * use_dt$prem_multiplier
    prem_ape <- round(prem_ape, 2)

    prem_ape <- pmax(use_dt$prem_min, prem_ape)

    rpsp_prop <- use_dt[, c(1 - prem_sp_prop, prem_sp_prop)]
    prem_type <- sample(c("RP", "SP"), N, prob = rpsp_prop, replace = TRUE)

    prod_data_dt <- data.table(prod_type       = iter_type
                              ,prem_type       = prem_type
                              ,policy_startage = policy_start_age
                              ,prem_ape        = prem_ape
                               )
    return(prod_data_dt)
}


### This function calculates the premium for the protection policy
calculate_protection_premium <- function(applicant_age, duration, isjointlife, islifeonly, acttables, interest_rate) {

    if(isjointlife) {
        life_factor <- Axyzn(acttables[c('m_lifetable', 'f_lifetable')]
                             ,x = rep(applicant_age, 2)
                             ,n = duration
                             ,i = interest_rate) /
                       axyzn(acttables[c('m_lifetable', 'f_lifetable')]
                             ,x = rep(applicant_age, 2)
                             ,n = duration
                             ,i = interest_rate)

        morb_factor <- Axyzn(acttables[c('mn_morb_lifetable', 'fn_morb_lifetable')]
                             ,x = rep(applicant_age, 2)
                             ,n = duration
                             ,i = interest_rate) /
                       axyzn(acttables[c('mn_morb_lifetable', 'fn_morb_lifetable')]
                             ,x = rep(applicant_age, 2)
                             ,n = duration
                             ,i = interest_rate)
    } else {
        life_factor <- Axn(acttables[['m_lifetable']]
                           ,x = applicant_age
                           ,n = duration
                           ,i = interest_rate) /
                       axn(acttables[['m_lifetable']]
                           ,x = applicant_age
                           ,n = duration
                           ,i = interest_rate)

         morb_factor <- Axn(acttables[['mn_morb_lifetable']]
                            ,x = applicant_age
                            ,n = duration
                            ,i = interest_rate) /
                        axn(acttables[['mn_morb_lifetable']]
                            ,x = applicant_age
                            ,n = duration
                            ,i = interest_rate)
    }


    prem <- ifelse(islifeonly, life_factor, life_factor + morb_factor)

    return(prem)
}


### Calculate the mortality rating based on some parameters
calculate_mort_rating <- function(N, shape, scale) {
    rating <- round(rgamma(N, shape, scale), 0)

    rating <- pmax(rating,  4)
    rating <- pmin(rating, 20)

    return(rating * 25)
}


### This function creates a function that generates a lapse month
create_lapse_calculation_function <- function(risk_baseline
                                             ,risk_gender
                                             ,risk_smoker
                                             ,risk_mortstat
                                             ,risk_logpremape
                                             ,risk_highpremape
                                             ,risk_cluster
                                             ,dt_calendar            = NULL
                                             ,min_monthly_lapse_prob = min_monthly_lapse_prob
                                             ,premape_high           = premape_high
                                               ) {

    calculate_lapse_time <- function(dt_data, verbose = FALSE) {
        stopifnot(nrow(dt_data) == 1)

        ### Create the monthly hazard probabilities
        lapse_prob <- risk_baseline +
                      risk_gender * (dt_data$gender_life1 == 'M') +
                      risk_smoker  [[dt_data$smoker_life1]] +
                      risk_mortstat[[dt_data$mortgage_status]] +
                      risk_logpremape * log(dt_data$prem_ape) +
                      risk_highpremape * (dt_data$prem_ape > premape_high) +
                      risk_cluster [[dt_data$cluster_id]]

        max_month  <- length(lapse_prob)

        ### If we have included calendar effects we need to calculate and include them
        if(!is.null(dt_calendar)) {
            month_date <- seq(dt_data$policy_startdate, length.out = max_month, by = 'month')

            dt_month <- data.table(calc_date = month_date)
            setkey(dt_month, calc_date)

            calendar_prob <- dt_calendar[dt_month, roll = TRUE]$lapse_prob
        } else {
            calendar_prob <- rep(0, max_month)
        }

        lapse_prob <- lapse_prob + calendar_prob


        lapse_prob <- pmax(min_monthly_lapse_prob, lapse_prob)

        month_idx <- which(runif(max_month) < lapse_prob)

        if(length(month_idx) == 0) {
            lapse_month <- max_month
        } else {
            lapse_month <- month_idx[1]
        }

        if(verbose) print(lapse_prob[1:60])

        return(lapse_month)
    }

    return(calculate_lapse_time)
}


generate_random_startdates <- function(N
                                      ,data_startdate
                                      ,data_snapshotdate
                                      ,datestr_exclude
                                      ,year_weight_dt
                                      ,month_weight_dt
                                      ,dow_weight_dt) {

    print(data_startdate); print(data_snapshotdate);
    startdate_universe <- seq(data_startdate, data_snapshotdate, by = 'day')
    startdate_universe <- startdate_universe[!format(startdate_universe, "%m%d") %in% datestr_exclude]

    year_prob_dt <- data.table(year      = year_weight_dt$year
                              ,yearstr   = year_weight_dt$year
                              ,year_prop = rdirichlet(1, 20 * year_weight_dt$weight)[1,]
                               )

    n_years <- year_prob_dt %>% nrow

    month_prob_dt <- rdirichlet(n_years, 10 * month_weight_dt$weight) %>%
        melt %>%
        mutate(year       = year_prob_dt$year[Var1]
              ,month      = month_weight_dt$month[Var2]
              ,monthstr   = sprintf('%02d', Var2)
              ,month_prop = value
               )
    setDT(month_prob_dt)

    yearmonth_props_dt <- year_prob_dt %>%
        left_join(month_prob_dt, by = 'year') %>%
        mutate(yearmon = paste0(yearstr, '-', monthstr)
              ,prob    = year_prop * month_prop)

    startdate_ym <- sample(yearmonth_props_dt$yearmon, N, replace = TRUE, prob = yearmonth_props_dt$prob)

    ym_count_dt <- data.table(startdate_ym = startdate_ym) %>%
        group_by(startdate_ym) %>%
        summarise(count = n()) %>%
        arrange(startdate_ym) %>%
        ungroup

    startdates_dt <- map2(ym_count_dt$startdate_ym
                         ,ym_count_dt$count
                         ,function(ym, count) generate_random_month_dates(ym, count, dow_dt)) %>%
        rbindlist %>%
        arrange(dates)


    return(startdates_dt)
}


generate_random_month_dates <- function(ym_datestr, N, dow_prob_dt) {

    cat(paste0("Generating start dates for ", ym_datestr, "\n"))
    start_datestr <- paste0(ym_datestr, "-01")

    month_count <- Hmisc::monthDays(start_datestr)

    valid_dates <- seq(as.Date(start_datestr), as.Date(paste0(ym_datestr, "-", month_count)), by = 'day')

    dates_dt <- data.table(date    = valid_dates
                          ,weeknum = format(valid_dates, "%u")
                          ,dow     = format(valid_dates, "%a")
                           ) %>%
        filter(weeknum %in% 1:5
              ,!format(date, "%m%d") %in% holiday_datestr) %>%
        inner_join(dow_prob_dt, by = 'dow') %>%
        mutate(prob = weight / sum(weight)) %>%
        arrange(date)

    dates <- sample(dates_dt$date, N, replace = TRUE, prob = dates_dt$prob)

    startdates_dt <- data.table(dates = sort(dates))

    return(startdates_dt)
}
