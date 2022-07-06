
#============================================================
# SystemicDCCroll: Function to compute VaR and CoVaR forecasts for rolling DCC-GARCH models
# ===========================================================

# Inputs:
# data         = a tsibble containing the financial data
# DCC_spec     = a DCC model specification from the rmgarch package
# length_IS    = the length of the in-sample estimation length of the rolling windows in days
# refit_freq   = the model refit frequency in days
# alpha        = (scalar) risk level for CoVaR
# beta         = (possibly multiple) risk levels for VaR and MES
# H_sqrt       = the method to compute the decomposition of the covariance matrix H. It has to be either "Cholesky" or "symmetric"
# DCC_FC_dist  = the method to compute the VaR and CoVaR forecasts from the model innovations.
#                  It has to be either "std" for VaR and CoVaR obtained from a multivariate Student t-distribution or
#                  "np" for nonparametrically estimated VaR and CoVaR from the empirical residual distribution

SystemicDCCroll <- function(data, DCC_spec,
                            length_IS=1000, refit_freq=100,
                            beta, alpha, H_sqrt="Cholesky", DCC_FC_dist=c("std","np")){

  if (H_sqrt=="symmetric" & any(DCC_FC_dist=="std")){
    warning("The symmetric square root decomposition does not (yet) work with the forecasts based on Student-t distributed residuals. Hence, the nonparametric technique is used.")
    DCC_FC_dist <- "np"
  }

  # Check inputs:
  if (!is_tsibble(data)) stop("Error: Please enter a 'tsibble' object for the argument 'data'.")


  TT <- dim(data)[1]
  length_OOS <- TT - length_IS

  # The refit_points correspond to the first FC target date. Maybe not optimal.
  refit_points <- seq(length_IS+1, TT, by=refit_freq)

  H_FC <- array(NA, dim=c(2,2,length_OOS))
  resid_DCC_refit <- array(NA, dim=c(length_IS,2,length(refit_points)))
  shape_param_DCC_refit <- array(NA, dim=c(length(refit_points)))
  eps_DCC_refit <- array(NA, dim=c(length_IS,2,length(refit_points)))
  data_IS_list <- list()

  tt_counter <- 1
  for (tt in refit_points){
    # tt is the first OOS time point in each refit window!
    start_dfIS <- tt-length_IS
    end_dfIS <- min(tt+refit_freq-2,TT-1) # Use the min as the last equally spaced window might end after TT
    df_IS <- data[start_dfIS:end_dfIS,] %>%
      tibble::as_tibble() %>%
      dplyr::select(x,y)
    n_roll_tt <- min(refit_freq-1, TT-tt) # Only roll for the remaining steps after the last refit, might be less than refit_freq-1

    # DCC fit
    DCC_fit = dccfit(DCC_spec, data = df_IS, out.sample = n_roll_tt,  # out.sample: the hold-out observations for dccforecast to work
                     fit.control = list(eval.se = FALSE, stationarity = TRUE, scale = FALSE))
    H <- rcov(DCC_fit)

    # DCC forecasts
    DCC_FC <- dccforecast(DCC_fit, n.ahead = 1, n.roll = n_roll_tt)
    H_FC_list <- rcov(DCC_FC)
    H_FC[,,(tt-length_IS):(tt-length_IS+n_roll_tt)] <- do.call(abind::abind, c(H_FC_list, along = 3)) # Convert list to array
    # rho_FC <- H_FC[1,2,]

    # Either symmetric or Cholesky decompositions for both, the fitted and forecasted covariance matrices H and H_FC
    if (H_sqrt=="symmetric"){
      Sigma <- sapply(1:length_IS, function(i) pracma::sqrtm(H[,,i])$B, simplify = "array")
    } else if (H_sqrt=="Cholesky"){
      Sigma <- sapply(1:length_IS, function(i) t(chol(H[,,i])), simplify = "array") # base::chol() returns the upper triangular matrix; needs to be transposed!
    } else {
      warning("Please specify a valid method for 'H_sqrt'. The symmetric decomposition is used instead.")
      Sigma <- sapply(1:length_IS, function(i) pracma::sqrtm(H[,,i])$B, simplify = "array")
    }

    # Standardized residuals
    Sigma_inv <- sapply(1:length_IS, function(i) solve(Sigma[,,i]), simplify = "array")
    eps_DCC_refit[,,tt_counter] <- sapply(1:length_IS, function(i) Sigma_inv[,,i] %*% t(as.matrix(df_IS))%>%.[,i])
    shape_param_DCC_refit[tt_counter] <- rshape(DCC_fit) # scalar
    # resid_DCC_refit[,,tt_counter] <- residuals(DCC_fit) # length_IS x 2

    # Collect x,y, and the corresponding standardized residuals eps in a list of data.frames;
    # This is used later on for the nonparametric version of DCCtoCoVaR
    data_IS_list <- rlist::list.append(data_IS_list,
                                       df_IS[start_dfIS:(end_dfIS-n_roll_tt),] %>%
                                         dplyr::mutate(eps_x = eps_DCC_refit[,1,tt_counter],
                                                       eps_y = eps_DCC_refit[,2,tt_counter]))

    tt_counter <- tt_counter + 1
  }

  # Either symmetric or Cholesky decompositions for H_FC
  if (H_sqrt=="symmetric"){
    Sigma_FC <- sapply(1:dim(H_FC)[3], function(i) pracma::sqrtm(H_FC[,,i])$B, simplify = "array")
  } else if (H_sqrt=="Cholesky"){
    Sigma_FC <- sapply(1:dim(H_FC)[3], function(i) t(chol(H_FC[,,i])), simplify = "array") # base::chol() returns the upper triangular matrix; needs to be transposed!
  } else {
    # Note: A warning is already given above!
    Sigma_FC <- sapply(1:dim(H_FC)[3], function(i) pracma::sqrtm(H_FC[,,i])$B, simplify = "array")
  }


  # Transform DCC Covariance forecasts to system risk measure forecasts
  FC_df <- tibble()

  # Student t-distribution
  if (any(DCC_FC_dist=="std")){
    SRM_FCs_tmp <- sapply(1:(TT-length_IS), function(i){
      DCC_to_CoVaR(Sigma=Sigma_FC[,,i],
                   nu=shape_param_DCC_refit[ceiling(i/refit_freq)] %>% round(), #See that the correct nu is called here!
                   data_IS=NULL,
                   alpha=alpha,
                   beta=beta) %>% as.numeric()})

    FC_df <- bind_rows(FC_df,
                       tibble::tibble(Date=data$Date[(length_IS+1):TT],
                                      VaR=SRM_FCs_tmp[1,],
                                      CoVaR=SRM_FCs_tmp[3,],
                                      MES=SRM_FCs_tmp[2,],
                                      model="DCC-GARCH-std",
                                      row.names = NULL))

  }


  # Nonparametric estimates
  if (any(DCC_FC_dist=="np")) {
    SRM_FCs_tmp <- sapply(1:(TT-length_IS), function(i){
      DCC_to_CoVaR(Sigma=Sigma_FC[,,i],
                   nu=NULL,
                   data_IS=data_IS_list[[ceiling(i/refit_freq)]],
                   alpha=alpha,
                   beta=beta) %>% as.numeric()})

    FC_df <- dplyr::bind_rows(FC_df,
                              tibble::tibble(Date=data$Date[(length_IS+1):TT],
                                             VaR=SRM_FCs_tmp[1,],
                                             CoVaR=SRM_FCs_tmp[3,],
                                             MES=SRM_FCs_tmp[2,],
                                             model="DCC-GARCH-np",
                                             row.names = NULL))
  }



  FC_df <- FC_df %>%
    dplyr::left_join(tibble::as_tibble(data), by="Date") %>%
    dplyr::select(Date, x, y, VaR, CoVaR, MES, model)

  obj <- list(SRM_FC=FC_df,
              Sigma_FC=Sigma_FC,
              refit_points=refit_points,
              shape_param=shape_param_DCC_refit,
              beta=beta,
              alpha=alpha)
  class(obj) <- "SystemicDCCroll"

  return(obj)
}



print.SystemicDCCroll <- function(obj){
  print(obj$SRM_FC)
}


autoplot.SystemicDCCroll <- function(obj, SRM_plot="CoVaR", model_plot="DCC-GARCH-np", facet_names=c("X / VaR","Y / CoVaR")){

  df_tmp <- obj$SRM_FC %>%
    dplyr::filter(model==model_plot) %>%
    dplyr::mutate(VaR_violation=(x > VaR))

  # Only plot the CoVaR for now!
  df_long <- left_join(
    df_tmp %>%
      dplyr::select(Date, VaR_violation, x, y) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="NegativeReturns"),
    df_tmp %>%
      dplyr::select(Date, VaR_violation, VaR, SRM_plot) %>%
      dplyr::rename(x=VaR, y=SRM_plot) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="SRMForecasts"),
    by=c("Date", "VaR_violation", "Symbol")) %>%
    as_tibble()

  levels(df_long$Symbol) <- facet_names

  p <- ggplot2::ggplot(df_long %>% arrange(VaR_violation)) +
    ggplot2::geom_point(aes(x=Date, y=NegativeReturns, color=VaR_violation)) +
    ggplot2::scale_colour_manual(values = c("grey", "black")) +
    ggnewscale::new_scale_color() + # For two color scales!!!
    ggplot2::geom_line(aes(x=Date, y=SRMForecasts, color=Symbol)) +
    ggplot2::scale_colour_manual(values = c("red", "blue")) +
    ggplot2::facet_wrap(~Symbol, ncol=1) +
    ggplot2::theme_bw()

  p
}


plot.SystemicDCCroll <- function(obj, ...){
  p <- autoplot(obj, ...)
  print(p)
}





#============================================================
# DCC_to_CoVaR: Function to compute non-parametric VaR, MES and CoVaR forecasts for DCC-GARCH models
# ===========================================================

# Inputs:
# Sigma_FC     = forecasted (2 x 2) conditional 'square-root' decomposition of covariance matrix (for time n+1)
# nu           = scalar degrees of freedom of multivariate t-distribution (= NULL, if data_IS is supplied)
# data_IS      = a data.frame containing the in-sample values of x,y, and the standardized in-sample residuals (= NULL, if nu is supplied)
# alpha        = risk/probability level for CoVaR
# beta         = risk/probability levels for VaR and MES

DCC_to_CoVaR <- function(Sigma_FC, nu, data_IS, alpha, beta){

  if( is.null(data_IS) ){
    ### (A) Forecasts for t-distributed innovations
    # 1. Calculate VaR, ES and MES of eps
    VaR.hat.b <- MES.hat.b <- ES.hat.b <- numeric(length(beta))   # Note that CoVaR depends on $beta$ through the beta-quantile of eps_X
    sd.wt <- sqrt(nu / (nu-2) )
    ScaleMatrix_tdist <- (nu - 2) / nu * matrix(c(1, 0, 0, 1), nrow=2) # Correct for scale != Variance
    f.MES <- function(x, nu) { x[1] * dmvt(x, sigma = ScaleMatrix_tdist, df = nu, log = FALSE) } # "x" is vector

    for(b in 1:length(beta)){
      VaR.hat.b[b] <- qt(p=beta[b], df = nu) / sd.wt
      MES.hat.b[b] <- 1/(1-beta[b]) * adaptIntegrate(f.MES, lowerLimit = c(-Inf, VaR.hat.b[b]), upperLimit = c(Inf, Inf), nu)$integral
    }
    ES.hat.b <- 1/sd.wt * (dt(qt(beta, df=nu) , df=nu) / (1 - beta)) * ((nu + (qt(beta, df=nu))^2) / (nu - 1))  # see [MFE15, Example 2.15]

    # 2. Compute conditional VaR, CoVaR and MES
    # The following formula ONLY apply to the Cholesky decomposition! (This could be implemented easier given we have Sigma_FC already, right?)
    H_FC <- Sigma_FC %*% t(Sigma_FC) # Covariance forecast
    rho.np1         <- H_FC[2,1] / sqrt( H_FC[1,1] * H_FC[2,2] )
    theta.hat.n.b.1 <- sqrt(H_FC[2,2]) * sqrt(1 - rho.np1^2) * MES.hat.b      # MES part
    theta.hat.n.b.2 <- sqrt(H_FC[2,2]) * rho.np1 *  ES.hat.b                  # ES part
    MES_FC     <- theta.hat.n.b.1 + theta.hat.n.b.2                      # this is the MES forecast
    VaR_FC     <- sqrt(H_FC[1,1]) * VaR.hat.b                            # this is the VaR forecast
    CoVaR_FC   <- numeric( length(beta))
    for(b in 1:length(beta)){
      CoVaR_FC[b] <- CoVaR.true.t(alpha, beta[b], nu, H_FC)
    }

  } else {
    ### (B) Forecasts based on nonparametrically computing the (systemic) risk measures from past DCC-innovations
    # 1. Compute VaR forecasts
    n <- nrow(data_IS)

    # u_hat := Sigma_Forecast * eps_InSample, a combination of the in sample residuals and forecastsed covariance matrix
    u_hat <- Sigma_FC %*% (data_IS %>% dplyr::select(eps_x, eps_y) %>% as.matrix() %>% t())
    data_IS$u_hat_x <- u_hat[1,]
    data_IS$u_hat_y <- u_hat[2,]

    VaR_FC <- quantile(data_IS$u_hat_x, probs = beta)

    # 2. Compute conditional risk measures
    MES_FC <- CoVaR_FC <- numeric( length(beta) ) # Note that CoVaR depends on $beta$ through the beta-quantile of eps_X

    for(b in 1:length(beta)){
      u_hat_Y_VaRviolation <- data_IS %>%
        dplyr::filter(u_hat_x > VaR_FC[b]) %>%
        dplyr::pull(u_hat_y)
      CoVaR_FC[b]  <- quantile(u_hat_Y_VaRviolation, probs=alpha)
      MES_FC[b]    <-     mean(u_hat_Y_VaRviolation)
    }
  }

  return(list(VaR = VaR_FC, MES = MES_FC, CoVaR = CoVaR_FC))
}





# ## Test
# data_sim <- tsibble(Date=1:5000,
#                     x=rnorm(5000,0,1),
#                     y=rnorm(5000,0,1.5),
#                     index=Date)
#
# garch11.spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean=FALSE),
#                           variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
#                           distribution.model = "std")
# dcc.garch11.spec = dccspec(uspec = multispec( replicate(2, garch11.spec) ),   # use standard univariate GARCH models for marginals
#                            dccOrder = c(1,1),
#                            distribution = "mvt")
#
# # 3. DCC GARCH forecasts
# asdf <- SystemicDCCroll(data=data_sim,
#                         DCC_spec=dcc.garch11.spec,
#                         length_IS=3000, refit_freq=1000,
#                         DCC_FC_dist=c("std","np"), beta=0.95, alpha=0.95)
#



