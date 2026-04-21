### After transformation
applytrans = function(value, num){
  eps = 1e-05;
  if(num == 1){ fx = value;  }
  if(num == 2){ fx = diff(value);}
  if(num == 3){ fx = diff(diff(value));}
  if(num == 4){ fx = log(value + eps);}
  if(num == 5){ fx = diff(log(value+ eps));}
  if(num == 6){ fx = diff(diff(log(value+ eps)));}
  if(num == 7){ 
    TT = length(value)
    y = (imputedata[-1] / imputedata[-TT] - 1)
    fx = diff(y);
  }
  return(fx)
}



## Estimate Factors
getFactor = function(Y, r, isCorr){
  ## Input data is dim*length
  if(missing(isCorr)){isCorr = TRUE}
  if(!isCorr){ 
    ssd = apply(Y, 1, sd);
    Y = Y/ssd;
  }
  vv = svd(cov(t(Y))); 
  ##  vv = svd(cor(t(Y))); 
  q=dim(Y)[1];
  Lam = vv$u[,1:r]*sqrt(q);
  Fh = 1/q*t(Lam)%*%Y;
  return(Fh)
}



ICselect = function(Y, maxr, takeLog, sigmaE){
  q =dim(Y)[1]; n =ncol(Y);
  if(missing(sigmaE)){ sigmaE = FALSE};
  if(missing(takeLog)){ takeLog = FALSE};
  if(missing(maxr)){ maxr = q};
  
  
  vv = svd(cov(t(Y))); 
  if(sigmaE == TRUE){
    Lam = vv$u[,1:maxr]*sqrt(q);
    Fh = 1/q*t(Lam)%*%Y;
    Yhat = Lam%*%Fh;
    sigmamax = sum((Y-Yhat)^2)/(q*n);
  } else{
    sigmamax=1;
  }
  
  psum1 =psum2 = psum3 = psum4 = err= numeric(maxr);
  for(r in 1:maxr){
    Lam = vv$u[,1:r]*sqrt(q);
    Fh = 1/q*t(Lam)%*%Y;
    Yhat = Lam%*%Fh;
    if(takeLog == TRUE){
      err[r]= log(sum((Y-Yhat)^2)/(q*n));
    } else{
      err[r]= sum((Y-Yhat)^2)/(q*n);
    }
    
    psum1[r] = err[r] + r*sigmamax*(q+n)/(q*n)*log(n*q/(n+q));
    psum2[r] = err[r] + r*sigmamax*(q+n)/(q*n)*log(min(n,q));
    psum3[r] = err[r] + r*sigmamax*(log(min(n,q))/min(n,q));
    psum4[r] = err[r] + r*sigmamax*(n+q-r)*(log(n*q))/(n*q);
  }
  id1 = which.min(psum1);
  id2 = which.min(psum2);
  id3 = which.min(psum3);
  id4 = which.min(psum4);
  out = list();
  out$err = err;
  out$psum1 = psum1;
  out$psum2 = psum2;
  out$psum3 = psum3;
  out$psum4 = psum4;
  
  out$order = c(id1, id2, id3, id4);
  return(out)
}


# x: 숫자 벡터(예: 분기값이 3회씩 복제된 월별 데이터)
# method_last: 마지막 구간 처리("extrapolate"=이전 기울기로 외삽, "carry"=마지막 값 유지)
interpolate_duplicated <- function(x, method_last = c("extrapolate", "carry")) {
  stopifnot(is.numeric(x))
  method_last <- match.arg(method_last)
  
  # 동일값이 연속된 구간(run) 단위로 요약
  r <- rle(x)
  qval <- r$values         # 각 구간의 대표값(분기값)
  qlen <- r$lengths        # 각 구간 길이(예: 3개월)
  
  nseg <- length(qval)
  if (nseg == 0L) return(x)
  if (nseg == 1L) {
    # 분기값이 하나뿐이면: 외삽은 0 기울기, carry도 동일
    return(rep(qval, qlen))
  }
  
  # 인접 분기값 차이 및 구간별 월간(혹은 반복단위) 기울기
  dq <- diff(qval)                           # 길이 nseg-1
  slope <- dq / qlen[-length(qlen)]          # 각 구간 내 한 스텝 증가량
  
  # 결과 벡터 구성
  out <- numeric(length(x))
  idx_start <- 1L
  
  for (i in seq_len(nseg)) {
    len_i <- qlen[i]
    if (i < nseg) {
      step <- slope[i]
      # 다음 분기로 향해 len_i-1번만 증가(마지막 값은 다음 분기의 첫 값이 됨)
      seg <- qval[i] + step * (0:(len_i - 1))
    } else {
      # 마지막 구간 처리
      if (method_last == "extrapolate") {
        step_last <- slope[length(slope)]       # 직전 구간의 기울기 재사용
        seg <- qval[i] + step_last * (0:(len_i - 1))
      } else { # "carry"
        seg <- rep(qval[i], len_i)
      }
    }
    out[idx_start:(idx_start + len_i - 1)] <- seg
    idx_start <- idx_start + len_i
  }
  out
}



get_mode <- function(x) {
  tab <- table(x)
  mode_val <- names(tab)[which.max(tab)]
  return(as.numeric(mode_val))
}



FAVAR_IRF = function(XX, YY, r=4, bp=.25, impact=NULL, nahead=48, ncores=25){
  ## Four factors
  out = getFactor(Y=XX, r=r);
  data_var = data.frame(t(out), YY);
  ## The last variable is the impulse; (usually FYFF)
  if(is.null(impact)){
  ylast = names(data_var)[ncol(data_var)]
  } else{
  ylast = impact;
  }
  # Shock size of 25 basis points
  var = vars::VAR(data_var, p = 13);
  irf_point = vars::irf(var, n.ahead = nahead, impulse = ylast, response = ylast, boot = FALSE)
  impulse_sd = bp/sd(data_var[,ylast])
  scale = impulse_sd/(irf_point$irf[[1]][1]) # position of FYFF response at step 0
  
  # Computing Loading Factors
  ZZ = data.frame(t(XX), YY); ZZ = as.matrix(ZZ);
  reg_loadings = lm(ZZ ~ 0 + ., data=data_var)
  loadings = reg_loadings$coefficients
  
  #### BOOTSTRAPING ########
  # --- User-defined inputs ---
  # var      : fitted VAR model using lineVar
  # loadings : (K x nvars) matrix
  # scale    : scalar multiplier
  # nahead   : forecast horizon
  # nvars    : number of variables
  # R        : bootstrap replications
  # nsteps   : same as nahead
  
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
  set.seed(123) 
  nrep=1000;
  var = lineVar(data_var, lag = 13, include = "const")
  # --- Parallel bootstrap: returns list of (nahead × nvars) matrices ---
  IRF_list <- foreach(j = 1:nrep, .packages = c("vars", "tsDyn"), .inorder = TRUE) %dopar% {
    data_boot <- VAR.boot(var, boot.scheme = "resample")
    var_boot  <- vars::VAR(data_boot, p = 13)
    irf1      <- vars::irf(var_boot, n.ahead = nahead, impulse = ylast, boot = FALSE)
    (irf1$irf[[1]] %*% loadings) * scale
  }
  stopCluster(cl)
  
  IRF_array <- simplify2array(IRF_list)  # shape: [nsteps, nvars, R]
  # --- Vectorized quantile computation ---
  qfun <- function(x, p) apply(x, c(1, 2), quantile, probs = p)
  Upper <- qfun(IRF_array, 0.9)
  Lower <- qfun(IRF_array, 0.1)
  IRF   <- qfun(IRF_array, 0.5)
  # Clean up memory
  rm(IRF_list, IRF_array)
  
  return(list(Upper = Upper, Lower = Lower, IRF = IRF))
}

plotIRFbands <- function(IRFobj, tcode, variables,
                         IRFobj2 = NULL,
                         IRFobj3 = NULL,
                         model1_name = "Model 1",
                         model2_name = "Model 2",
                         model3_name = "Model 3",
                         save_png = FALSE, ncol = 3,
                         file = NULL, ct = 1.5,
                         width  = 11, height = 8,
                         units  = "in", res = 300) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2'가 필요합니다. install.packages('ggplot2')로 설치하세요.")
  }
  
  IRF1   <- IRFobj$IRF
  Upper1 <- IRFobj$Upper
  Lower1 <- IRFobj$Lower
  
  # Optional second/third objects
  IRF2 <- if (!is.null(IRFobj2)) IRFobj2$IRF else NULL
  IRF3 <- if (!is.null(IRFobj3)) IRFobj3$IRF else NULL
  
  idx    <- match(variables, colnames(IRF1))
  idx_na <- idx[!is.na(idx)]
  valid  <- colnames(IRF1)[idx_na]
  
  if (length(valid) == 0L) {
    stop("`variables`에 해당하는 열이 IRF에 없습니다.")
  }
  
  ida <- match(valid, names(tcode))
  idt <- tcode[ida]
  
  nT1 <- nrow(IRF1)
  
  df_list1 <- vector("list", length(idx_na))
  df_list2 <- vector("list", length(idx_na))  # may remain empty if IRFobj2=NULL
  df_list3 <- vector("list", length(idx_na))  # may remain empty if IRFobj3=NULL
  
  for (k in seq_along(idx_na)) {
    j1 <- idx_na[k]
    
    ## (1) 차분/성장률 계열이면 누적 IRF (첫 번째 IRF)
    if (idt[k] %in% c(2L, 5L)) {
      y1_main  <- cumsum(IRF1[, j1])
      y1_upper <- cumsum(Upper1[, j1])
      y1_lower <- cumsum(Lower1[, j1])
    } else {
      y1_main  <- IRF1[, j1]
      y1_upper <- Upper1[, j1]
      y1_lower <- Lower1[, j1]
    }
    
    # Second series
    has2 <- FALSE
    if (!is.null(IRF2)) {
      j2 <- match(valid[k], colnames(IRF2))
      if (!is.na(j2)) {
        has2 <- TRUE
        if (idt[k] %in% c(2L, 5L)) {
          y2_main <- cumsum(IRF2[, j2])
        } else {
          y2_main <- IRF2[, j2]
        }
      }
    }
    
    # Third series
    has3 <- FALSE
    if (!is.null(IRF3)) {
      j3 <- match(valid[k], colnames(IRF3))
      if (!is.na(j3)) {
        has3 <- TRUE
        if (idt[k] %in% c(2L, 5L)) {
          y3_main <- cumsum(IRF3[, j3])
        } else {
          y3_main <- IRF3[, j3]
        }
      }
    }
    
    ## (2) 변수별로 최대 절대값 계산 (모든 시리즈를 함께 고려)
    vals_for_ct <- c(y1_main, y1_upper, y1_lower)
    if (has2) vals_for_ct <- c(vals_for_ct, y2_main)
    if (has3) vals_for_ct <- c(vals_for_ct, y3_main)
    max_abs_k <- max(abs(vals_for_ct), na.rm = TRUE)
    
    ## (3) 이 변수에서만 ±ct를 넘으면, 그 변수만 잘라서 사용
    if (max_abs_k > abs(ct)) {
      y1_main_plot  <- pmax(pmin(y1_main,  ct), -ct)
      y1_upper_plot <- pmax(pmin(y1_upper, ct), -ct)
      y1_lower_plot <- pmax(pmin(y1_lower, ct), -ct)
      if (has2) y2_main_plot <- pmax(pmin(y2_main, ct), -ct)
      if (has3) y3_main_plot <- pmax(pmin(y3_main, ct), -ct)
    } else {
      y1_main_plot  <- y1_main
      y1_upper_plot <- y1_upper
      y1_lower_plot <- y1_lower
      if (has2) y2_main_plot <- y2_main
      if (has3) y3_main_plot <- y3_main
    }
    
    # 첫 번째 IRF 데이터
    df_list1[[k]] <- data.frame(
      step     = seq_len(nT1),
      variable = valid[k],
      IRF      = y1_main_plot,
      Upper    = y1_upper_plot,
      Lower    = y1_lower_plot
    )
    
    # 두 번째 IRF 데이터 (bands 없음, 있으면)
    if (has2) {
      df_list2[[k]] <- data.frame(
        step     = seq_len(length(y2_main_plot)),
        variable = valid[k],
        IRF      = y2_main_plot
      )
    } else {
      df_list2[[k]] <- NULL
    }
    
    # 세 번째 IRF 데이터 (bands 없음, 있으면)
    if (has3) {
      df_list3[[k]] <- data.frame(
        step     = seq_len(length(y3_main_plot)),
        variable = valid[k],
        IRF      = y3_main_plot
      )
    } else {
      df_list3[[k]] <- NULL
    }
  }
  
  df1_all <- do.call(rbind, df_list1)
  df1_all$variable <- factor(df1_all$variable, levels = valid)
  
  have_df2 <- any(vapply(df_list2, function(x) !is.null(x), logical(1L)))
  have_df3 <- any(vapply(df_list3, function(x) !is.null(x), logical(1L)))
  
  if (have_df2) {
    df2_all <- do.call(rbind, df_list2)
    df2_all$variable <- factor(df2_all$variable, levels = valid)
  } else {
    df2_all <- NULL
  }
  
  if (have_df3) {
    df3_all <- do.call(rbind, df_list3)
    df3_all$variable <- factor(df3_all$variable, levels = valid)
  } else {
    df3_all <- NULL
  }
  
  library(ggplot2)
  
  # ---- (A) IRFobj2, IRFobj3 모두 없으면 기존과 동일 ----
  if ((is.null(IRFobj2) || is.null(df2_all)) &&
      (is.null(IRFobj3) || is.null(df3_all))) {
    
    p <- ggplot(df1_all, aes(x = step, y = IRF)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper),
                  fill = "#c4d7f2", alpha = 0.4) +
      geom_line(size = 0.7, colour = "#1f4e79") +
      geom_hline(yintercept = 0, size = 0.4, colour = "grey40") +
      facet_wrap(~ variable, ncol = ncol, scales = "free_y") +
      labs(x = "Months", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        axis.title.x     = element_text(margin = margin(t = 6)),
        plot.margin      = margin(5.5, 7, 5.5, 5.5)
      )
    
  } else {
    # ---- (B) 두 번째, 세 번째 IRF를 다른 색으로 오버레이 ----
    
    # 실제로 사용되는 색/레이블만 legend에 표시
    breaks_used <- c("model1",
                     if (have_df2) "model2" else NULL,
                     if (have_df3) "model3" else NULL)
    labels_used <- c(model1_name,
                     if (have_df2) model2_name else NULL,
                     if (have_df3) model3_name else NULL)
    
    p <- ggplot() +
      # Bands + first IRF
      geom_ribbon(
        data = df1_all,
        aes(x = step, ymin = Lower, ymax = Upper),
        fill = "#c4d7f2", alpha = 0.4
      ) +
      geom_line(
        data = df1_all,
        aes(x = step, y = IRF, colour = "model1"),
        size = 0.7
      ) +
      # Second IRF (no bands)
      { if (have_df2)
        geom_line(
          data = df2_all,
          aes(x = step, y = IRF, colour = "model2"),
          size = 0.7
        )
      } +
      # Third IRF (no bands)
      { if (have_df3)
        geom_line(
          data = df3_all,
          aes(x = step, y = IRF, colour = "model3"),
          size = 0.7
        )
      } +
      geom_hline(yintercept = 0, size = 0.4, colour = "grey40") +
      facet_wrap(~ variable, ncol = ncol, scales = "free_y") +
      scale_colour_manual(
        values = c("model1" = "#1f4e79",
                   "model2" = "#d55e00",
                   "model3" = "#009e73"),
        breaks = breaks_used,
        labels = labels_used,
        name   = NULL
      ) +
      labs(x = "Months", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        axis.title.x     = element_text(margin = margin(t = 6)),
        plot.margin      = margin(5.5, 7, 5.5, 5.5),
        legend.position  = "bottom"
      )
  }
  
  if (save_png) {
    if (is.null(file)) {
      obj_name <- deparse(substitute(IRFobj))
      file <- paste0(obj_name, ".pdf")
    }
    ggsave(
      filename = file,
      plot     = p, device = "pdf",
      width    = width,
      height   = height,
      units    = units,
      dpi      = res
    )
    invisible(file)
  } else {
    print(p)
    invisible(p)
  }
}




plotIRF <- function(IRFobj, dat, tcode, variables, ylim) {
  # IRFobj: list containing IRF, Upper, Lower matrices (nsteps x nvars)
  # dat    : data frame with column names matching variable names
  # tcode  : transformation codes (named vector)
  # variables: character vector of variable names to plot
  
  IRF   <- IRFobj$IRF
  Upper <- IRFobj$Upper
  Lower <- IRFobj$Lower
  
  idx <- match(variables, colnames(dat))
  idt <- tcode[variables]

  par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
  
    # cumulative if tcode == 2 || 5 (differenced)
    if (idt[i] == 2 || 5) {
      y_main   <- cumsum(IRF[, idx])
      y_upper  <- cumsum(Upper[, idx])
      y_lower  <- cumsum(Lower[, idx])
    } else {
      y_main   <- IRF[, idx]
      y_upper  <- Upper[, idx]
      y_lower  <- Lower[, idx]
    }
    
   if(missing(ylim)){ ylim = range(y_lower, y_upper, na.rm = TRUE);}
  
    plot(y_main, type = "l", lwd = 2, main = variables,
         ylab = "", xlab = "Steps",
         ylim = ylim, cex.main = 1.5, cex.axis = 1.3)
    
    lines(y_upper, lty = 2, col = "red")
    lines(y_lower, lty = 2, col = "red")
    abline(h = 0)
}




