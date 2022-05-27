alasso.cv <- function(x, y){

	p <- ncol(x)
	alasso <- lse <- numeric(p + 1)
	lm.out <- lm(y ~ x)
 	lse <- lm.out$coef
	xtilde <- sweep(x, 2, abs(lse[-1]), "*")

 	cv.out <- cv.glmnet(xtilde, y)
	lasso <- as.numeric(
		cv.out$glmnet.fit$beta[, which(cv.out$lambda == 
		cv.out$lambda.1se)])

	alasso[-1] <- lasso * abs(lse[-1])
  	alasso[1] <- mean(y - x %*% alasso[-1])
	alasso
}

stopes <- function(x, y, m = 20, prop_split = 0.50, prop_trim = 0.20, 
	q_tail = 0.90){

	n <- nrow(x)
	p <- ncol(x)
	full_data <- cbind(y, x)
	quan_NA <- 0

	p_univ_full <- rep(NA, p)
	p_univ <- matrix(NA, nrow = m, ncol = p)
	max_pval_splits <- min_pval_splits <- rep(NA, m)
	beta_stopes <- beta_pelt <- numeric(p + 1)

	cutpoints <- rep(NA, m * p)
	SS_sq <- matrix(NA, nrow = m * p, ncol = m)
	SS <- rep(NA, m * p)
	SS_splits <- matrix(NA, nrow = m, ncol = p + 1)

	for (j in 1 : p){

		p_univ_full[j] <- summary(lm(y ~ x[ , j]))$coefficients[2, 4]

	}

	for (jsim in 1 : m) {

		split <- sample(1 : n, round(prop_split * n), replace = FALSE)
		train_data <- full_data[split, ]
		test_data <- full_data[ - split, ]
		y_train <- train_data[ , 1]
		x_train <- train_data[ , 2 : (p + 1)]

		y_test <- test_data[ , 1]
		x_test <- test_data[ , 2 : (p + 1)]

		for (j in 1 : p){
		
			p_univ[jsim, j] <- summary(lm(y_train ~ x_train[ , j]))$coefficients[2, 4]

		}
    
		max_pval_splits[jsim] <- max(p_univ[jsim, ])
		min_pval_splits[jsim] <- min(p_univ[jsim, ])

		cutpoints[((jsim - 1) * p + 1) : (jsim * p)] <- p_univ[jsim, ]
		thresholds <- c(p_univ[jsim, order(p_univ[jsim, ])], 1.01) 

		fit_multiv <- summary(lm(y_train ~ 1))
		linearPred_splits <- rep(1, nrow(x_test)) * fit_multiv$coefficients[1, 1]
		SS_splits[jsim, 1] <- sum((y_test - linearPred_splits)^2)

		for (k in 2 : (p + 1)){

			signif <- which(p_univ[jsim, ] < thresholds[k])

			x_train_splits <- x_train[, signif, drop = F]
				fit_multiv <- summary(lm(y_train ~ x_train_splits))
   			linearPred_splits <- cbind(rep(1, nrow(x_test)), x_test[, signif]) %*% 
				fit_multiv$coefficients[ , 1]
			SS_splits[jsim, k] <- sum((y_test - linearPred_splits)^2, na.rm = TRUE)
		}
	}

	order <- order(cutpoints)

	for (jsim in 1 : m){

		cuts <- p_univ[jsim, order(p_univ[jsim, ])]
		
		if (cuts[1] == min(cutpoints)){

			SS_sq[1, jsim] <- SS_splits[jsim, 1]
			SS_sq[which(cutpoints[order] > cuts[1] & cutpoints[order] < cuts[2]), jsim] <- SS_splits[jsim, 2]

			for (k in 2 : (p - 1)) SS_sq[which(cutpoints[order] >= cuts[k] & cutpoints[order] < cuts[k + 1]), jsim] <- 
				SS_splits[jsim, k + 1]
		} else {

			SS_sq[which(cutpoints[order] < cuts[1]), jsim] <- SS_splits[jsim, 1]

			for (k in 1 : (p - 1)) SS_sq[which(cutpoints[order] >= cuts[k] & cutpoints[order] < cuts[k + 1]), jsim] <- 
				SS_splits[jsim, k + 1]
 		}

		SS_sq[which(cutpoints[order] >= cuts[p]), jsim] <- SS_splits[jsim, p + 1]
	}

	mean_trim <- function(x){
		p <- mean(x, trim = prop_trim, na.rm = TRUE)
		p
	}

	lower_quant <- quantile(min_pval_splits, probs = q_tail)
	upper_quant <- quantile(max_pval_splits, probs = 1 - q_tail)
	quan <- which(cutpoints[order] >= lower_quant & cutpoints[order] <= upper_quant)

	if(length(quan) > 1){
	
		SSsq <- SS_sq[quan, ]
		SS <- apply(SSsq, 1, mean_trim)
		minSS <- min(SS)
  
		SS_0_25percent <- minSS + .25 * sd(SSsq[which(SS == minSS)[1], ])
		final_cutpoints <- cutpoints[order][quan[which(SS <= SS_0_25percent)]][1]

		changepointPELT <- cpt.meanvar(SS, method = "PELT") 
		cptsPELT <- cpts(changepointPELT)
		cpts <- cptsPELT[1]
		cpts_neigh <- union(union(cpts, cpts + 1), cpts - 1)
		cpts_neigh = union(union(cpts - 1, cpts), cpts + 1)

		minSSchangepoint <- which(SS == min(SS[cpts_neigh]))
		final_cutpoints_PELT <- mean(cutpoints[order][quan[minSSchangepoint]])

	} else {

		quan_NA <- 1
		SS <- apply(SS_sq, 1, mean_trim)
		minSS <- min(SS)

		SS_0_25prop <- minSS + .25 * sd(SS_sq[which(SS == minSS)[1], ])  
		final_cutpoints <- cutpoints[order][which(SS <= SS_0_25prop)][1]

		changepointPELT <- cpt.meanvar(SS, method = "PELT")
		cptsPELT <- cpts(changepointPELT)
		cpts <- cptsPELT[1]
		cpts_neigh <- union(union(cpts, cpts + 1), cpts - 1)
		cpts_neigh <- union(union(cpts - 1, cpts), cpts + 1)

		minSSchangepoint <- which(SS == min(SS[cpts_neigh])) 

		final_cutpoints_PELT <- mean(cutpoints[order][minSSchangepoint])
	}

	J_stopes <- p_univ_full < final_cutpoints

	if(sum(J_stopes) >= 1){

		x_stopes <- x[ , J_stopes]
		stopes_fit <- lm(y ~ x_stopes)
		beta_stopes[c(TRUE, J_stopes)] <- coef(stopes_fit)
		beta_stopes[c(FALSE, !J_stopes)] <- 0 
	} else{

		stopes_fit <- lm(y ~ 1)
		beta_stopes[c(TRUE, J_stopes)] <- coef(stopes_fit)
		beta_stopes[c(FALSE, !J_stopes)] <- 0 
	}

	J_pelt <- p_univ_full < final_cutpoints_PELT

	if(sum(J_pelt) >= 1){

		x_pelt <- x[ , J_pelt]
		pelt_fit <- lm(y ~ x_pelt)
		beta_pelt[c(TRUE, J_pelt)] <- coef(pelt_fit)
		beta_pelt[c(FALSE, !J_pelt)] <- 0 
	} else{

		pelt_fit <- lm(y ~ 1)
		beta_pelt[c(TRUE, J_pelt)] <- coef(pelt_fit)
		beta_pelt[c(FALSE, !J_pelt)] <- 0
	}

	list(beta_stopes = beta_stopes, J_stopes = J_stopes, 
		beta_pelt = beta_pelt, J_pelt = J_pelt, 
		final_cutpoints = final_cutpoints,
		final_cutpoints_PELT = final_cutpoints_PELT,
		quan_NA = quan_NA
	)
}

opts_th <- function(X, Y, m, crit = "aic", type = "binseg", prop_split = 0.5, prop_trim = 0.2, q_tail = 0.5, ...) {

  n <- nrow(X)
  p <- ncol(X)
  betahat <- SE <- rep(NA, p + 1)
  Jhat <- rep(NA, p)
  pval <- rep(NA, p)
  pval_split <- matrix(NA, nrow = m, ncol = p)
  pval_max <- pval_min <- rep(NA, m)
  
  for(j in 1 : p){
    
    glm.outj <- glm(Y ~ X[ , j], ...)
    pval[j] <- summary(glm.outj)$coef[2, 4]
  }

  if(crit == "aic"){

    aic_split <- matrix(NA, nrow = m, ncol = p + 1)
    aic_split_int <- matrix(NA, nrow = m, ncol = m * p)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pval_split[i, j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pval_min[i] <- min(pval_split[i, ])
      pval_max[i] <- max(pval_split[i, ])
      pval_sort <- sort(pval_split[i, ])
    
      glm.out0 <- glm(Y_split ~ 1, ...)
      aic_split[i, 1] <- AIC(glm.out0)
      
      for (k in 1 : p) {

        Jhatk <- pval_split[i, ] <= pval_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        aic_split[i, k + 1] <- AIC(glm.outk)
      }
    
    }
    
    cutpoints <- sort(pval_split)
    
    for(i in 1:m){
      
      cuts <- sort(pval_split[i, ])
      
      if(min(cutpoints) == cuts[1]){
        
        aic_split_int[i, cutpoints == cuts[1]] <- aic_split[i, 1]
        aic_split_int[i, cuts[1] < cutpoints & cutpoints < cuts[2]] <- aic_split[i, 2]
        
        for(k in 2 : (p - 1)){
          
          aic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            aic_split[i, k + 1]
        }
      } else{
        
        aic_split_int[i, cutpoints < cuts[1]] <- aic_split[i, 1]
        for(k in 1 : (p - 1)){
          
          aic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            aic_split[i, k + 1]
        }
        
      }
      
      aic_split_int[i, cuts[p] <= cutpoints] <- aic_split[i, p + 1]
    }
    
    lower_quant <- quantile(pval_min, probs = q_tail)
    upper_quant <- quantile(pval_max, probs = 1 - q_tail)
    subset_quant <- lower_quant <= cutpoints & cutpoints <= upper_quant
    
    if(sum(subset_quant) >= 1){
      
      cutpoints <- cutpoints[subset_quant]
      aic_split_int <- aic_split_int[, subset_quant]
    }
    
    aic_mean <- apply(aic_split_int, 2, mean, trim = prop_trim, na.rm = TRUE)
    minaic.out <- sort(aic_mean, index.return = TRUE)
      
    # min AIC rule
    if(type == "min"){
      cuthat <- cutpoints[minaic.out$ix[1]]
    }
      
    # min AIC + 0.5 SD rule
    if(type == "sd"){
      minsd <- min(aic_mean) + 0.5 * sd(aic_split_int[, minaic.out$ix[1]])
      cuthat <- cutpoints[aic_mean <= minsd][1]
    }
    
    # pelt
    if(type == "pelt"){
      pelt.out <- cpt.meanvar(aic_mean, method = "PELT")
      cpts_pelt <- cpts(pelt.out)
      cpts_neigh <- c(cpts_pelt[1], cpts_pelt[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    # binseg
    if(type == "binseg"){
      bs.out <- cpt.meanvar(aic_mean, method = "BinSeg")
      cpts_bs <- cpts(bs.out)
      cpts_neigh <- c(cpts_bs[1], cpts_bs[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    # amoc
    if(type == "amoc"){
      amoc.out <- cpt.meanvar(aic_mean, method = "AMOC")
      cpts_amoc <- cpts(amoc.out)
      cpts_neigh <- c(cpts_amoc, cpts_amoc + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
    
    Jhat <- pval <= cuthat
    
    if(sum(Jhat) >= 1){
      
      XB <- X[, Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    } else{
      
      glm.out <- glm(Y ~ 1, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
    output <- list(betahat = betahat, Jhat = Jhat, SE = SE, cuthat = cuthat, 
      pval = pval, cutpoints = cutpoints, aic_mean = aic_mean)
    
  }
  
  if(crit == "bic"){

    bic_split <- matrix(NA, nrow = m, ncol = p + 1)
    bic_split_int <- matrix(NA, nrow = m, ncol = m * p)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pval_split[i, j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pval_min[i] <- min(pval_split[i, ])
      pval_max[i] <- max(pval_split[i, ])
      pval_sort <- sort(pval_split[i, ])
    
      glm.out0 <- glm(Y_split ~ 1, ...)
      bic_split[i, 1] <- BIC(glm.out0)
      
      for (k in 1 : p) {

        Jhatk <- pval_split[i, ] <= pval_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        bic_split[i, k + 1] <- BIC(glm.outk)
      }
    
    }
    
    cutpoints <- sort(pval_split)
    
    for(i in 1:m){
      
      cuts <- sort(pval_split[i, ])
      
      if(min(cutpoints) == cuts[1]){
        
        bic_split_int[i, cutpoints == cuts[1]] <- bic_split[i, 1]
        bic_split_int[i, cuts[1] < cutpoints & cutpoints < cuts[2]] <- bic_split[i, 2]
        
        for(k in 2 : (p - 1)){
          
          bic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            bic_split[i, k + 1]
        }
      } else{
        
        bic_split_int[i, cutpoints < cuts[1]] <- bic_split[i, 1]
        for(k in 1 : (p - 1)){
          
          bic_split_int[i, cuts[k] <= cutpoints & cutpoints < cuts[k + 1]] <- 
            bic_split[i, k + 1]
        }
        
      }
      
      bic_split_int[i, cuts[p] <= cutpoints] <- bic_split[i, p + 1]
    }
    
    lower_quant <- quantile(pval_min, probs = q_tail)
    upper_quant <- quantile(pval_max, probs = 1 - q_tail)
    subset_quant <- lower_quant <= cutpoints & cutpoints <= upper_quant
    
    if(sum(subset_quant) >= 1){
      
      cutpoints <- cutpoints[subset_quant]
      bic_split_int <- bic_split_int[, subset_quant]
    }
    
    bic_mean <- apply(bic_split_int, 2, mean, trim = prop_trim, na.rm = TRUE)
    minbic.out <- sort(bic_mean, index.return = TRUE)
      
    # min BIC rule
    if(type == "min"){
      cuthat <- cutpoints[minbic.out$ix[1]]
    }
      
    # min BIC + 0.5 SD rule
    if(type == "sd"){
      minsd <- min(bic_mean) + 0.5 * sd(bic_split_int[, minbic.out$ix[1]])
      cuthat <- cutpoints[bic_mean <= minsd][1]
    }
    
    # pelt
    if(type == "pelt"){
      pelt.out <- cpt.meanvar(bic_mean, method = "PELT")
      cpts_pelt <- cpts(pelt.out)
      cpts_neigh <- c(cpts_pelt[1], cpts_pelt[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
    
    # binseg
    if(type == "binseg"){
      bs.out <- cpt.meanvar(bic_mean, method = "BinSeg")
      cpts_bs <- cpts(bs.out)
      cpts_neigh <- c(cpts_bs[1], cpts_bs[1] + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    # amoc
    if(type == "amoc"){
      amoc.out <- cpt.meanvar(bic_mean, method = "AMOC")
      cpts_amoc <- cpts(amoc.out)[1]
      cpts_neigh <- c(cpts_amoc, cpts_amoc + 1)
      cuthat <- mean(cutpoints[cpts_neigh])
    }
      
    Jhat <- pval <= cuthat
    
    if(sum(Jhat) >= 1){
      
      XB <- X[, Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    } else{
      
      glm.out <- glm(Y ~ 1, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[,1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
    output <- list(betahat = betahat, Jhat = Jhat, SE = SE, cuthat = cuthat, 
      pval = pval, cutpoints = cutpoints, bic_mean = bic_mean)
  }
  
  output
}

opts <- function(X, Y, m, crit = "aic", prop_split = 0.5, cutoff = 0.75, ...) {

  n <- nrow(X)
  p <- ncol(X)
  betahat <- SE <- rep(NA, p + 1)
  pvals <- rep(NA, p)
  Jhat_split <- matrix(NA, nrow = m, ncol = p)

  if(crit == "aic"){

    aic_split <- rep(NA, p + 1)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pvals[j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pvals_sort <- sort(pvals)
      
      glm.out0 <- glm(Y_split ~ 1, ...)
      aic_split[1] <- AIC(glm.out0)
        
      for (k in 1 : p) {

        Jhatk <- pvals <= pvals_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        aic_split[k + 1] <- AIC(glm.outk)
      }
    
      aic_sort <- sort(aic_split[-1], index.return = TRUE)
      idx_min <- aic_sort$ix[1]
      aic_sort0 <- sort(aic_split, index.return = TRUE)
      idx_min0 <- aic_sort0$ix[1]
      
      if(idx_min0 != 1){
        
        Jhat_split[i, ] <- pvals <= pvals_sort[idx_min]
        
      } else {
        
        Jhat_split[i, ] <- rep(FALSE, p)
      }

    }

    freqs <- apply(Jhat_split, 2, mean)
    Jhat <- freqs >= cutoff

    if(sum(Jhat) >= 1){
    
      XB <- X[ , Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
    
    } else{
    
      glm.out <- glm(Y ~ 1, ...)
      Jhat <- rep(FALSE, p)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
  }
  
  if(crit == "bic"){

    bic_split <- rep(NA, p + 1)

    for (i in 1 : m) {

      sind <- sample(1 : n, round(prop_split * n), replace = FALSE)
      Y_split <- Y[sind]
      X_split <- X[sind, , drop = FALSE]

      for (j in 1 : p){ 
      
        glm.outj <- glm(Y_split ~ X_split[ , j], ...)
        pvals[j] <- summary(glm.outj)$coef[2, 4]
      }
    
      pvals_sort <- sort(pvals)
      
      glm.out0 <- glm(Y_split ~ 1, ...)
      bic_split[1] <- BIC(glm.out0)
        
      for (k in 1 : p) {

        Jhatk <- pvals <= pvals_sort[k]
        glm.outk <- glm(Y_split ~ X_split[ , Jhatk, drop = FALSE], ...)
        bic_split[k + 1] <- BIC(glm.outk)
      }
    
      bic_sort <- sort(bic_split[-1], index.return = TRUE)
      idx_min <- bic_sort$ix[1]
      bic_sort0 <- sort(bic_split, index.return = TRUE)
      idx_min0 <- bic_sort0$ix[1]
      
      if(idx_min0 != 1){
        
        Jhat_split[i, ] <- pvals <= pvals_sort[idx_min]
        
      } else {
        
        Jhat_split[i, ] <- rep(FALSE, p)
      }

    }

    freqs <- apply(Jhat_split, 2, mean)
    Jhat <- freqs >= cutoff

    if(sum(Jhat) >= 1){
    
      XB <- X[ , Jhat, drop = FALSE]
      glm.out <- glm(Y ~ XB, ...)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
    
    } else{
    
      glm.out <- glm(Y ~ 1, ...)
      Jhat <- rep(FALSE, p)
      JB <- c(TRUE, Jhat)
      betahat[JB] <- summary(glm.out)$coef[, 1]
      betahat[!JB] <- 0
      SE[JB] <- summary(glm.out)$coef[, 2]
      SE[!JB] <- 0
      
    }
    
  }
    
  list(betahat = betahat, Jhat = Jhat, SE = SE, freqs = freqs)
}

