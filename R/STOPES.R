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
