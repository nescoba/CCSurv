library(candisc)
library(survey)
library(dplyr)

# Define a function named csdcanon
csdcanon <- function(OgX, OgY, finwgt, n_cc, rawcoef_var1, rawcoef_var2) {
    # Create diagonal matrix with survey weights
    W <- finwgt
    diag_W <- diag(W)

    X <- OgX
    Y <- OgY

    # Standardize matrices X and Y
    X <- scale(X)
    Y <- scale(Y)

    # Calculate canonical correlations using the 'matlib' package
    pvals <- calcpval(OgX, OgY, X, Y, diag_W, W, N = 5, rawcoef_var1, rawcoef_var2)
    OGStataUVW <- pvals$OGStataUVW

    # Define weight index
    weightindex <- 2 * n_cc + 1

    reg1 <- list()
    reg2 <- list()

    # Perform survey regression using the 'survey' package
    design <- svydesign(weights = OGStataUVW[, weightindex], data = OGStataUVW, ids = ~1)

    names_OGStataUVW <- names(OGStataUVW)
    for (i in 1:n_cc) {
        secondindex <- i + n_cc
        formula_str_1 <- paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])
        formula_str_2 <- paste(names_OGStataUVW[secondindex], "~", names_OGStataUVW[i])
        reg1[[i]] <- svyglm(as.formula(formula_str_1), design)
        reg2[[i]] <- svyglm(as.formula(formula_str_2), design)
    }

    return(list(reg1 = reg1, reg2 = reg2, pvals = pvals))
}


calcpval <- function(OgX, OgY, X, Y, diag_W, W, N, rawcoef_var1, rawcoef_var2) {
    return_list <- list()

    myX <- X
    myY <- Y
    mydiag_W <- diag_W
    myW <- W

    myrawcoef_var1 <- rawcoef_var1
    myrawcoef_var2 <- rawcoef_var2

    myOgX <- OgX
    myOgY <- OgY

    U <- myOgX %*% myrawcoef_var1
    V <- myOgY %*% myrawcoef_var2

    OGStataU <- U
    OGStataV <- V

    UV <- cbind(U, V)
    UVW <- cbind(UV, myW)
    OGStataUVW <- as.data.frame(UVW)

    return_list$OGStataUVW <- OGStataUVW

    # Using Eq. (8) and (9) to find the correlations and variances of the canonical matrices
    # This calculations do not lead to correct canonical correlations or variances of the canonical variates as it did when using all my process above
    weighted_varU <- numeric(ncol(U))
    weighted_varV <- numeric(ncol(V))
    weighted_rho <- numeric(ncol(U))

    # Subtract the mean from each column of U and V
    U <- sweep(U, 2, colMeans(U), "-")
    V <- sweep(V, 2, colMeans(V), "-")

    for (i in 1:ncol(U)) {
        weighted_varU[i] <- (t(U[, i]) %*% mydiag_W %*% U[, i]) / sum(mydiag_W)
        weighted_varV[i] <- (t(V[, i]) %*% mydiag_W %*% V[, i]) / sum(mydiag_W)
        weighted_rho[i] <- ((t(U[, i]) %*% mydiag_W %*% V[, i]) / sum(mydiag_W)) / sqrt(weighted_varU[i] * weighted_varV[i])
    }

    return_list$weighted_rho <- weighted_rho


    # Roy's Greatest Root

    Lambda <- numeric(length(weighted_rho))
    df <- numeric(length(weighted_rho))
    ChiSq_Wilks_Lambda <- numeric(length(weighted_rho))
    p_val_Wilks_Lambda <- numeric(length(weighted_rho))
    # p_val_Wilks_LambdaFREQ <- rep(1, length(weighted_rho)) # commented out

    for (i in 1:length(weighted_rho)) {
        Lambda[i] <- (1 - (weighted_rho[i]^2))
        for (j in (i + 1):length(weighted_rho)) {
            Lambda[i] <- Lambda[i] * (1 - (weighted_rho[j]^2))
        }
        df[i] <- (ncol(myX) + 1 - i) * (ncol(myY) + 1 - i)
        ChiSq_Wilks_Lambda[i] <- (((nrow(myX) - 1) - 0.5 * (ncol(myX) + ncol(myY) + 1)) * log(Lambda[i])) * (-1)
        p_val_Wilks_Lambda[i] <- pchisq(ChiSq_Wilks_Lambda[i], df[i], lower.tail = FALSE)
    }

    return_list$Lambda <- Lambda
    return_list$p_val_Wilks_Lambda <- p_val_Wilks_Lambda

    p <- ncol(myX)
    q <- ncol(myY)
    pq <- matrix(c(ncol(myX), ncol(myY)), nrow = 1)
    s <- min(pq)
    m <- (abs(p - q) - 1) / 2
    n <- (nrow(myX) - p - q - 2) / 2
    v1 <- s + (2 * m) + 1
    v2 <- s + (2 * n) + 1
    largest_root <- weighted_rho[1]^2
    Roys_Greatest_Root_stat <- numeric(length(weighted_rho))
    Roys_Greatest_Root_stat[1] <- (v2 * largest_root) / (v1 * (1 - largest_root))
    Roys_Greatest_Root_p_value <- numeric(length(weighted_rho))
    Roys_Greatest_Root_p_value[1] <- pf(Roys_Greatest_Root_stat[1], v1, v2, lower.tail = FALSE)

    # For the roots k where k>=2 we use the formulas at the bottom of page 61 of Gittins
    for (j in 2:length(weighted_rho)) {
        s <- min(pq) - j - 1
        m <- (abs(p - q - j) - 1) / 2
        n <- (nrow(myX) - p - q - j) / 2
        v1 <- s + (2 * m) + 1
        v2 <- s + (2 * n) + 1
        largest_root <- weighted_rho[j]^2
        Roys_Greatest_Root_stat[j] <- (v2 * largest_root) / (v1 * (1 - largest_root))
        Roys_Greatest_Root_p_value[j] <- pf(Roys_Greatest_Root_stat[j], v1, v2, lower.tail = FALSE)
    }

    return_list$Roys_Greatest_Root_stat <- Roys_Greatest_Root_stat
    return_list$Roys_Greatest_Root_p_value <- Roys_Greatest_Root_p_value

    # Pillai's Trace

    V <- 0
    for (i in 1:length(weighted_rho)) {
        V <- V + (weighted_rho[i]^2)
    }
    # From https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_introreg_sect038.htm#statug_introreg001883
    s <- min(pq)
    m <- (abs(p - q) - 1) / 2
    n <- (nrow(myX) - p - q - 2) / 2 # In the SAS help n=(v-p-1)/2 where v= error degrees of freedom
    df1 <- s * ((2 * m) + s + 1)
    df2 <- s * ((2 * n) + s + 1)
    # The degrees of freedom are the same as in SAS 9 & 48
    Pillais_Trace_stat <- (((2 * n) + s + 1) / ((2 * m) + s + 1)) * V / (s - V)
    Pillais_Trace_p_value <- pf(Pillais_Trace_stat, df1 = df1, df2 = df2, lower.tail = FALSE)

    return_list$Pillais_Trace_stat <- Pillais_Trace_stat
    return_list$Pillais_Trace_p_value <- Pillais_Trace_p_value


    # Hotelling-Lawley Trace

    U <- 0
    for (i in 1:length(weighted_rho)) {
        U <- U + (weighted_rho[i]^2) / (1 - (weighted_rho[i]^2))
    }

    if (n > 0) {
        b <- (p + 2 * n) * (q + 2 * n) / (2 * (2 * n + 1) * (n - 1))
        c <- (2 + (p * q + 2) / (b - 1)) / (2 * n)
        mydf1 <- p * q
        mydf2 <- 4 + (p * q + 2) / (b - 1)
        Hotelling_Lawley_Trace_stat <- (U / c) * ((4 + (p * q + 2) / (b - 1)) / (p * q))
        Hotelling_Lawley_Trace_p_value <- pf(Hotelling_Lawley_Trace_stat, mydf1, mydf2, lower.tail = FALSE)
    } else {
        mydf1 <- s * (2 * m + s + 1)
        mydf2 <- 2 * (s * n + 1)
        Hotelling_Lawley_Trace_stat <- (2 * (s * n + 1) * U) / ((s^2) * (2 * m + s + 1))
        Hotelling_Lawley_Trace_p_value <- pf(Hotelling_Lawley_Trace_stat, mydf1, mydf2, lower.tail = FALSE)
    }

    return_list$Hotelling_Lawley_Trace_stat <- Hotelling_Lawley_Trace_stat
    return_list$Hotelling_Lawley_Trace_p_value <- Hotelling_Lawley_Trace_p_value

    return(return_list)
}



data <- read.csv("BCC/toy/fit.csv")
cc <- cancor(cbind(Chins, Situps, Jumps) ~ Weight + Waist + Pulse, data = data, weights = Weight)

csd <- csdcanon(as.matrix(data[, 1:3]), as.matrix(data[, 4:6]), data$Weight, cc$ndim, cc$coef$X, cc$coef$Y)

print(csd$reg1)
print(csd$reg2)
print(csd$pvals$weighted_rho)
print(csd$pvals$Lambda)
print(csd$pvals$p_val_Wilks_Lambda)
print(csd$pvals$Roys_Greatest_Root_stat)
print(csd$pvals$Roys_Greatest_Root_p_value)
print(csd$pvals$Pillais_Trace_stat)
print(csd$pvals$Pillais_Trace_p_value)
print(csd$pvals$Hotelling_Lawley_Trace_stat)
print(csd$pvals$Hotelling_Lawley_Trace_p_value)
