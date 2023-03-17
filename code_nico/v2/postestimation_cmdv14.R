library(survey)
library(stargazer)
library(candisc)

#' Computes Canonical correlation p-values from regression coefficients
#'
#' @param data A data frame containing the variables to be used in the canonical correlation analysis
#' @param vars1 A vector of variable names for the first set of variables
#' @param vars2 A vector of variable names for the second set of variables
#' @param finwgt Name of the column that contains the survey weights
#' @param n_cc Number of canonical correlations to be computed
#' @param rawcoef_var1 A matrix containing the raw coefficients for the first set of variables
#' @param rawcoef_var2 A matrix containing the raw coefficients for the second set of variables
#' @ccorr A matrix containing the canonical correlations
#'
#' @return A data frame with the p-values for the canonical correlations

csdcanon <- function(data, vars1, vars2, finwgt, n_cc, rawcoef_var1, rawcoef_var2) {
    # determine which of vars1 and vars2 is the longest
    if (length(vars1) >= length(vars2)) {
        longvars <- vars1
        shortvars <- vars2
    } else {
        longvars <- vars2
        shortvars <- vars1
    }

    # Define OgX as the columns of data corresponding to longvars
    OgX <- data[, longvars]
    # Define OgY as the columns of data corresponding to shortvars
    OgY <- data[, shortvars]

    # Define finwgt as the column of data corresponding to finwgt
    finwgt <- data[, finwgt]
    diag_W <- diag(finwgt)

    # Standarize the colums of OgX and store the result in X
    X <- scale(OgX)
    # Standarize the colums of OgY and store the result in Y
    Y <- scale(OgY)

    OgX <- as.matrix(OgX)
    OgY <- as.matrix(OgY)
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    pvals <- calcpval(OgX, OgY, X, Y, diag_W, finwgt, n_cc, rawcoef_var1, rawcoef_var2)

    Results <- pvals$Results

    colnames(Results) <- c("CCIndex", "CC", "WilksLambda", "WilksLambdaFREQ", "RoysGreatestRoot", "PillaisTrace", "HotellingLawleyTrace", "WeightedReg", "CSDReg", "WeightedReg2")

    OGStataUVW <- pvals$OGStataUVW

    # calculate weightindex and set survey design using svydesign function from survey package
    weightindex <- (2 * n_cc) + 1

    des <- svydesign(weights = OGStataUVW[[weightindex]], ids = ~1, data = OGStataUVW)

    names_OGStataUVW <- names(OGStataUVW)
    # loop through 1 to n_cc
    for (i in 1:n_cc) {
        secondindex <- i + n_cc
        formula_str_1 <- paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])
        formula_str_2 <- paste(names_OGStataUVW[secondindex], "~", names_OGStataUVW[i])
        # run svyglm function from survey package
        model1 <- svyglm(as.formula(formula_str_1), design = des)
        # save coefficients to Results matrix
        Results[i, 8] <- summary(model1)$coefficients[2, 4]

        model2 <- svyglm(as.formula(formula_str_2), design = des)
        Results[i, 10] <- summary(model2)$coefficients[2, 4]
    }

    # create an empty rowfmt string in R
    rowfmt <- ""

    for (i in 1:n_cc) {
        secondindex <- i + n_cc
        formula_str_1 <- paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])
        # run svyglm function from survey package
        model <- svyglm(as.formula(formula_str_1), design = des)
        # save coefficients to Results matrix
        Results[i, 9] <- summary(model)$coefficients[2, 4]

        # adjust rowfmt accordingly
        if (i == n_cc) {
            rowfmt <- paste0(rowfmt, "-")
        } else {
            rowfmt <- paste0(rowfmt, "&")
        }
    }


    Results <- Results[, 1:9]

    # set column names of Results matrix
    colnames(Results) <- c("CC Index", "CC Magnitude", "Wilk's Lambda", "Wilk's Lambda (n)", "Roy's Greatest Root", "Pillai's Trace", "Hotelling-Lawley Trace", "Weighted Regression", "CSD Regression")

    # print matrix using 'stargazer' package

    return(Results)
}


#' Computes the p-values for the canonical correlations
#'
#' @param OgX A matrix containing the first set of variables
#' @param OgY A matrix containing the second set of variables
#' @param X A matrix containing the standardized first set of variables
#' @param Y A matrix containing the standardized second set of variables
#' @param diag_W A matrix containing the diagonal of the survey weights
#' @param finwgt A vector containing the survey weights
#' @param n_cc Number of canonical correlations to be computed
#' @param rawcoef_var1 A matrix containing the raw coefficients for the first set of variables
#' @param rawcoef_var2 A matrix containing the raw coefficients for the second set of variables
#'
#' @return A list containing the p-values for the canonical correlations
calcpval <- function(OgX, OgY, X, Y, diag_W, finwgt, n_cc, rawcoef_var1, rawcoef_var2) {
    output <- list()
    myResults <- matrix(0, nrow = n_cc, ncol = 10)

    myX <- X
    myY <- Y
    mydiag_W <- diag_W
    myOgX <- OgX
    myOgY <- OgY
    myrawcoef_var1 <- rawcoef_var1
    myrawcoef_var2 <- rawcoef_var2
    myW <- finwgt


    # calculate U and V in R
    U <- myOgX %*% myrawcoef_var1
    V <- myOgY %*% myrawcoef_var2

    # store U, V, and UVW matrices in R
    OGStataU <- U
    OGStataV <- V
    UV <- cbind(U, V)
    UVW <- cbind(UV, myW)
    output$OGStataUVW <- as.data.frame(UVW)

    # calculate weighted variance and correlation coefficient
    weighted_varU <- 97 * rep(1, ncol(U))
    weighted_varV <- 98 * rep(1, ncol(V))
    weighted_rho <- 99 * rep(1, ncol(U))

    for (i in 1:ncol(U)) {
        meanU <- mean(U[, i])
        meanV <- mean(V[, i])
        U[, i] <- U[, i] - meanU
        V[, i] <- V[, i] - meanV
        meanU <- mean(U[, i])
        meanV <- mean(V[, i])

        weighted_varU[i] <- t(U[, i]) %*% mydiag_W %*% U[, i] / sum(diag(mydiag_W))
        weighted_varV[i] <- t(V[, i]) %*% mydiag_W %*% V[, i] / sum(diag(mydiag_W))
        weighted_rho[i] <- (t(U[, i]) %*% mydiag_W %*% V[, i] / sum(diag(mydiag_W))) / sqrt(weighted_varU[i] * weighted_varV[i])

        myResults[i, 1] <- i
        myResults[i, 2] <- weighted_rho[i]
    }

    sumFREQ <- sum(trunc(myW))
    Lambda <- numeric(length(weighted_rho))
    df <- numeric(length(weighted_rho))
    ChiSq_Wilks_Lambda <- numeric(length(weighted_rho))
    ChiSq_Wilks_LambdaFREQ <- numeric(length(weighted_rho))
    p_val_Wilks_Lambda <- numeric(length(weighted_rho))
    p_val_Wilks_LambdaFREQ <- numeric(length(weighted_rho))

    for (i in 1:length(weighted_rho)) {
        Lambda[i] <- (1 - (weighted_rho[i]^2))
        for (j in (i + 1):length(weighted_rho)) {
            Lambda[i] <- Lambda[i] * (1 - (weighted_rho[j]^2))
        }

        df[i] <- (ncol(myX) + 1 - i) * (ncol(myY) + 1 - i)
        ChiSq_Wilks_Lambda[i] <- (((nrow(myX) - 1) - 0.5 * (ncol(myX) + ncol(myY) + 1)) * log(Lambda[i])) * -1
        ChiSq_Wilks_LambdaFREQ[i] <- (((sumFREQ - 1) - 0.5 * (ncol(myX) + ncol(myY) + 1)) * log(Lambda[i])) * -1
        p_val_Wilks_Lambda[i] <- pchisq(ChiSq_Wilks_Lambda[i], df[i], lower.tail = FALSE)
        p_val_Wilks_LambdaFREQ[i] <- pchisq(ChiSq_Wilks_LambdaFREQ[i], df[i], lower.tail = FALSE)
        myResults[i, 3] <- p_val_Wilks_Lambda[i]
        myResults[i, 4] <- p_val_Wilks_LambdaFREQ[i]
    }

    p <- ncol(myX)
    q <- ncol(myY)
    pq <- c(ncol(myX), ncol(myY))
    s <- min(pq)
    m <- (abs(p - q) - 1) / 2
    n <- (nrow(myX) - p - q - 2) / 2
    v1 <- s + (2 * m) + 1
    v2 <- s + (2 * n) + 1
    largest_root <- weighted_rho[1]^2
    Roys_Greatest_Root_stat <- rep(0, length(weighted_rho))
    Roys_Greatest_Root_stat[1] <- (v2 * largest_root) / (v1 * (1.0 - largest_root))
    Roys_Greatest_Root_p_value <- rep(0, length(weighted_rho))
    Roys_Greatest_Root_p_value[1] <- pf(Roys_Greatest_Root_stat[1], v1, v2, lower.tail = FALSE)

    myResults[1, 5] <- Roys_Greatest_Root_p_value[1]

    # Loop over roots k where k >= 2
    for (k in 2:length(weighted_rho)) {
        s <- min(pq) - k + 1
        m <- (abs(p - q - k + 1) - 1) / 2
        n <- (nrow(myX) - p - q - (k + 1)) / 2
        v1 <- s + (2 * m) + 1
        v2 <- s + (2 * n) + 1
        largest_root <- weighted_rho[k]^2
        Roys_Greatest_Root_stat[k] <- (v2 * largest_root) / (v1 * (1.0 - largest_root))
        Roys_Greatest_Root_p_value[k] <- pf(Roys_Greatest_Root_stat[k], v1, v2, lower.tail = FALSE)

        myResults[k, 5] <- Roys_Greatest_Root_p_value[k]
    }

    V <- 0
    for (i in 1:length(weighted_rho)) {
        V <- V + (weighted_rho[i]^2)
    }

    s <- min(pq)
    m <- (abs(p - q) - 1) / 2
    n <- (nrow(myX) - p - q - 2) / 2
    df1 <- s * ((2 * m) + s + 1)
    df2 <- s * ((2 * n) + s + 1)

    Pillais_Trace_stat <- rep(0, length(weighted_rho))
    Pillais_Trace_stat[1] <- (((2 * n) + s + 1) / ((2 * m) + s + 1)) * V / (s - V)
    Pillais_Trace_p_value <- rep(0, length(weighted_rho))
    Pillais_Trace_p_value[1] <- pf(Pillais_Trace_stat[1], df1, df2, lower.tail = FALSE)

    myResults[1, 6] <- Pillais_Trace_p_value[1]

    for (j in 2:length(weighted_rho)) {
        V <- 0
        for (i in j:length(weighted_rho)) {
            V <- V + (weighted_rho[i]^2)
        }

        Pillais_Trace_stat[j] <- (((2 * n) + s + 1) / ((2 * m) + s + 1)) * V / (s - V)
        Pillais_Trace_p_value[j] <- pf(Pillais_Trace_stat[j], df1, df2, lower.tail = FALSE)
        myResults[j, 6] <- Pillais_Trace_p_value[j]
    }

    Hotelling_Lawley_Trace_stat <- rep(0, length(weighted_rho))
    Hotelling_Lawley_Trace_p_value <- rep(0, length(weighted_rho))

    U <- 0
    for (i in 1:length(weighted_rho)) {
        U <- U + (weighted_rho[i]^2) / (1 - (weighted_rho[i]^2))
    }

    if (n > 0) {
        b <- (p + 2 * n) * (q + 2 * n) / (2 * (2 * n + 1) * (n - 1))
        c <- (2 + (p * q + 2) / (b - 1)) / (2 * n)
        mydf1 <- p * q
        mydf2 <- 4 + ((p * q + 2) / (b - 1))
        Hotelling_Lawley_Trace_stat[1] <- (U / c) * ((4 + (p * q + 2) / (b - 1)) / (p * q))
        Hotelling_Lawley_Trace_p_value[1] <- pf(Hotelling_Lawley_Trace_stat[1], mydf1, mydf2, lower.tail = FALSE)
        myResults[1, 7] <- Hotelling_Lawley_Trace_p_value[1]
    } else {
        mydf1 <- s * (2 * m + s + 1)
        mydf2 <- 2 * (s * n + 1)
        Hotelling_Lawley_Trace_stat[1] <- (2 * (s * n + 1) * U) / ((s^2) * (2 * m + s + 1))
        Hotelling_Lawley_Trace_p_value[1] <- pf(Hotelling_Lawley_Trace_stat[1], mydf1, mydf2, lower.tail = FALSE)
        myResults[1, 7] <- Hotelling_Lawley_Trace_p_value[1]
    }

    for (j in 2:length(weighted_rho)) {
        U <- 0
        for (i in j:length(weighted_rho)) {
            U <- U + (weighted_rho[i]^2) / (1 - (weighted_rho[i]^2))
        }

        if (n > 0) {
            b <- (p + 2 * n) * (q + 2 * n) / (2 * (2 * n + 1) * (n - 1))
            c <- (2 + (p * q + 2) / (b - 1)) / (2 * n)
            mydf1 <- p * q
            mydf2 <- 4 + (p * q + 2) / (b - 1)
            Hotelling_Lawley_Trace_stat[j] <- (U / c) * ((4 + (p * q + 2) / (b - 1)) / (p * q))
            Hotelling_Lawley_Trace_p_value[j] <- pf(Hotelling_Lawley_Trace_stat[j], mydf1, mydf2, lower.tail = FALSE)
            myResults[j, 7] <- Hotelling_Lawley_Trace_p_value[j]
        } else {
            mydf1 <- s * (2 * m + s + 1)
            mydf2 <- 2 * (s * n + 1)
            Hotelling_Lawley_Trace_stat[j] <- (2 * (s * n + 1) * U) / ((s^2) * (2 * m + s + 1))
            Hotelling_Lawley_Trace_p_value[j] <- pf(Hotelling_Lawley_Trace_stat[j], mydf1, mydf2, lower.tail = FALSE)
            myResults[j, 7] <- Hotelling_Lawley_Trace_p_value[j]
        }
    }

    output$Results <- myResults

    return(output)
}

data <- read.csv("BCC/toy/fit.csv")
cc <- cancor(cbind(Chins, Situps, Jumps) ~ Weight + Waist + Pulse, data = data, weights = Weight)

Results <- csdcanon(data, c("Chins", "Situps", "Jumps"), c("Weight", "Waist", "Pulse"), "Weight", cc$ndim, cc$coef$X, cc$coef$Y)

stargazer(Results,
    title = "Complete List of Canonical Correlation p-values", type = "text",
    column.labels = c("CC Index", "CC Magnitude", "Wilk's Lambda", "Wilk's Lambda (n)", "Roy's Greatest Root", "Pillai's Trace", "Hotelling-Lawley Trace", "Weighted Regression", "CSD Regression"),
    align = TRUE, digits = 4, column.sep.width = "-",
    row.sep = "-",
    float = FALSE
)
