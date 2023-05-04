library(candisc)
library(survey)
library(stargazer)

plot_canon <- function(cc_object, dim1, dim2) {
    s1 <- diag(sqrt(diag(cov(cc_object$X))))
    mycoef1 <- s1 %*% cc_object$coef$X
    stdU <- cc_object$X %*% mycoef1
    s2 <- diag(sqrt(diag(cov(cc_object$Y))))
    mycoef2 <- s2 %*% cc_object$coef$Y
    stdV <- cc_object$Y %*% mycoef2

    graphcoef <- rbind(mycoef1, mycoef2)
    graphcoef <- cbind(graphcoef, c(rep(0, length(mycoef1)), rep(1, length(mycoef2))))

    # grid for two plots
    par(mfrow = c(1, 2))

    # plot a unit circle
    x <- seq(-1, 1, length = 100)
    y <- sqrt(1 - x^2)
    plot(x, y, type = "l", xlab = "", ylab = "", lty = 2, xlim = c(-2, 2), ylim = c(-2, 2))
    lines(x, -y, lty = 2)

    points(
        graphcoef[, dim2], graphcoef[, dim1],
        col = ifelse(graphcoef[, ncol(graphcoef) - 1] == 0, "red", "blue"),
        pch = ifelse(graphcoef[, ncol(graphcoef) - 1], 22, 4), xlab = paste0("CC", dim1),
        ylab = paste0("CC", dim2)
    )
    text(
        graphcoef[, dim2], graphcoef[, dim1],
        labels = c(cc_object$names$X, cc_object$names$Y)
    )

    stdV_dim1 <- stdV[dim1, ]
    stdU_dim1 <- stdU[dim1, ]

    plot(stdU_dim1, stdV_dim1, xlab = paste0("U", dim1), ylab = paste0("V", dim1))
    text(stdU_dim1, stdV_dim1, labels = 1:length(stdU_dim1))

    par(mfrow = c(1, 1))
}

csdcanon <- function(cc_object, howmany = NA, dim1, dim2, selection, command = "classic", pweight = NA, psu = NA, strata = NA) {
    plot_canon(cc_object, dim1, dim2)

    varcount1 <- length(cc_object$names$X)
    varcount2 <- length(cc_object$names$Y)

    if (varcount1 >= varcount2) {
        OgX <- cc_object$X
        OgY <- cc_object$Y
    } else {
        OgX <- cc_object$Y
        OgY <- cc_object$X
    }

    if (command != "classic") {
        diag_W <- diag(cc_object$weights)
    } else {
        diag_W <- diag(nrow = nrow(OgX))
    }

    # standarize the OgX and OgY
    X <- scale(OgX)
    Y <- scale(OgY)

    n_cc <- cc_object$ndim

    output <- calcpval(cc_object, selection, OgX, OgY, diag_W)

    Results <- output$Results
    OGStataUVW <- output$OGStataUVW

    colnames(Results) <- c("Statistic", "df1", "df2", "Chi-Sq/F", "p-val", "Index")

    weightindex <- (2 * n_cc) + 1 # Reading the weights as they were stored in this column in the mata function

    names_OGStataUVW <- names(OGStataUVW)

    if (command != "classic") {
        design_object <- svydesign(ids = ~1, weights = ~ OGStataUVW[, weightindex], data = OGStataUVW)
        # svyset(pweight = paste0("OGStataUVW", weightindex))

        # Notice coefficients of the simple linear regression are equal to the canonical correlations because the variances of the canonical variates are equal to 1
        for (i in 1:n_cc) {
            secondindex <- i + n_cc
            svyreg1 <- svyglm(as.formula(paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])), design = design_object)
            p1 <- coef(summary(svyreg1))[1, 4]

            svyreg2 <- svyglm(as.formula(paste(names_OGStataUVW[secondindex], "~", names_OGStataUVW[i])), design = design_object)
            p2 <- coef(summary(svyreg1))[1, 4]

            if (p1 >= p2) {
                Results[(7 * (i - 1)) + 6, 1] <- coef(summary(svyreg1))[1, 1]
                Results[(7 * (i - 1)) + 6, 2] <- svyreg1$df.residual + length(svyreg1$coefficients) - 1
                Results[(7 * (i - 1)) + 6, 3] <- svyreg1$df.null - svyreg1$df.residual
                aux <- summary(svyreg1)$coefficients[, 1]
                Results[(7 * (i - 1)) + 6, 4] <- aux[1]
                Results[(7 * (i - 1)) + 6, 5] <- summary(svyreg1)$coef[2, 4]
                Results[(7 * (i - 1)) + 6, 6] <- i
            } else {
                Results[(7 * (i - 1)) + 6, 1] <- coef(summary(svyreg2))[1, 1]
                Results[(7 * (i - 1)) + 6, 2] <- svyreg2$df.residual + length(svyreg2$coefficients) - 1
                Results[(7 * (i - 1)) + 6, 3] <- svyreg2$df.null - svyreg2$df.residual
                Results[(7 * (i - 1)) + 6, 4] <- summary(svyreg2)$fstatistic[1]
                Results[(7 * (i - 1)) + 6, 5] <- summary(svyreg2)$coef[2, 4]
                Results[(7 * (i - 1)) + 6, 6] <- i
            }
        }
    } else {
        design_object <- svydesign(ids = ~1, weights = ~1, data = OGStataUVW)
        for (i in 1:n_cc) {
            secondindex <- i + n_cc
            svyreg1 <- svyglm(as.formula(paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])), design = design_object)
            p1 <- coef(summary(svyreg))[1, 4]

            Results[(7 * (i - 1)) + 6, 1] <- coef(summary(svyreg1))[1, 1]
            Results[(7 * (i - 1)) + 6, 2] <- svyreg1$df.residual + length(svyreg1$coefficients) - 1
            Results[(7 * (i - 1)) + 6, 3] <- svyreg1$df.null - svyreg1$df.residual
            aux <- summary(svyreg1)$coefficients[, 1]
            Results[(7 * (i - 1)) + 6, 4] <- aux[1]
            # Results[(7 * (i - 1)) + 6, 4] <- summary(svyreg1)$fstatistic[1]
            Results[(7 * (i - 1)) + 6, 5] <- summary(svyreg1)$coef[2, 4]
            Results[(7 * (i - 1)) + 6, 6] <- i
        }
    }

    if (command != "classic") {
        design_object <- svydesign(ids = ~1, weights = ~ OGStataUVW[, weightindex], data = OGStataUVW)
        # svyset(pweight = paste0("OGStataUVW", weightindex))

        # Notice coefficients of the simple linear regression are equal to the canonical correlations because the variances of the canonical variates are equal to 1
        for (i in 1:n_cc) {
            secondindex <- i + n_cc
            svyreg1 <- svyglm(as.formula(paste(names_OGStataUVW[i], "~", names_OGStataUVW[secondindex])), design = design_object)
            p1 <- coef(summary(svyreg1))[1, 4]

            svyreg2 <- svyglm(as.formula(paste(names_OGStataUVW[secondindex], "~", names_OGStataUVW[i])), design = design_object)
            p2 <- coef(summary(svyreg1))[1, 4]

            if (p1 >= p2) {
                Results[(7 * (i - 1)) + 6, 1] <- coef(summary(svyreg1))[1, 1]
                Results[(7 * (i - 1)) + 6, 2] <- svyreg1$df.residual + length(svyreg1$coefficients) - 1
                Results[(7 * (i - 1)) + 6, 3] <- svyreg1$df.null - svyreg1$df.residual

                aux <- summary(svyreg1)$coefficients[, 1]
                Results[(7 * (i - 1)) + 6, 4] <- aux[1]

                # Results[(7 * (i - 1)) + 6, 4] <- summary(svyreg1)$fstatistic[1]
                Results[(7 * (i - 1)) + 6, 5] <- summary(svyreg1)$coef[2, 4]
                Results[(7 * (i - 1)) + 6, 6] <- i
            } else {
                Results[(7 * (i - 1)) + 6, 1] <- coef(summary(svyreg2))[1, 1]
                Results[(7 * (i - 1)) + 6, 2] <- svyreg2$df.residual + length(svyreg2$coefficients) - 1
                Results[(7 * (i - 1)) + 6, 3] <- svyreg2$df.null - svyreg2$df.residual
                aux <- summary(svyreg1)$coefficients[, 1]
                Results[(7 * (i - 1)) + 6, 4] <- aux[1]
                # Results[(7 * (i - 1)) + 6, 4] <- summary(svyreg2)$fstatistic[1]
                Results[(7 * (i - 1)) + 6, 5] <- summary(svyreg2)$coef[2, 4]
                Results[(7 * (i - 1)) + 6, 6] <- i
            }
        }
    }

    if (is.null(howmany)) {
        howmany <- n_cc
    }

    for (i in 1:howmany) {
        # Selecting rows in information for Canonical correlation i
        aux <- Results[(7 * (i - 1) + 1), 1]
        mag <- round(aux[1], 0.0001)
        ResultsAux <- Results[(7 * (i - 1) + 2):(7 * (i - 1) + 7), 1:5]
        rownames(ResultsAux) <- c("Wilks' Lambda", "Pillai's Trace", "Hotelling-Lawley Trace", "Roy's Greatest Root", "Weighted Reg", "Complex Survey Design Reg")
        print(ResultsAux, title = paste0("Statistics for Canonical Correlation: ", i), twidth = 30, format = "%10.4f", rowtitle = paste0("Canonical Correlation=", mag), border = c("top", "bottom"))
    }

    # return(Results)
}

calcpval <- function(cc_object, selection, OgX, OgY, diag_W) {
    output <- list()
    myResults <- matrix(NA, nrow = 7 * ncol(OgX), ncol = 6)

    if (length(cc_object$names$Y) > length(cc_object$names$X)) {
        myrawcoef_var2 <- cc_object$coef$X
        myrawcoef_var1 <- cc_object$coef$Y
    } else {
        myrawcoef_var1 <- cc_object$coef$X
        myrawcoef_var2 <- cc_object$coef$Y
    }

    U <- OgX %*% myrawcoef_var1
    V <- OgY %*% myrawcoef_var2

    output$OGStataU <- U
    output$OGStataV <- V

    UVW <- cbind(U, V, cc_object$weights)
    output$OGStataUVW <- as.data.frame(UVW)

    weighted_varU <- 97 * rep(1, ncol(U))
    weighted_varV <- 98 * rep(1, ncol(V))
    weighted_rho <- 99 * rep(1, ncol(U))

    for (i in 1:ncol(U)) {
        meanU <- mean(U[, i])
        meanV <- mean(V[, i])
        U[, i] <- U[, i] - meanU
        V[, i] <- V[, i] - meanV

        # Now the means should be zero
        meanU <- mean(U[, i])
        meanV <- mean(V[, i])
        weighted_varU[i] <- t(U[, i]) %*% diag_W %*% U[, i] / sum(diag_W)
        weighted_varV[i] <- t(V[, i]) %*% diag_W %*% V[, i] / sum(diag_W)
        weighted_rho[i] <- ((t(U[, i]) %*% diag_W %*% V[, i]) / sum(diag_W)) / sqrt(weighted_varU[i] * weighted_varV[i])
        myResults[(7 * (i - 1)) + 1, 1] <- weighted_rho[i] # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 2] <- NA # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 3] <- NA # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 4] <- NA # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 5] <- NA # row 1, 7,...
    }

    if (selection == "FREQ") {
        samplesize <- sum(trunc(cc_object$weights))
        samplesize
    } else {
        samplesize <- nrow(cc_object$X)
        samplesize
    }

    # Lawley's approximation term
    Lawley <- vector(length = length(weighted_rho))
    for (j in 1:length(weighted_rho)) {
        for (i in 1:j) { # Notice Eq. (8) starts in 1, i.e. when d=0 then i=0+1=1
            Lawley[j] <- Lawley[j] + (1 / (weighted_rho[i]^2))
        } # end for i loop
    } # end j for loop
    # End Lawley's approximation term

    Lambda <- vector(length = length(weighted_rho))
    df <- vector(length = length(weighted_rho))
    ChiSq_Wilks_Lambda <- vector(length = length(weighted_rho))
    p_val_Wilks_Lambda <- vector(length = length(weighted_rho))
    # Independence Likelihood ratio test Eq. (3.12) and (3.13) in page 57 when i=1
    # Dimensionality Likelihood ratio test Eq. (3.19) and (3.20) in page 60 when i=2,...,s
    for (i in 1:length(weighted_rho)) {
        Lambda[i] <- (1 - (weighted_rho[i]^2))
        for (j in (i + 1):length(weighted_rho)) {
            Lambda[i] <- Lambda[i] * (1 - (weighted_rho[j]^2)) # Eq. (3.19) of Gittins and Eq. (7) of Calinski
        } # end for j loop
        df[i] <- (ncol(OgX) + 1 - i) * (ncol(OgY) + 1 - i) # Same degrees of freedom for Batlett/Gittins and Lawley
        ChiSq_Wilks_Lambda[i] <- (((samplesize - 1) - 0.5 * (ncol(OgX) + ncol(OgY) + 1) + Lawley[i] - i) * log(Lambda[i])) * -1 # Eq. below Table 3.3 of Gittins and Eq. (12) of Calinski
        p_val_Wilks_Lambda[i] <- pchisq(ChiSq_Wilks_Lambda[i], df[i], lower.tail = FALSE)
        myResults[(7 * (i - 1)) + 2, 1] <- Lambda[i] # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 2] <- df[i] # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 3] <- NA # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 4] <- ChiSq_Wilks_Lambda[i] # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 5] <- p_val_Wilks_Lambda[i] # row 2, 10,...
    } # end for i loop

    V <- rep(0, length(weighted_rho))
    Pillais_Trace_stat <- rep(1, length(weighted_rho))
    Pillais_Trace_p_value <- rep(1, length(weighted_rho))
    for (j in 1:length(weighted_rho)) {
        for (i in j:length(weighted_rho)) {
            V[j] <- V[j] + (weighted_rho[i]^2)
        }
        Pillais_Trace_stat[j] <- (samplesize - 1 - (2 * j) + Lawley[j]) * V[j]
        Pillais_Trace_p_value[j] <- pchisq(Pillais_Trace_stat[j], (ncol(OgX) + 1 - j) * (ncol(OgY) + 1 - j), lower.tail = FALSE)
        myResults[(7 * (j - 1)) + 3, 1] <- V[j]
        myResults[(7 * (j - 1)) + 3, 2] <- (ncol(OgX) + 1 - j) * (ncol(OgY) + 1 - j) - 1
        myResults[(7 * (j - 1)) + 3, 3] <- NA
        myResults[(7 * (j - 1)) + 3, 4] <- Pillais_Trace_stat[j]
        myResults[(7 * (j - 1)) + 3, 5] <- Pillais_Trace_p_value[j]
    }

    U <- rep(0, length(weighted_rho))
    Hotelling_Lawley_Trace_stat <- rep(1, length(weighted_rho))
    Hotelling_Lawley_Trace_p_value <- rep(1, length(weighted_rho))
    for (j in 1:length(weighted_rho)) {
        for (i in j:length(weighted_rho)) {
            U[j] <- U[j] + (weighted_rho[i]^2) / (1 - (weighted_rho[i]^2))
        }
        Hotelling_Lawley_Trace_stat[j] <- (samplesize - ncol(OgX) - ncol(OgY) - 2 + Lawley[j]) * U[j]
        Hotelling_Lawley_Trace_p_value[j] <- pchisq(Hotelling_Lawley_Trace_stat[j], (ncol(OgX) + 1 - j) * (ncol(OgY) + 1 - j), lower.tail = FALSE)
        myResults[(7 * (j - 1)) + 4, 1] <- U[j]
        myResults[(7 * (j - 1)) + 4, 2] <- (ncol(OgX) + 1 - j) * (ncol(OgY) + 1 - j)
        myResults[(7 * (j - 1)) + 4, 3] <- NA
        myResults[(7 * (j - 1)) + 4, 4] <- Hotelling_Lawley_Trace_stat[j]
        myResults[(7 * (j - 1)) + 4, 5] <- Hotelling_Lawley_Trace_p_value[j]
    }

    p <- ncol(OgX)
    q <- ncol(OgY)
    largest_root <- rep(1, length(weighted_rho))
    v1 <- rep(1, length(weighted_rho))
    v2 <- rep(1, length(weighted_rho))
    Roys_Greatest_Root_stat <- rep(1, length(weighted_rho))
    Roys_Greatest_Root_p_value <- rep(1, length(weighted_rho))
    # First we have the so-called "Overall test"

    for (s in 0:(length(weighted_rho) - 1)) {
        par1 <- p
        par2 <- samplesize + s - q - 1
        par3 <- q - s

        fraks <- min(par1, par3)
        frakm <- (abs(par1 - par3) - 1) / 2
        frakn <- (par2 - par1 - 1) / 2

        v1[s + 1] <- fraks + (2 * frakm) + 1
        v2[s + 1] <- fraks + (2 * frakn) + 1
        largest_root[s + 1] <- weighted_rho[1]^2 # Remember that the the cycle starts with r_1 and ends at r_s-1
        # if(v1==0){
        #   v1=.000001
        # }
        Roys_Greatest_Root_stat[s + 1] <- (v2[s + 1] * largest_root[s + 1]) / (v1[s + 1] * (1.0 - largest_root[s + 1])) # We are storing
        Roys_Greatest_Root_p_value[s + 1] <- pf(Roys_Greatest_Root_stat[s + 1], v1[s + 1], v2[s + 1], lower.tail = FALSE)
        myResults[(7 * (s + 1 - 1)) + 5, 1] <- largest_root[s + 1] # row 5, 13,...
        myResults[(7 * (s + 1 - 1)) + 5, 2] <- v1[s + 1] # row 5, 13,...
        myResults[(7 * (s + 1 - 1)) + 5, 3] <- v2[s + 1] # row 5, 13,...
        myResults[(7 * (s + 1 - 1)) + 5, 4] <- Roys_Greatest_Root_stat[s + 1] # row 5, 13,...
        myResults[(7 * (s + 1 - 1)) + 5, 5] <- Roys_Greatest_Root_p_value[s + 1]
    }

    for (i in 1:length(weighted_rho)) { # i will control the rows
        # Row 1: CC, Row 2: Wilk's Lambda, Row 3: Pillai's Trace, Row 4: Hotelling Trace, Row 5: Roy's Greatest Root, Row 6: Wilk's Lambda FREQ, Row 7 Weighted Regression, Row 8: CSD Reg
        # Column 6: Just an index that keeps track of which correlation number each row belongs to
        myResults[(7 * (i - 1)) + 1, 6] <- i # row 1, 9,...
        myResults[(7 * (i - 1)) + 2, 6] <- i # row 2, 10,...
        myResults[(7 * (i - 1)) + 3, 6] <- i # row 3, 11,...
        myResults[(7 * (i - 1)) + 4, 6] <- i # row 4, 12,...
        myResults[(7 * (i - 1)) + 5, 6] <- i # row 5, 13,...
        myResults[(7 * (i - 1)) + 6, 6] <- i # row 6, 14,...
    } # end for loop


    output$Results <- myResults

    return(output)
}


# data <- read.csv("BCC/toy/fit.csv")
# cc <- cancor(cbind(Chins, Situps, Jumps) ~ Weight + Waist + Pulse, data = data, weights = Weight)

# Results <- csdcanon(cc_object = cc, howmany = 2, dim1 = 1, dim2 = 2, selection = "FREQ", cc$ndim, cc$coef$X, cc$coef$Y)

# stargazer(Results,
#     title = "Complete List of Canonical Correlation p-values", type = "text"
#     # column.labels = c("CC Index", "CC Magnitude", "Wilk's Lambda", "Wilk's Lambda (n)", "Roy's Greatest Root", "Pillai's Trace", "Hotelling-Lawley Trace", "Weighted Regression", "CSD Regression"),
#     # align = TRUE, digits = 4, column.sep.width = "-",
#     # row.sep = "-",
#     # float = FALSE
# )
