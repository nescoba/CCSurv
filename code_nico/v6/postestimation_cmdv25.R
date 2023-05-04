#' Description
#'
#' @param thing
#' @param other
#'
#' @return thing2

csdcanon <- function(data, varlist1, varlist2, finwgt, n, stdcoef_var1, stdcoef_var2, dim1 = NA, selection, svysetcommand = "classic", dim2 = NA) {
    varcount1 <- length(varlist1)
    varcount2 <- length(varlist2)

    # create a matrix with the columns of 'data' that are listed in 'varlist1'
    auxOgX <- data[, varlist1]
    # create a matrix with the columns of 'data' that are listed in 'varlist2'
    auxOgY <- data[, varlist2]
    # Capturing all of the variables names
    names <- c(varlist1, varlist2)

    # First set of variables is the longest
    totalvar <- 0
    if (varcount1 >= varcount2) {
        X <- varlist1
        Y <- varlist2
        totalvar <- varcount2 + 1
    } else {
        X <- varlist2
        Y <- varlist1
        totalvar <- varcount1 + 1
    }

    OgX <- data[, X]
    OgY <- data[, Y]

    n <- varcount1 + varcount2

    varname <- vector(length = n)

    # Create a new variable called counter and assign it the observation number of each row in the dataset
    counter <- 1:n

    # Iterate over a range of values from 1 to a specified number n
    for (i in 1:n) {
        # Get the ith name
        v1 <- names[i]
        # Assign name i to value i of the variable varname
        varname[counter == i] <- v1
    }

    if (!is.null(dim1)) {
        # Notice that the order of the std coeff is still the same as in canon hence the variables names are in the correct order

        mycoef1 <- stdcoef_var1
        stdU <- auxOgX %*% mycoef1
        stdU <- as.matrix(stdU, ncol = 1)
        mycoef1 <- cbind(mycoef1, matrix(0, nrow = nrow(mycoef1)))

        mycoef2 <- stdcoef_var2
        stdV <- auxOgY %*% mycoef2
        stdV <- as.matrix(stdV, ncol = 1)
        mycoef2 <- cbind(mycoef2, matrix(1, nrow = nrow(mycoef2)))

        # Variables

        graphcoef <- cbind(mycoef1, mycoef2)

        # Needs to add different color and shape to the markers depending on which set of variables they came from

        plot(graphcoef[, dim2], graphcoef[, dim1], col = ifelse(graphcoef[, totalvar] == 0, "red", "blue"), pch = ifelse(graphcoef[, totalvar] == 0, 22, 4), xlab = paste0("CC", dim1), ylab = paste0("CC", dim2), xlim = c(-1, 1), ylim = c(-1, 1))

        # Units

        stdV_dim1 <- stdV[dim1]
        stdU_dim1 <- stdU[dim1]

        plot(stdU_dim1, stdV_dim1, xlab = paste0("U", dim1), ylab = paste0("V", dim1))

        # Combine

        dev.copy(png, file = "variables.png")
        dev.off()

        dev.copy(png, file = "units.png")
        dev.off()
    } # end if graphs


    if (svysetcommand != "classic") {
        W <- finwgt
        # Creating diagonal matrix with the survey Weights
        diag_W <- diag(nrow(OgX))
        for (i in 1:nrow(OgX)) {
            diag_W[i, i] <- W[i, 1]
        }
    } else {
        # Creating a weight matrix with all weights equal to 1 for the classic case
        W <- matrix(1, nrow = nrow(OgX), ncol = 1)
        diag_W <- diag(nrow(OgX))
        for (i in 1:nrow(OgX)) {
            diag_W[i, i] <- 1
        }
    }

    # Standarizing X variables
    OgX <- apply(OgX, 2, scale)
    X <- as.matrix(OgX)

    # Standarizing Y variables
    OgY <- apply(OgY, 2, scale)
    Y <- as.matrix(OgY)

    # mkmat sends the matrices to mata so there is no need to send them as arguments, e(rawcoef_var1) and e(rawcoef_var1) are sent to mata automatically when they are created by canon
    # Capturing number of canonical correlations, e.g. e(n_cc) = 3 has the number of canonical correlations; it is a scalar
    n_cc <- e(n_cc)

    # Creating matrix where all new p-values will be stored
    Results <- matrix(NA, nrow = 7 * n_cc, ncol = 6)

    output <- calcpval(selection)

    Results <- output$Results
    OGStataUVW <- output$OGStataUVW

    colnames(Results) <- c("Statistic", "df1", "df2", "Chi-Sq/F", "p-val", "Index")

    weightindex <- (2 * n_cc) + 1 # Reading the weights as they were stored in this column in the mata function

    if (svysetcommand != "classic") {
        svyset(pweight = paste0("OGStataUVW", weightindex))

        # Notice coefficients of the simple linear regression are equal to the canonical correlations because the variances of the canonical variates are equal to 1
        for (i in 1:n_cc) {
            secondindex <- i + n_cc
            svyreg <- svyglm(as.formula(paste0("OGStataUVW", i, " ~ OGStataUVW", secondindex)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
            p1 <- coef(summary(svyreg))[1, 4]

            svyreg <- svyglm(as.formula(paste0("OGStataUVW", secondindex, " ~ OGStataUVW", i)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
            p2 <- coef(summary(svyreg))[1, 4]

            if (p1 >= p2) {
                svyreg <- svyglm(as.formula(paste0("OGStataUVW", i, " ~ OGStataUVW", secondindex)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
                Results[(7 * (i - 1)) + 6, 1] <- coef(summary(svyreg))[1, 1]
                Results[(7 * (i - 1)) + 6, 2] <- svyreg$df.residual + length(svyreg$coefficients) - 1
                Results[(7 * (i - 1)) + 6, 3] <- svyreg$df.null - svyreg$df.residual
                Results[(7 * (i - 1)) + 6, 4] <- summary(svyreg)$fstatistic[1]
                Results[(7 * (i - 1)) + 6, 5] <- summary(svyreg)$coef[2, 4]
                Results[(7 * (i - 1)) + 6, 6] <- i
            } else {
                svyreg <- svyglm(as.formula(paste0("OGStataUVW", secondindex, " ~ OGStataUVW", i)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
                Results[(7 * (i - 1)) + 6, 1] <- coef(summary(svyreg))[1, 1]
                Results[(7 * (i - 1)) + 6, 2] <- svyreg$df.residual + length(svyreg$coefficients) - 1
                Results[(7 * (i - 1)) + 6, 3] <- svyreg$df.null - svyreg$df.residual
                Results[(7 * (i - 1)) + 6, 4] <- summary(svyreg)$fstatistic[1]
                Results[(7 * (i - 1)) + 6, 5] <- summary(svyreg)$coef[2, 4]
                Results[(7 * (i - 1)) + 6, 6] <- i
            }
        }
    }

    if (svysetcommand != "classic") {
        svysetcommand

        for (i in 1:n_cc) {
            secondindex <- i + n_cc
            svyreg <- svyglm(as.formula(paste0("OGStataUVW", i, " ~ OGStataUVW", secondindex)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
            p1 <- coef(summary(svyreg))[2, 4]

            svyreg <- svyglm(as.formula(paste0("OGStataUVW", secondindex, " ~ OGStataUVW", i)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
            p2 <- coef(summary(svyreg))[2, 4]

            if (p1 >= p2) {
                svyreg <- svyglm(as.formula(paste0("OGStataUVW", i, " ~ OGStataUVW", secondindex)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
                Results[(7 * (i - 1)) + 7, 1] <- coef(summary(svyreg))[1, 1]
                Results[(7 * (i - 1)) + 7, 2] <- svyreg$df.residual + length(svyreg$coefficients) - 1
                Results[(7 * (i - 1)) + 7, 3] <- svyreg$df.null - svyreg$df.residual
                Results[(7 * (i - 1)) + 7, 4] <- summary(svyreg)$fstatistic[1]
                Results[(7 * (i - 1)) + 7, 5] <- summary(svyreg)$coef[2, 4]
                Results[(7 * (i - 1)) + 7, 6] <- i
            } else {
                svyreg <- svyglm(as.formula(paste0("OGStataUVW", secondindex, " ~ OGStataUVW", i)), design = svydesign(ids = ~1, weights = ~ OGStataUVW[[paste0("OGStataUVW", weightindex)]]))
                Results[(7 * (i - 1)) + 7, 1] <- coef(summary(svyreg))[1, 1]
                Results[(7 * (i - 1)) + 7, 2] <- svyreg$df.residual + length(svyreg$coefficients) - 1
                Results[(7 * (i - 1)) + 7, 3] <- svyreg$df.null - svyreg$df.residual
                Results[(7 * (i - 1)) + 7, 4] <- summary(svyreg)$fstatistic[1]
                Results[(7 * (i - 1)) + 7, 5] <- summary(svyreg)$coef[2, 4]
                Results[(7 * (i - 1)) + 7, 6] <- i
            }
        }
    }

    if (howmany == "") {
        howmany <- n_cc
    }

    for (i in 1:howmany) {
        # Selecting rows in information for Canonical correlation i
        aux <- Results[(7 * (i - 1) + 1), 1]
        mag <- round(aux[1, 1], 0.0001)
        ResultsAux <- Results[(7 * (i - 1) + 2):(7 * (i - 1) + 7), 1:5]
        rownames(ResultsAux) <- c("Wilks' Lambda", "Pillai's Trace", "Hotelling-Lawley Trace", "Roy's Greatest Root", "Weighted Reg", "Complex Survey Design Reg")
        print(ResultsAux, title = paste0("Statistics for Canonical Correlation: ", i), twidth = 30, format = "%10.4f", rowtitle = paste0("Canonical Correlation=", mag), border = c("top", "bottom"))
    }
}

calcpval <- function(selection) {
    if (nrow(myrawcoef_var2) > nrow(myrawcoef_var1)) {
        aux <- myrawcoef_var2
        myrawcoef_var2 <- myrawcoef_var1
        myrawcoef_var1 <- aux
    }

    U <- myOgX %*% myrawcoef_var1
    V <- myOgY %*% myrawcoef_var2
    # Now U and V are not complex matrices

    # Sending U and V from Mata to R, in R they will be called OGStataU and OGStataV
    OGStataU <- as.matrix(U)
    OGStataV <- as.matrix(V)
    UV <- cbind(U, V)
    UVW <- cbind(UV, myW)
    OGStataUVW <- as.matrix(UVW)

    # Calculating Canonical Correlations
    # Using Eq. (8) and (9) to find the correlations and variances of the canonical matrices
    # This calculation does not lead to correct canonical correlations or variances of the canonical variates as it did when using all my process above
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
        weighted_varU[i] <- t(U[, i]) %*% mydiag_W %*% U[, i] / sum(mydiag_W)
        weighted_varV[i] <- t(V[, i]) %*% mydiag_W %*% V[, i] / sum(mydiag_W)
        weighted_rho[i] <- ((t(U[, i]) %*% mydiag_W %*% V[, i]) / sum(mydiag_W)) / sqrt(weighted_varU[i] * weighted_varV[i])
        myResults[(7 * (i - 1)) + 1, 1] <- weighted_rho[i] # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 2] <- NA # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 3] <- NA # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 4] <- NA # row 1, 7,...
        myResults[(7 * (i - 1)) + 1, 5] <- NA # row 1, 7,...
    }

    if (selection == "FREQ") {
        samplesize <- sum(trunc(myW))
        samplesize
    } else {
        samplesize <- nrow(myX)
        samplesize
    }
    # End Finding sample size
    # Lawley's approximation term
    Lawley <- matrix(0, nrow = 1, ncol = ncol(weighted_rho))
    for (j in 1:ncol(weighted_rho)) {
        for (i in 1:j) { # Notice Eq. (8) starts in 1, i.e. when d=0 then i=0+1=1
            Lawley[1, j] <- Lawley[1, j] + (1 / (weighted_rho[1, i]^2))
        } # end for i loop
    } # end j for loop
    # End Lawley's approximation term

    Lambda <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    df <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    ChiSq_Wilks_Lambda <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    p_val_Wilks_Lambda <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    # Independence Likelihood ratio test Eq. (3.12) and (3.13) in page 57 when i=1
    # Dimensionality Likelihood ratio test Eq. (3.19) and (3.20) in page 60 when i=2,...,s
    for (i in 1:ncol(weighted_rho)) {
        Lambda[1, i] <- (1 - (weighted_rho[1, i]^2))
        for (j in (i + 1):ncol(weighted_rho)) {
            Lambda[1, i] <- Lambda[1, i] * (1 - (weighted_rho[1, j]^2)) # Eq. (3.19) of Gittins and Eq. (7) of Calinski
        } # end for j loop
        df[1, i] <- (ncol(myX) + 1 - i) * (ncol(myY) + 1 - i) # Same degrees of freedom for Batlett/Gittins and Lawley
        ChiSq_Wilks_Lambda[1, i] <- (((samplesize - 1) - 0.5 * (ncol(myX) + ncol(myY) + 1) + Lawley[1, i] - i) * log(Lambda[1, i])) * -1 # Eq. below Table 3.3 of Gittins and Eq. (12) of Calinski
        p_val_Wilks_Lambda[1, i] <- pchisq(ChiSq_Wilks_Lambda[1, i], df[1, i], lower.tail = FALSE)
        myResults[(7 * (i - 1)) + 2, 1] <- Lambda[1, i] # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 2] <- df[1, i] # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 3] <- NA # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 4] <- ChiSq_Wilks_Lambda[1, i] # row 2, 10,...
        myResults[(7 * (i - 1)) + 2, 5] <- p_val_Wilks_Lambda[1, i] # row 2, 10,...
    } # end for i loop

    V <- matrix(0, nrow = 1, ncol = ncol(weighted_rho))
    Pillais_Trace_stat <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    Pillais_Trace_p_value <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    for (j in 1:ncol(weighted_rho)) {
        for (i in j:ncol(weighted_rho)) { # Notice Eq. (8) starts in 1, i.e. when d=0 then i=0+1=1
            V[1, j] <- V[1, j] + (weighted_rho[1, i]^2) # Eq(8) of Calinski
        } # end for i loop
        # The degrees of freedom are the same as in SAS 9 & 48, in Eq (8) of of All formulas p-values.pdf
        # Pillais_Trace_stat[1,j] = (samplesize - 1 * V[1,j] # Eq. (11) of Calinski
        Pillais_Trace_stat[1, j] <- (samplesize - 1 - (2 * j) + Lawley[1, j]) * V[1, j] # Eq. (14) of Calinski
        Pillais_Trace_p_value[1, j] <- pchisq(Pillais_Trace_stat[1, j], (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j), lower.tail = FALSE)
        myResults[(7 * (j - 1)) + 3, 1] <- V[1, j] # row 3, 11,...
        myResults[(7 * (j - 1)) + 3, 2] <- (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j) - 1 # Notice the -1 to transform the i into d # row 3, 11,...
        myResults[(7 * (j - 1)) + 3, 3] <- NA # row 3, 11,...
        myResults[(7 * (j - 1)) + 3, 4] <- Pillais_Trace_stat[1, j] # row 3, 11,...
        myResults[(7 * (j - 1)) + 3, 5] <- Pillais_Trace_p_value[1, j] # row 3, 11,...
    } # end for j loop

    U <- matrix(0, nrow = 1, ncol = ncol(weighted_rho))
    Hotelling_Lawley_Trace_stat <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    Hotelling_Lawley_Trace_p_value <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    for (j in 1:ncol(weighted_rho)) {
        for (i in j:ncol(weighted_rho)) {
            U[1, j] <- U[1, j] + (weighted_rho[1, i]^2) / (1 - (weighted_rho[1, i]^2)) # Eq(6) of Calinski
        } # end for i loop
        # Hotelling_Lawley_Trace_stat[1,j] = (samplesize - ncol(myX) - ncol(myY) - 2) * U[1,j] # Eq. (9) of Calinski Bartlett's approximation
        Hotelling_Lawley_Trace_stat[1, j] <- (samplesize - ncol(myX) - ncol(myY) - 2 + Lawley[1, j]) * U[1, j] # Eq. (13) of Calinski Fujikoshi approximation
        Hotelling_Lawley_Trace_p_value[1, j] <- pchisq(Hotelling_Lawley_Trace_stat[1, j], (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j), lower.tail = FALSE)
        myResults[(7 * (j - 1)) + 4, 1] <- U[1, j] # row 4, 12,...
        myResults[(7 * (j - 1)) + 4, 2] <- (ncol(myX) + 1 - j) * (ncol(myY) + 1 - j) # row 4, 12,...
        myResults[(7 * (j - 1)) + 4, 3] <- NA # row 4, 12,...
        myResults[(7 * (j - 1)) + 4, 4] <- Hotelling_Lawley_Trace_stat[1, j] # row 4, 12,...
        myResults[(7 * (j - 1)) + 4, 5] <- Hotelling_Lawley_Trace_p_value[1, j] # row 4, 12,...
    } # end for j loop

    p <- ncol(myX)
    q <- ncol(myY)
    pq <- matrix(c(p, q), nrow = 1, ncol = 2)
    largest_root <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    v1 <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    v2 <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    Roys_Greatest_Root_stat <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    Roys_Greatest_Root_p_value <- matrix(1, nrow = 1, ncol = ncol(weighted_rho))
    # First we have the so-called "Overall test"
    s <- min(pq)
    m <- (abs(p - q) - 1) / 2
    n <- (samplesize - p - q - 2) / 2
    v1[1, 1] <- s + (2 * m) + 1
    v2[1, 1] <- s + (2 * n) + 1
    largest_root[1, 1] <- weighted_rho[1, 1]^2
    Roys_Greatest_Root_stat[1, 1] <- (v2[1, 1] * largest_root[1, 1]) / (v1[1, 1] * (1 - largest_root[1, 1]))
    Roys_Greatest_Root_p_value[1, 1] <- pf(Roys_Greatest_Root_stat[1, 1], v1[1, 1], v2[1, 1], lower.tail = FALSE)

    myResults[5, 1] <- largest_root[1, 1] # row 5, 13,...
    myResults[5, 2] <- v1[1, 1] # row 5, 13,...
    myResults[5, 3] <- v2[1, 1] # row 5, 13,...
    myResults[5, 4] <- Roys_Greatest_Root_stat[1, 1] # row 5, 13,...
    myResults[5, 5] <- Roys_Greatest_Root_p_value[1, 1]

    for (k in 1:(cols(weighted_rho) - 1)) {
        # formula at the bottom of page 61 of Gittins
        s <- min(pq) - k - 1
        # formulas to match Eq(3.17) when k=1
        # s <- min(pq) - k + 1
        # formula at the bottom of page 61 of Gittins
        m <- (abs(p - q - k) - 1) / 2
        # formulas to match Eq(3.17) when k=1
        # m <- (abs(p-q-k+1)-1)/2
        # formula at the bottom of page 61 of Gittins
        n <- (samplesize - p - q - k) / 2
        # formulas to match Eq(3.17) when k=1
        # n <- (samplesize-p-q-(k+1))/2
        v1[1, k + 1] <- s + (2 * m) + 1
        v2[1, k + 1] <- s + (2 * n) + 1
        largest_root[1, k + 1] <- weighted_rho[1, k]^2 # Remember that the the cycle starts with r_1 and ends at r_s-1
        # if(v1==0){
        #   v1=.000001
        # }
        Roys_Greatest_Root_stat[1, k + 1] <- (v2[1, k + 1] * largest_root[1, k + 1]) / (v1[1, k + 1] * (1.0 - largest_root[1, k + 1])) # We are storing the results in k+1 or we will delete the independence union-intersection test
        Roys_Greatest_Root_p_value[1, k + 1] <- Ftail(v1[1, k + 1], v2[1, k + 1], Roys_Greatest_Root_stat[1, k + 1])
        myResults[(7 * (k + 1 - 1)) + 5, 1] <- largest_root[1, k + 1] # row 5, 13,...
        myResults[(7 * (k + 1 - 1)) + 5, 2] <- v1[1, k + 1] # row 5, 13,...
        myResults[(7 * (k + 1 - 1)) + 5, 3] <- v2[1, k + 1] # row 5, 13,...
        myResults[(7 * (k + 1 - 1)) + 5, 4] <- Roys_Greatest_Root_stat[1, k + 1] # row 5, 13,...
        myResults[(7 * (k + 1 - 1)) + 5, 5] <- Roys_Greatest_Root_p_value[1, k + 1] # row 5, 13,...
    } # end for loop

    for (i in 1:cols(weighted_rho)) { # i will control the rows
        # Row 1: CC, Row 2: Wilk's Lambda, Row 3: Pillai's Trace, Row 4: Hotelling Trace, Row 5: Roy's Greatest Root, Row 6: Wilk's Lambda FREQ, Row 7 Weighted Regression, Row 8: CSD Reg
        # Column 6: Just an index that keeps track of which correlation number each row belongs to
        myResults[(7 * (i - 1)) + 1, 6] <- i # row 1, 9,...
        myResults[(7 * (i - 1)) + 2, 6] <- i # row 2, 10,...
        myResults[(7 * (i - 1)) + 3, 6] <- i # row 3, 11,...
        myResults[(7 * (i - 1)) + 4, 6] <- i # row 4, 12,...
        myResults[(7 * (i - 1)) + 5, 6] <- i # row 5, 13,...
        myResults[(7 * (i - 1)) + 6, 6] <- i # row 6, 14,...
    } # end for loop

    return myResults
}
