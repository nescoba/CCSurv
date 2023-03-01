# load data
fit <- read.csv("BCC/toy/fit.csv")

comp <- cancor(fit[, c("Weight", "Waist", "Pulse")], fit[, c("Chins", "Situps", "Jumps")])$cor

# center and scale the data
library(car)
fit_scaled <- data.frame(cbind(
    scale(fit[, c("Weight", "Waist", "Pulse")]),
    scale(fit[, c("Chins", "Situps", "Jumps")])
))

# create matrices X and Y
X <- as.matrix(fit_scaled[, c("Weight", "Waist", "Pulse")])
Y <- as.matrix(fit_scaled[, c("Chins", "Situps", "Jumps")])

# create matrices OgX and OgY
OgX <- as.matrix(fit[, c("Weight", "Waist", "Pulse")])
OgY <- as.matrix(fit[, c("Chins", "Situps", "Jumps")])


# create diagonal matrix which represents survey weights, but for now we will use the identity matrix
diag_W <- diag(nrow(X))



factorYX <- t(Y) %*% diag_W %*% X
factorInvXX <- solve(t(X) %*% diag_W %*% X)
factorXY <- t(X) %*% diag_W %*% Y
factorInvYY <- solve(t(Y) %*% diag_W %*% Y)

A <- factorYX %*% factorInvXX %*% factorXY %*% factorInvYY

# compute eigenvectors and eigenvalues
library(Matrix)
eigen_A <- eigen(A)
bT <- eigen_A$vectors
lambdaSq <- eigen_A$values
lambda <- sqrt(lambdaSq)
aT <- bT %*% factorYX %*% factorInvXX

# In question
for (i in 1:length(lambda)) {
    aT[i, ] <- (1 / lambda[i]) * aT[i, ]
}

# compute canonical variates
UT <- aT %*% t(X)
U <- t(UT)
VT <- bT %*% t(Y)
V <- t(VT)

# compute canonical correlations and variances
weighted_varU <- numeric(length = length(lambda))
weighted_varV <- numeric(length = length(lambda))
weighted_rho <- numeric(length = length(lambda))
for (i in 1:length(lambda)) {
    weighted_varU[i] <- UT[i, ] %*% diag_W %*% U[, i] / sum(diag_W)
    weighted_varV[i] <- VT[i, ] %*% diag_W %*% V[, i] / sum(diag_W)
    weighted_rho[i] <- (UT[i, ] %*% diag_W %*% V[, i] / sum(diag_W)) / sqrt(weighted_varU[i] * weighted_varV[i])
}



weighted_varU
weighted_varV
weighted_rho










# Now using survey weights

# Diagonal matrix whose entries is the Weights column of fit
diag_W <- diag(fit$Weight)

factorYX <- t(Y) %*% diag_W %*% X
factorInvXX <- solve(t(X) %*% diag_W %*% X)
factorXY <- t(X) %*% diag_W %*% Y
factorInvYY <- solve(t(Y) %*% diag_W %*% Y)

A <- factorYX %*% factorInvXX %*% factorXY %*% factorInvYY
eigen_result <- eigen(A)

bT <- eigen_result$vectors
lambdaSq <- eigen_result$values

lambda <- sqrt(lambdaSq)

aT <- bT %*% factorYX %*% factorInvXX

# In question
for (i in 1:length(lambda)) {
    aT[, i] <- (1 / lambda[i]) * aT[, i]
}

UT <- aT %*% t(X)
U <- t(UT)
VT <- bT %*% t(Y)
V <- t(VT)

weighted_varU <- numeric(length = length(lambda))
weighted_varV <- numeric(length = length(lambda))
weighted_rho <- numeric(length = length(lambda))

for (i in 1:length(lambda)) {
    weighted_varU[i] <- UT[i, ] %*% diag_W %*% U[, i] / sum(diag_W)
    weighted_varV[i] <- VT[i, ] %*% diag_W %*% V[, i] / sum(diag_W)
    weighted_rho[i] <- (UT[i, ] %*% diag_W %*% V[, i] / sum(diag_W)) / sqrt(weighted_varU[i] * weighted_varV[i])
}

# Output
StataU <- U
StataV <- V

data <- cbind(StataU, StataV)
# data as data.frame with columns StataU1, StataU2, StataU3, StataV1, StataV2, StataV3
dat <- data.frame(StataU1 = StataU[, 1], StataU2 = StataU[, 2], StataU3 = StataU[, 3], StataV1 = StataV[, 1], StataV2 = StataV[, 2], StataV3 = StataV[, 3])

# Create a weighted correlation matrix
cor_mat <- cov.wt(data, cor = TRUE, wt = fit$Weight)

# Create a design object with the sampling weights
library(survey)
des <- svydesign(
    ids = ~1,
    data = dat,
    weights = ~ fit$Weight
)

# Run weighted regression on the first canonical variable
fit_u1_v1 <- svyglm(StataU1 ~ StataV1, design = des)
coef_u1_v1 <- coef(summary(fit_u1_v1))
rawcoef_var1 <- coef_u1_v1[1:2, ]

# Run weighted regression on the second canonical variable
fit_u2_v2 <- svyglm(StataU2 ~ StataV2, design = des)
coef_u2_v2 <- coef(summary(fit_u2_v2))
rawcoef_var2 <- coef_u2_v2[1:2, ]

# Run weighted regression on the third canonical variable
fit_u3_v3 <- svyglm(StataU3 ~ StataV3, design = des)
coef_u3_v3 <- coef(summary(fit_u3_v3))
rawcoef_var3 <- coef_u3_v3[1:2, ]

# Calculate the canonical variates in R
canon_variates <- cancor(fit[, c("Weight", "Waist", "Pulse")], fit[, c("Chins", "Situps", "Jumps")], weights = fit$Weight)$cor

canon_variates$cor

U <- X %*% rawcoef_var1
V <- Y %*% rawcoef_var2

OGStataU <- U
OGStataV <- V

weighted_varU <- numeric(length = ncol(U))
weighted_varV <- numeric(length = ncol(V))
weighted_rho <- numeric(length = ncol(U))

# Substract the column mean from each column of U and V
U <- sweep(U, 2, colMeans(U), "-")
V <- sweep(V, 2, colMeans(V), "-")

for (i in 1:ncol(U)) {
    weighted_varU[i] <- t(U)[i, ] %*% diag_W %*% U[, i] / sum(diag_W)
    weighted_varV[i] <- t(V)[i, ] %*% diag_W %*% V[, i] / sum(diag_W)
    weighted_rho[i] <- (t(U)[i, ] %*% diag_W %*% V[, i] / sum(diag_W)) / sqrt(weighted_varU[i] * weighted_varV[i])
}

OGSdata <- data.frame(OGStataU1 = OGStataU[, 1], OGStataU2 = OGStataU[, 2], OGStataU3 = OGStataU[, 3], OGStataV1 = OGStataV[, 1], OGStataV2 = OGStataV[, 2], OGStataV3 = OGStataV[, 3])

des <- svydesign(
    ids = ~1,
    data = OGSdata,
    weights = ~ fit$Weight
)

fit_u1_v1 <- svyglm(OGStataU1 ~ OGStataV1, design = des)
coef_u1_v1 <- coef(summary(fit_u1_v1))
rawcoef_var1 <- coef_u1_v1[1:2, ]

fit_u2_v2 <- svyglm(OGStataU2 ~ OGStataV2, design = des)
coef_u2_v2 <- coef(summary(fit_u2_v2))
rawcoef_var2 <- coef_u2_v2[1:2, ]

fit_u3_v3 <- svyglm(OGStataU3 ~ OGStataV3, design = des)
coef_u3_v3 <- coef(summary(fit_u3_v3))
rawcoef_var3 <- coef_u3_v3[1:2, ]




# Load required packages
library(MASS)

Lambda <- numeric(length = length(weighted_rho))
df <- numeric(length = length(weighted_rho))
ChiSq_Wilks_Lambda <- numeric(length = length(weighted_rho))
p_val_Wilks_Lambda <- numeric(length = length(weighted_rho))

for (i in 1:length(weighted_rho)) {
    Lambda[i] <- 1 - weighted_rho[i]^2
    for (j in (i + 1):length(weighted_rho)) {
        Lambda[i] <- Lambda[i] * (1 - weighted_rho[j]^2)
    }
    df[i] <- (ncol(X) + 1 - i) * (ncol(Y) + 1 - i)
    ChiSq_Wilks_Lambda[i] <- (((nrow(X) - 1) - .5 * (ncol(X) + ncol(Y) + 1)) * log(Lambda[i])) * (-1)
    p_val_Wilks_Lambda[i] <- pchisq(ChiSq_Wilks_Lambda[i], df[i], lower.tail = FALSE)
}

p <- ncol(X)
q <- ncol(Y)
pq <- c(p, q)
s <- min(pq)
m <- (abs(p - q) - 1) / 2
n <- (nrow(X) - p - q - 2) / 2
v1 <- s + 2 * m + 1
v2 <- s + 2 * n + 1


# Compute Roy's greatest root statistic and p-value
largest_root <- weighted_rho[1, 1]^2
Roys_Greatest_Root_stat <- numeric(length = length(weighted_rho))
Roys_Greatest_Root_p_value <- numeric(length = length(weighted_rho))
Roys_Greatest_Root_stat[1] <- (v2 * largest_root) / (v1 * (1.0 - largest_root))
Roys_Greatest_Root_p_value[2] <- pf(Roys_Greatest_Root_stat, df1 = v1, df2 = v2, lower.tail = FALSE)

# Compute the remaining canonical correlations
for (j in 2:ncol(weighted_rho)) {
    # Compute the jth canonical correlation
    s <- min(pq) - j - 1
    m <- (abs(p - q - j) - 1) / 2
    n <- (nrow(X) - p - q - j) / 2
    v1 <- s + 2 * m + 1
    v2 <- s + 2 * n + 1
    largest_root <- weighted_rho[j]^2
    Roys_Greatest_Root_stat[j] <- (v2 * largest_root) / (v1 * (1.0 - largest_root))
    Roys_Greatest_Root_p_value[j] <- pf(Roys_Greatest_Root_stat, df1 = v1, df2 = v2, lower.tail = FALSE)
}

V <- 0
for (i in 1:length(weighted_rho)) {
    V <- V + (weighted_rho[i]^2)
}

s <- min(p, q)
m <- (abs(p - q) - 1) / 2
n <- (nrow(XX) - p - q - 2) / 2

df1 <- s * ((2 * m) + s + 1)
df2 <- s * ((2 * n) + s + 1)

Pillais_Trace_stat <- (((2 * n) + s + 1) / ((2 * m) + s + 1)) * V / (s - V)
Pillais_Trace_p_value <- pf(Pillais_Trace_stat, df1, df2, lower.tail = FALSE)

U <- 0
for (i in 1:ncol(weighted_rho)) {
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
