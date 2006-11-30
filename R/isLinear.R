isLinear <- function(object, ...)
  UseMethod("isLinear")

isLinear.default <- function(object, ...)
  stop("no linearity tests available for this model")

isLinear <- function(object, s_t) {

    sampleSize <- length(object$yy);
    T <- object$n.used - object$m;  # The number of lagged samples

    # Build the regressand vector
    y_t <- object$yy;

    # Build the regressors matrix
    x_t <- object$xx;

    # "1. Regress y_t on x_t and compute the residual sum of squares"
    regression1 <- lm(y_t ~ ., data=data.frame(x_t));
    SSR0 <- sum(regression1$residuals^2);

    # "2. Regress y_t (or regression1$resid) on x_t and x_t * s_t
    #      (first order) and compute the residual sum of squares"

    aux_data1 <- data.frame(y_t = y_t, a = x_t, b = x_t * s_t);
    aux_regression1 <- lm(y_t ~ ., data=aux_data1);
    SSR1 <- sum(aux_regression1$residuals^2);

    # 3. Compute the first order statistic
    n <- object$m + 1;
    m <- dim(aux_data1)[2] - n;
    F_1 <- ((SSR0 - SSR1) / m) / (SSR1 / (T - n - m));

    # Look up the statistic in the table, get the p-value
    lmStatTaylor1 <- pf(F_1, m, T - m - n, lower.tail = FALSE);

    # Regress y_t on the restrictions and compute the RSS
    aux_data2 <- data.frame(y_t = y_t, a = x_t, b = x_t * s_t,
                            c = x_t * s_t^2, d = x_t * s_t^3)
    aux_regression2 <- lm(y_t ~ ., data=aux_data2)
    SSR2 <- sum(aux_regression2$residuals^2);

    # Compute the third order statistic
    n <- object$m + 1;
    m <- dim(aux_data2)[2] - n;
    F_2 = ((SSR0 - SSR2) / m) / (SSR2 / (T - m - n));

    # Look up the statistic in the table, get the p-value
    lmStatTaylor3 <- pf(F_2, m, T - m - n, lower.tail = FALSE);

    c(firstOrderTest = lmStatTaylor1, thirdOrderTest = lmStatTaylor3)

}