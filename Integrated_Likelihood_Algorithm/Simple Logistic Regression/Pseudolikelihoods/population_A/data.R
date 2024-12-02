seed <- 3409836

set.seed(seed)

alpha_0 <- -1

beta_0 <- 1

n <- 100

x <- rnorm(n, 0, 1)

y <- rbinom(n, 1, sigmoid::sigmoid(alpha_0 + beta_0 * x))

fit <- glm(y ~ x, family = "binomial")

fit$fitted.values

fit$coefficients

plot(x, y)
points(sort(x), type = "l", sigmoid(fit$coefficients[1] + fit$coefficients[2] * sort(x)))
