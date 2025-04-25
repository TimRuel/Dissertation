# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../utils.R")

seed <- 7835

set.seed(seed)

# True Parameter Generation -----------------------------------------------

J <- 6 # number of levels of response variable

X1_levels <- c("A", "B", "C") # Levels of categorical predictor

p <- ncol(model.matrix( ~ factor(X1_levels)[1] * J - 1)) # number of terms in model after dummy encoding

entropy_ranges <- get_entropy_ranges(J, X1_levels)

X2_intervals <- list("A" = c(0, 1),
                     "B" = c(0, 1),
                     "C" = c(0, 1))

X2_shape_params <- list("A" = c(2, 5),
                        "B" = c(3, 3),
                        "C" = c(5, 2))

num_vals <- 1000

Beta_0 <- get_Beta_0(X1_levels, p, J, X2_intervals, num_vals)

X1_ref_level <- X1_levels[1]

n_samples <- 1e6

X1 <- X1_levels |> 
  rep(each = n_samples) |> 
  factor() |> 
  relevel(X1_ref_level)

X2_arg_list <- list(X2_intervals, X2_shape_params, n_samples) |> 
  modify_if(is.numeric, ~ setNames(rep(., length(X1_levels)), X1_levels)) |> 
  unlist(recursive = FALSE) |> 
  split(X1_levels) |> 
  map(unname) |> 
  map(~ setNames(., c("interval", "shape", "n_samples")))

X2 <- generate_X2_samples(X2_arg_list)

X_design <- model.matrix(~ X1*X2 - 1)

true_probs <- X_design %*% cbind(0, Beta_0) |> 
  apply(1, softmax) |> 
  t() |> 
  data.frame() |> 
  rename_with( ~ paste0("Y", 1:J)) |> 
  mutate(X1 = rep(X1_levels, each = n_samples),
         X2 = X2) |> 
  select(X1, X2, everything())

theta_0 <- true_probs |> 
  select(-X2) |> 
  group_by(X1) |> 
  summarise(across(everything(), \(x) mean(x))) |> 
  data.frame()

H_0 <- theta_0 |> 
  group_by(X1) |> 
  rowwise() |> 
  mutate(entropy = entropy(c_across(everything()))) |> 
  select(X1, entropy) |> 
  data.frame()

base_plot <- ggplot() +   
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#444444", color = NA),  # Dark panel
    plot.background = element_rect(fill = "#2A2A2A", color = NA),  # Dark background
    panel.grid.major = element_line(color = "#4C4C4C"),  # Darker grid
    panel.grid.minor = element_line(color = "#333333"),  
    axis.ticks = element_line(color = "white"),
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    strip.text = element_text(color = "white"),
    plot.title = element_text(color = "white", face = "bold"),
    legend.background = element_rect(fill = "#444444"),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

entropy_df <- true_probs |> 
  mutate(X2 = plyr::round_any(X2, 0.01)) |> 
  group_by(X1) |> 
  distinct(X2, .keep_all = TRUE) |> 
  arrange(X2) |> 
  rowwise() |>
  mutate(entropy = entropy(c_across(-X2))) |> 
  select(X1, X2, entropy) |> 
  data.frame()

base_plot + 
  geom_line(aes(x = X2, y = entropy, color = X1),
            size = 1,
            data = entropy_df) +
  geom_hline(yintercept = unlist(entropy_ranges), color = "black", size = 0.6, linetype = "dashed") +
  labs(title = "True Entropy of Probabilities of Y Classes by Theoretical X1 and X2", x = "X2", y = "Entropy", color = "X1 Level")

# Data Generation ---------------------------------------------------------

m <- c(60, 60, 60)  # number of observations at each level of categorical predictor

names(m) <- X1_levels

C <- length(X1_levels) # number of levels of categorical predictor

n <- sum(m) # total number of observations

X1 <- X1_levels |> 
  rep(times = m) |> 
  factor() |> 
  relevel(X1_ref_level)

X2_arg_list <- list(X2_intervals, X2_shape_params, m) |> 
  unlist(recursive = FALSE) |> 
  split(X1_levels) |> 
  map(unname) |> 
  map(~ setNames(., c("interval", "shape", "n_samples")))

X2 <- generate_X2_samples(X2_arg_list)

X_design <- model.matrix(~ X1*X2 - 1)

Y_probs <- X_design %*% cbind(0, Beta_0) |> 
  apply(1, softmax) |> 
  t() |> 
  data.frame() |> 
  rename_with( ~ paste0("Y", 1:J)) |> 
  mutate(X1 = rep(X1_levels, times = m),
         X2 = X2) |> 
  select(X1, X2, everything())

Y <- Y_probs |>
  select(-c(X1, X2)) |>
  apply(1, \(prob) sample(1:J, size = 1, prob = prob)) |>
  unlist() |>
  unname() |>
  factor(levels = 1:J)

Y_design <- model.matrix(~ Y)[,-1]

data <- data.frame(X1 = X1,
                   X2 = X2,
                   Y)

formula <- Y ~ .^2 - 1

model <- fit_multinomial_logistic_model(data, formula)

Beta_MLE <- get_Beta_MLE(model)

observed_entropy_df <- Y_probs |> 
  group_by(X1) |> 
  rowwise() |> 
  mutate(entropy = entropy(c_across(-X2))) |> 
  select(X1, X2, entropy) |> 
  data.frame()

Y_probs_long <- Y_probs |> 
  pivot_longer(cols = paste0("Y", 1:J),
               names_to = "Y_class",
               names_prefix = "Y",
               values_to = "Probability")

base_plot + 
  geom_point(aes(x = X2, y = Probability, color = Y_class),
             size = 1,
             data = Y_probs_long) +
  facet_wrap(~X1, 
             ncol = 1,
             labeller = labeller(X1 = function(x) paste("X1 =", x)),
             scales = "free_x") +
  scale_x_continuous(breaks = seq(0, 4, 0.5), expand = c(0.05, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.1, 0)) +
  labs(title = "True Probabilities of Y Classes by X1 and X2", 
       x = "X2", 
       y = "Probability", 
       color = "Class of Y")

base_plot + 
  geom_line(aes(x = X2, y = entropy, color = X1),
            size = 1,
            data = observed_entropy_df) +
  geom_hline(yintercept = unlist(entropy_ranges), color = "black", size = 0.6, linetype = "dashed") +
  labs(title = "True Entropy of Probabilities of Y Classes by Observed X1 and X2", x = "X2", y = "Entropy", color = "X1 Level")

base_plot + 
  geom_histogram(aes(x = X2, fill = Y),
                 color = "black",
                 binwidth = 0.1, 
                 position = "stack",
                 data = data) +
  facet_wrap(~X1, 
             ncol = 1,
             labeller = labeller(X1 = function(x) paste("X1 =", x)),
             scales = "free_x") +
  scale_x_continuous(breaks = seq(-2, 4, 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Observed Y Values by X1 and X2", x = "X2", y = "Count", fill = "Observed Y")
