library(progressr)
handlers(global = TRUE)
handlers("progress", "beepr")

plan(multisession)

data <- c(1, 1, 2, 4, 7, 10)

step_size <- 0.01

psi_grid <- data |> 
  length() |> 
  log() |> 
  round_any(step_size, ceiling) |> 
  seq(0, to = _, step_size)

R <- 250

n_sims <- 1000

sims <- rmultinom(n_sims, length(data), data) |> 
  data.frame() |> 
  as.list() 

test <- 
  mclapply(get_multinom_entropy_IL, 0.1, 10) 

test[[1]] |>
  apply(2, mean) |>
  log() |>
  as.double()

do.call()


my_fcn <- function(xs) {
  p <- function(...) message(...)
  y <- lapply(xs, inner_func, p = p)
}

inner_func <- function(x, p) {
  p(sprintf("x=%g", x))
  sqrt(x)
}

my_fcn(1:5)


xs <- 1:5

with_progress({
  p <- progressor(along = xs)
  y <- future_lapply(xs, function(x, ...) {
    p(sprintf("x=%g", x))
    Sys.sleep(6.0-x)
    sqrt(x)
  })
})

slow_sum <- progressr::slow_sum
print(slow_sum)

x <- 1:10

## Without progress updates
y <- slow_sum(x)

handlers("txtprogressbar")  ## default
with_progress({
  y <- slow_sum(x)
})

if (requireNamespace("progress", quietly = TRUE)) {
  handlers("progress")
  with_progress({
    y <- slow_sum(x)
  })
}

