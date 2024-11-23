library(tidyverse)

x <- iris |> 
  filter(Species == "setosa") |> 
  pull(Petal.Width)

y <- iris |> 
  filter(Species == "setosa") |> 
  pull(Petal.Length)

seed <- 398724