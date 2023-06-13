# Simulation Dataset
S_true <- matrix(runif(5*5), nrow=5, ncol=5)
S <- matrix(runif(5*5), nrow=5, ncol=5)

out <- CorrIndex(cor(S_true, S))

expect_true(is.numeric(out))
