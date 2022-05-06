library(filling)

# transform 5% of entries into missing
A <- aux.rndmissing(lena128, x = 0.05)

# apply the method with 3 lambda values
fill <- fill.SoftImpute(A, lambdas = c(500, 100, 50))$X
fill2 <- soft_impute(A, lambdas = c(500, 100, 50))

all.equal(fill, fill2)

bench::mark(
  "filling" = fill.SoftImpute(A, lambdas = c(500, 100, 50))$X,
  "our" = soft_impute(A, lambdas = c(500, 100, 50))
)
