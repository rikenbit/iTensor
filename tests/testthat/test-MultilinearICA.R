# Simulation Dataset
arrX <- array(runif(10*20*30), dim=c(10,20,30))
X <- as.tensor(arrX)

## Perform MICA against Simulation Datasets
Js <- c(2,3,4)
out <- MultilinearICA(X, Js=Js)

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(out), TRUE)

### Test O-2: Object Names
expect_identical(names(out),
    c("As", "S", "Js", "algorithm"))

### Test 0-3: As
expect_identical(is.list(out$As), TRUE)

expect_equal(dim(out$As[[1]]), c(Js[1], dim(X)[1]))
expect_equal(dim(out$As[[2]]), c(Js[2], dim(X)[2]))
expect_equal(dim(out$As[[3]]), c(Js[3], dim(X)[3]))

### Test 0-4: S
expect_identical(is(out$S)[1], "Tensor")

expect_equal(dim(out$S)[1], Js[1])
expect_equal(dim(out$S)[2], Js[2])
expect_equal(dim(out$S)[3], Js[3])

### Test 0-5: Js
expect_identical(out$Js, Js)

### Test 0-6: algorithm
expect_identical(out$algorithm, "FastICA")

## Test Error
### Test E-1: X
expect_error(MultilinearICA(arrX, Js=Js))

### Test E-2: Js
expect_error(MultilinearICA(X, Js=c("2", "3", "4")))
expect_error(MultilinearICA(X, Js=c(2, 3, 4, 5)))
expect_error(MultilinearICA(X, Js=c(2, 3)))

### Test E-3: modes
expect_error(MultilinearICA(X, modes=c("1", "2", "3")))
expect_error(MultilinearICA(X, modes=1:4))
expect_error(MultilinearICA(X, modes=c(1, 3)))

### Test E-4: algorithm
expect_error(MultilinearICA(X, algorithm="FastICAA"))

## All Combination of J and modes
expect_error(expect_error(MultilinearICA(X, Js=c(2,3,4), modes=c(1,2,3))))
expect_error(expect_error(MultilinearICA(X, Js=c(2,4), modes=c(1,3))))
expect_error(expect_error(MultilinearICA(X, Js=c(3,4), modes=c(2,3))))
expect_error(expect_error(MultilinearICA(X, Js=2, modes=1)))
expect_error(expect_error(MultilinearICA(X, Js=3, modes=2)))
expect_error(expect_error(MultilinearICA(X, Js=4, modes=3)))
