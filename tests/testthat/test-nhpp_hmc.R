capture_output(hmc_fit <- nhpp_hmc(y~1,data=data.frame(y=sample(1:9,100,replace = T)),warm_up = 5,iter_max = 10))
test_that("nhpp errors correctly", {
    expect_error(nhpp_hmc(y~1,data=data.frame(y=sample(1:9,100,replace = T)),warm_up = 10,iter_max = 5))
})

test_that("nhpp prints correctly",{
    expect_output(print(nhpp_hmc(y~1,data=data.frame(y=sample(1:9,100,replace=T)),warm_up = 5, iter_max = 10)),
                  regexp = "Regression")
})

test_that("nhpp returns correct type",{
    expect_is(hmc_fit$beta,"mcmc")
})

