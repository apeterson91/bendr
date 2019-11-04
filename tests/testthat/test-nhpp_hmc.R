test_that("nhpp errors correctly", {
    expect_error(nhpp_hmc(y~1,data=data.frame(y=sample(1:9,100)),warm_up = 10,iter_max = 5))
})

