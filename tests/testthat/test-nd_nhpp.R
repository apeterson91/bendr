ITER <- 10
WARM <- 5
test_that("nd_nhpp function errors if given negative distances", {
    expect_error(nd_nhpp(c(-1,1,1),n_j = matrix(1,2,2),
                         iter_max = ITER, warm_up = WARM,),
                         regexp = "must be positive")
    expect_error(nd_nhpp(c(-1,1,1),n_j = matrix(1,2,2),
                         iter_max = WARM, warm_up = ITER),
                 regexp = "must be <")

})
