test_that("NSM works", {
    sol <- rdsys_examples('schnakenberg') %>%
        nsm(length.out = 100, verbose = FALSE)
    expect_equal(length(sol$t), 100)
    expect_equal(length(sol$u), 100)
})