test_that('ISSA works', {
    sol <- rdsys_examples('schnakenberg') %>%
        issa(length.out = 100, verbose = FALSE)
    expect_equal(length(sol$t), 100)
    expect_equal(length(sol$u), 100)
})