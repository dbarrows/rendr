test_that('RRE works', {
    sol <- rsys_examples('mm') %>%
        rre(length.out = 100) %>%
        .$sol
    expect_equal(nrow(sol), 100)
    expect_equal(ncol(sol), 5)
})