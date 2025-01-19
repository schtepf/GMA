test_that("signed logarithm (base 2)", {
  expect_equal(signed.log(c(-1, 0, 1, 3, 7), base=2), 
               c(-1, 0, 1, 2, 3))
})
