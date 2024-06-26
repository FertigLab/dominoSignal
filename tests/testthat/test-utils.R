test_that("bool conversion function works",{
  df <- data.frame(list(c1 = c("True", "False"),
                        c2 = c("False", "True"),
                        c3 = c(1, 2)),
                        c4 = c("a", "b"))

  c_df <- conv_py_bools(df)

  expect_equal(class(c_df$c1), "logical")
  expect_equal(class(c_df$c2), "logical")
  expect_equal(class(c_df$c3), "numeric")
  expect_equal(class(c_df$c4), "character")

})

test_that("read if char tries to read a file", {
  expect_error(read_if_char("./file_that_not_exists.csv",
                            "cannot open the connection"))
  expect_error(read_if_char(c('a', 'b')), "Length of obj must be one of: 1")
})

test_that("mandatory field absence yields error, presence does not", {
  expect_error(check_arg(arg = data.frame(a = c(1, 2), b = c(3, 4)),
                         allow_class = "data.frame",
                         need_vars = c("c")))
  expect_silent(check_arg(arg = data.frame(a = c(1, 2), b = c(3, 4)),
                          allow_class = "data.frame",
                          need_vars = c("a", "b")))
})

test_that("range checker works", {
  expect_error(check_arg(1, allow_class = "numeric", allow_range = c(2, 5)),
               "All values in 1 must be between 2 and 5")
  expect_silent(check_arg(3, allow_class = "numeric", allow_range = c(2, 5)))
})

test_that("check_arg works for class", {
  expect_error(check_arg(1, allow_class = "data.frame"),
               "Class of 1 must be one of: data.frame")
  expect_silent(check_arg(data.frame(a = c(1, 2), b = c(3, 4)),
                          allow_class = "data.frame"))
  expect_silent(check_arg(data.frame(a = c(1, 2), b = c(3, 4)),
                          allow_class = c("data.frame", "matrix")))
})
