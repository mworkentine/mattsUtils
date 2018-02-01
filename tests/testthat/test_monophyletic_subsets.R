library(mattsUtils)
context("Test monophyletic subsets")

set.seed(9981)
tt = rtree(10)

test_sets = list(
  all = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10"),
  mono1 = c("t5", "t9", "t10"),
  mono2 = c("t8", "t2"),
  all_but_one = c("t10", "t9", "t5", "t6", "t4", "t7", "t8"),
  all_but_one2 = c("t1", "t7", "t4"),
  none = c("t10", "t5", "t4", "t1", "t2"),
  none2 = c("t5", "t4", "t1", "t2"),
  all_but_two = c("t7", "t6", "t9", "t10"),
  two_not_one = c("t10", "t9", "t2", "t8", "t6"),
  two_not_two = c("t10", "t9", "t5", "t1", "t3", "t2", "t7")
)

expect_sets = list(
  all = list(c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10")),
  mono1 = list(c("t9", "t10")),
  mono2 = list(c("t8", "t2")),
  all_but_one = list(c("t10", "t9", "t5", "t6", "t4", "t7")),
  all_but_one2 = c("t1", "t7", "t4"),
  none = c("t10", "t5", "t4", "t1", "t2"),
  none2 = c("t5", "t4", "t1", "t2"),
  all_but_two = c("t7", "t6", "t9", "t10"),
  two_not_one = c("t10", "t9", "t2", "t8", "t6"),
  two_not_two = c("t10", "t9", "t5", "t1", "t3", "t2", "t7")
)


test_that("get_monophyletic_subsets finds groups ", {
  hits = get_monophyletic_subsets(tt, test_sets$all)

  expect_equal(sort(hits[[1]]), sort(test_sets$all))
  expect_equal(len)
})
