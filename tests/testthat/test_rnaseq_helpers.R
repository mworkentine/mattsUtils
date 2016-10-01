library(mattsUtils)
library(DESeq2)

context("rnaseq helpers")


load("rnaseq_test_data.rda")


test_that("run_spia runs with reactome", {
  spia_res_reactome = run_spia(res, dds, "dmelanogaster", "reactome", key = key)
  expect_equal_to_reference(spia_res_reactome,
                            "spia_res_dmelanogaster_reactome.rds")

})

test_that("run_spia runs with kegg", {
  spia_res_kegg = run_spia(res, dds, "dmelanogaster", "kegg", key = key)
  expect_equal_to_reference(spia_res_kegg,
                            "spia_res_dmelanogaster_kegg.rds")

})
