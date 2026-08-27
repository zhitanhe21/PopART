# Printed output ------------------------------------------------------------

make_test_fit <- function() {
  format_results <- getFromNamespace(".format_popart_results", "popart")
  contribution <- cbind(c(0.2, 0.4, 0.3), c(0.1, 0.5, 0.6))
  fit <- format_results(list(
    gformula_naive = contribution,
    gformula_proposed = contribution,
    ipw_naive = contribution,
    ipw_proposed = contribution,
    aipw_naive = contribution,
    aipw_proposed = contribution
  ))
  class(fit) <- "popart_fit"
  fit
}


test_that("print shows one row per estimator and version", {
  text <- capture.output(print(make_test_fit()))

  expect_true(any(grepl("PopART causal estimates", text, fixed = TRUE)))
  expect_true(any(grepl("eta(0)", text, fixed = TRUE)))
  expect_true(any(grepl("AIPW", text, fixed = TRUE)))
})


test_that("summary returns estimates, standard errors, and intervals", {
  result <- summary(make_test_fit())

  expect_equal(nrow(result), 24L)
  expect_named(result, c(
    "estimator", "version", "parameter", "estimate", "std_error", "95% CI"
  ))
})
