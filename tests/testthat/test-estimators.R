# Estimator output ----------------------------------------------------------

test_that("all six estimators are returned", {
  format_results <- getFromNamespace(".format_popart_results", "popart")

  contribution <- cbind(c(0.2, 0.4, 0.3), c(0.1, 0.5, 0.6))
  results <- format_results(list(
    gformula_naive = contribution,
    gformula_proposed = contribution,
    ipw_naive = contribution,
    ipw_proposed = contribution,
    aipw_naive = contribution,
    aipw_proposed = contribution
  ))

  expect_equal(nrow(results$estimates), 24L)
  expect_equal(nrow(results$variance), 24L)
  expect_equal(length(results$covariance), 6L)
  expect_setequal(
    unique(results$estimates$estimator),
    c("G-Formula", "IPW", "AIPW")
  )
  expect_setequal(unique(results$estimates$version), c("Naive", "Proposed"))
})


test_that("risk difference and risk ratio use the two arm means", {
  format_results <- getFromNamespace(".format_popart_results", "popart")

  contribution <- cbind(c(0.2, 0.3, 0.4), c(0.3, 0.5, 0.7))
  estimates <- format_results(
    list(gformula_naive = contribution)
  )$estimates

  expect_equal(
    estimates$estimate[estimates$parameter == "risk_difference"],
    0.2
  )
  expect_equal(
    estimates$estimate[estimates$parameter == "risk_ratio"],
    5 / 3
  )
})


test_that("variance is grouped by cluster", {
  format_results <- getFromNamespace(".format_popart_results", "popart")

  contribution <- cbind(c(0.1, 0.2, 0.4, 0.5), c(0.2, 0.3, 0.5, 0.8))
  cluster <- c(1, 1, 2, 2)
  result <- format_results(
    list(gformula_naive = contribution), cluster
  )

  centered <- sweep(contribution, 2L, colMeans(contribution))
  cluster_sum <- rowsum(centered, cluster)
  expected <- crossprod(cluster_sum) * 2 / (1 * nrow(contribution)^2)

  expect_equal(unname(result$covariance$gformula_naive), unname(expected))
})
