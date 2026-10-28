# Wilcoxon rank-sum p-value with the pre-R 4.6 rule (exact only for n < 50 and no ties, else normal
# approximation), as used in the original study; fixed so results do not depend on the R version.
wilcox_p <- function(x, y, alternative = "two.sided") {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  exact <- length(x) < 50 && length(y) < 50 && !anyDuplicated(c(x, y))
  wilcox.test(x, y, alternative = alternative, exact = exact)$p.value
}
