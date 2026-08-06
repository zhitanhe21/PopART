###############################################################################
# Lightweight reporting helpers
###############################################################################


#' Convert a data frame to a Markdown table.
#'
#' @param df Data frame to print.
#' @param digits Integer number of decimal places for numeric columns.
#'
#' @return Character scalar containing a Markdown table.
md_table <- function(df, digits = 3) {
  if (is.null(df) || nrow(df) == 0L) return("_No rows available._")
  out <- as.data.frame(df, stringsAsFactors = FALSE)
  out[] <- lapply(out, function(x) {
    if (is.numeric(x)) formatC(x, format = "f", digits = digits) else as.character(x)
  })
  header <- paste(names(out), collapse = " | ")
  sep <- paste(rep("---", ncol(out)), collapse = " | ")
  rows <- apply(out, 1, paste, collapse = " | ")
  paste(c(
    paste0("| ", header, " |"),
    paste0("| ", sep, " |"),
    paste0("| ", rows, " |")
  ), collapse = "\n")
}
