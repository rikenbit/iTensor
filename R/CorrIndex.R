#' CorrIndex
#'
#' Calculate the CorrIndex of the cross-correlation matrix of S_true and estimated S.
#' The closer the value is to 0, the closer estimated S is to S_true.
#' @param cross_correlation_matrix Cross-correlation matrix
#' @return CorrIndex, which means the closeness between S and S_true , is returned.
#' @examples
#' S_true <- matrix(runif(5*5), nrow=5, ncol=5)
#' S <- matrix(runif(5*5), nrow=5, ncol=5)
#' CorrIndex(cor(S_true, S))
#' @export
CorrIndex <- function(cross_correlation_matrix){
    # 1. Calculate the absolute value of the elements.
    abs_cross_correlation_matrix <- abs(cross_correlation_matrix)
    # 2. Calculate the maximum value for each row.
    max_row <- apply(abs_cross_correlation_matrix, 1, max)
    # 3. Calculate the maximum value for each column.
    max_col <- apply(abs_cross_correlation_matrix, 2, max)
    # 4. Divide by the maximum value for each row.
    row_divided <- abs_cross_correlation_matrix / max_row
    # 5. Divide by the maximum value for each column.
    col_divided <- abs_cross_correlation_matrix / max_col
    # 6. Calculate the sum for each row.
    row_sum <- apply(row_divided, 1, sum)
    # 7. Calculate the sum for each column.
    col_sum <- apply(col_divided, 2, sum)
    # 8. Subtract 1 from the sum for each row.
    row_subtracted <- row_sum - 1
    # 9. Subtract 1 from the sum for each column.
    col_subtracted <- col_sum - 1
    # 10. Calculate the sum of the difference between the sum for each row and 1 and the sum for each column and 1.
    sum(row_subtracted) + sum(col_subtracted)
}
