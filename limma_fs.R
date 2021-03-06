#' Author: David Pinyeiro
#' 
#' Select CpGs by differential methylation using limma
#'
#' This funtion uses limma to select top differentially methylated probes
#' (features) between groups (sample labels factor). The number of class labels
#' can be 2 or more, but no paired designs are implemented. It is equivalent to
#' anova_fs, but much faster. Only applicable to methylation array data.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the number of features (top differentially methylated CpGs
#'   to select).
#'
#' @return [character] A character vector with the \code{n} top differentially
#'   methylated (sorted by p-value) feature names, for all contrasts tested.
#'   For more than two groups in \code{sample_labels}, each contrast is designed
#'   as one Vs the others.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- limma_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export

limma_fs <- function(input_data, sample_labels, n) {
  # First, make sure sample_labels has valid names
  sample_labels <- factor(make.names(sample_labels))
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else if (sum(input_data == 0) | sum(input_data == 1)) {
    # Check whether absolute 0 or 1 values exists. So it is usually a sign of
    # sequencing data.
    warning("0 and/or 1 are found in your data. Is this data generated by
            bisulfite sequencing?. Limma filtering is not applied.")
    return(colnames(input_data))
  } else {
    # Calculate M-values using M=log2(Beta/(1-Beta)). All statistics will be
    # performed on M-values. Transpose to have rows as features and columns as
    # samples.
    
    sample_names <- rownames(input_data)
    input_data <- t(input_data)
    m_vals <- log2(input_data / (1-input_data))
    # Create a targets dataframe.
    pheno_data <- data.frame(sample_names, sample_labels)
    rownames(pheno_data) <- sample_names
    targets <- stats::model.frame(sample_names ~ sample_labels, pheno_data)
    # Design matrix (only unpaired test supported).
    design <- stats::model.matrix(~0+sample_labels, data=targets)
    colnames(design) <- levels(sample_labels)
    # Contrast matrix (one vs the others).
    if (nlevels(sample_labels) == 2) {
      contr <- paste0(levels(sample_labels)[1], "-", levels(sample_labels)[2])
      contMatrix <- limma::makeContrasts(contrasts = contr, levels = design)
    } else {
      # More than 2 groups. Using One Vs the Others contrasts.
      i <- 1
      contr <- character()
      while (i <= nlevels(sample_labels)) {
        one <- levels(sample_labels)[i]
        the_others <- levels(sample_labels)[-i]
        # We will make contrasts such as "A-(B+C)/2".
        the_others_mean <- paste0("(", paste(the_others, collapse = "+"),
                                  ")/", length(the_others))
        contr <- c(contr, paste0(one, "-", the_others_mean))
        i <- i + 1
      }
      contMatrix <- limma::makeContrasts(contrasts = contr, levels = design)
    }
    # fit the linear model
    fit <- limma::lmFit(m_vals, design)
    # fit the contrasts
    fit2 <- limma::contrasts.fit(fit, contMatrix)
    fit2 <- limma::eBayes(fit2)
    # Toptable of all contrasts, selected by p-value.
    list_of_results <- lapply(1:(nlevels(sample_labels) - 1), function(x) {
      limma::topTable(fit2, coef = x, number=Inf, sort.by = "P")
    })
    # Selecting the "best" CpGs of each contrast. If
    # n/nlevels(sample_labels) has no integer result, round
    # approximation is taken.
    each_contrast_n <- round(n/(nlevels(sample_labels) - 1), 0)
    feature_selection <- character()
    for (r in list_of_results) {
      feature_selection <- c(feature_selection, rownames(r)[1:each_contrast_n])
    }
    return(feature_selection)
  }
}
