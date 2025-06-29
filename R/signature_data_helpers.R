#' Universal Data Access Helpers for Signature Data
#'
#' These functions provide consistent access to signature data columns
#' that may have different names depending on whether they're from
#' pan-cluster aggregations or individual signature results.

#' Get signature strength value from data frame
#'
#' @param df Data frame with signature data
#' @return Vector of signature strength values
#' @export
get_signature_strength <- function(df) {
  if ("signature_strength" %in% colnames(df)) {
    return(df$signature_strength)
  } else if ("mean_signature_strength" %in% colnames(df)) {
    return(df$mean_signature_strength)
  } else if ("max_signature_strength" %in% colnames(df)) {
    return(df$max_signature_strength)
  } else {
    warning("No signature strength column found")
    return(rep(NA, nrow(df)))
  }
}

#' Get cluster information from data frame
#'
#' @param df Data frame with signature data
#' @return Vector of cluster information
#' @export
get_cluster_info <- function(df) {
  if ("cluster" %in% colnames(df)) {
    return(df$cluster)
  } else if ("cluster_count" %in% colnames(df)) {
    return(df$cluster_count)
  } else {
    warning("No cluster information column found")
    return(rep(NA, nrow(df)))
  }
}

#' Get signature metric safely from data frame
#'
#' @param data Data frame with signature data
#' @param metric Character, metric to extract ("strength", "cluster_info", "gene_pair")
#' @return Vector of values for the requested metric
#' @export
get_signature_metric <- function(data, metric) {
  switch(metric,
    "strength" = get_signature_strength(data),
    "cluster_info" = get_cluster_info(data),
    "gene_pair" = {
      if ("gene_pair" %in% colnames(data)) {
        data$gene_pair
      } else {
        warning("No gene_pair column found")
        rep(NA, nrow(data))
      }
    },
    stop("Unknown metric: ", metric)
  )
}

#' Validate signature data structure
#'
#' @param data Data frame to validate
#' @param required_cols Character vector of required column names
#' @return Logical, TRUE if valid
#' @export
validate_signature_data <- function(data, required_cols = c("gene_pair")) {
  if (is.null(data) || nrow(data) == 0) {
    stop("Empty or NULL signature data")
  }
  
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  return(TRUE)
}

#' Safe maximum calculation for signature strength
#'
#' @param data Data frame with signature data
#' @return Numeric, maximum signature strength or 0 if not available
#' @export
safe_max_signature_strength <- function(data) {
  strength_values <- get_signature_strength(data)
  if (all(is.na(strength_values))) {
    return(0)
  }
  return(max(strength_values, na.rm = TRUE))
}

#' Safe column access that handles both individual and aggregated signature data
#'
#' @param signature_obj Signature object (could be individual signature or aggregated)
#' @param column_name Character, name of column to access
#' @return Value from column or appropriate fallback
#' @export
safe_signature_access <- function(signature_obj, column_name) {
  switch(column_name,
    "signature_strength" = {
      if (!is.null(signature_obj$signature_strength)) {
        signature_obj$signature_strength
      } else if (!is.null(signature_obj$mean_signature_strength)) {
        signature_obj$mean_signature_strength
      } else if (!is.null(signature_obj$max_signature_strength)) {
        signature_obj$max_signature_strength
      } else {
        NA
      }
    },
    "cluster" = {
      if (!is.null(signature_obj$cluster)) {
        signature_obj$cluster
      } else if (!is.null(signature_obj$cluster_count)) {
        signature_obj$cluster_count
      } else {
        NA
      }
    },
    # Default: try direct access
    signature_obj[[column_name]]
  )
}