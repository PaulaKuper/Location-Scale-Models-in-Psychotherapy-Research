# Function to extract estimate, se and df of models
getfit <- function(model) {
  if (model$test == "z") {
    cbind(estimate = c(coef(model)$beta,
                       coef(model)$alpha), 
          se = c(sqrt(diag(vcov(model)$beta)),
                 sqrt(diag(vcov(model)$alpha))), 
          df = Inf)
  }
  else {
    cbind(estimate = c(coef(model)$beta,
                       coef(model)$alpha), 
          se = c(sqrt(diag(vcov(model)$beta)),
                 sqrt(diag(vcov(model)$alpha))),
          df = model$k -model$p)
  }
}

# Function to check for missing values in a data frame
check_missing <- function(df) {
  df <- as.data.frame(df)  # Convert matrix to data frame
  missing_values <- df %>% summarise_all(~ any(is.na(.)))
  return(missing_values)
}

# Function to check if any value in a data frame is TRUE
check_any_true <- function(df) {
  any(unlist(df))
}

# Function to fit all possible models in multimodel selection
fit_all_models <- function(formula) {
  rma(yi=yi, sei=sei, mods=~1,scale=formula, data=data,
      method= "REML", test= "knha", control= list(optimizer="BFGS", rel.tol=1e-8))
}


# Function to compute Chi-square test and Cramér's V for each pair of variables
compute_associations <- function(data) {
  # Initialize an empty matrix to store results
  association_matrix <- matrix(NA, nrow = ncol(data), ncol = ncol(data), dimnames = list(names(data), names(data)))
  
  # Loop through pairs of variables
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      if(i == j) {
        association_matrix[i, j] <- 1
      } else if(i < j) {
        # Compute Chi-square test
        chi_square <- chisq.test(table(data[, i], data[, j]))
        # Compute Cramér's V
        cramer_v <- assocstats(table(data[, i], data[, j]))$cramer 
        # Store results in the matrix
        association_matrix[i, j] <- paste0(round(cramer_v, 3)) %>% as.numeric()
      } else {
        association_matrix[i, j] <- association_matrix[j, i]
      }
    }
  }
  # Return the association matrix
  return(association_matrix)
}

