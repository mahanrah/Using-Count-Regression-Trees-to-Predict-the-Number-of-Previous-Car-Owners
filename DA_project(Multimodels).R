# Install the packages if not already installed
if (!require(MASS)) install.packages("MASS")
if (!require(pscl)) install.packages("pscl")
if (!require(glmmTMB)) install.packages("glmmTMB")
install.packages
# Package
library(MASS)
library(pscl)
library(glmmTMB)
# Load in the dataset
data <- read.csv("imputed_cars.csv")
data <- subset(data, select = -c(Color, model, Seller.Type, make, Location)) # Remove unnecessary columns
data <- na.omit(data)
data <- data[!apply(data == "", 1, any), ]
data_small <- data[1:1000,]

#Initialise list of final nodes - the bottom leaves of the tree
nodes_list <- list()

#Set up models 
#======== Set Up Models Function ========
fit_model <- function(formula, data, family = gaussian(), zero_inflated = FALSE, hurdle = FALSE) {
  if (zero_inflated || hurdle) {
    if (any(data$owners == 0)) {
      if (hurdle) {
        model <- glmmTMB(as.formula(formula), data = data, family = nbinom2, ziformula = ~1, dispformula = ~1)
      } else {
        model <- zeroinfl(as.formula(formula), data = data, dist = "negbin", link = "logit")
      }
    } else {
      stop("Zero-inflated or hurdle models are not applicable: no zero counts present.")
    }
  } else {
    if (family == "nbinom") {
      model <- glm.nb(as.formula(formula), data = data)
    } else {
      model <- glm(as.formula(formula), data = data, family = poisson())
    }
  }
  return(model)
}

#======== Choose Best Model ========
select_model <- function(data, formula) {
  
  
  model_specs <- list(
    list(name = "Poisson", family = "poisson"),
    list(name = "Negative Binomial", family = "nbinom"),
    list(name = "Zero-Inflated Poisson", zero_inflated = TRUE),
    list(name = "Zero-Inflated Negative Binomial", zero_inflated = TRUE, family = "nbinom"),
    list(name = "Hurdle Poisson", zero_inflated = TRUE, hurdle = TRUE),
    list(name = "Hurdle Negative Binomial", zero_inflated = TRUE, hurdle = TRUE, family = "nbinom")
  )
  
  best_model <- NULL
  lowest_aic <- Inf
  best_name <- NULL
  
for (spec in model_specs) {
  tryCatch({
    # Remove 'name' from the specs before passing to fit_model
    params <- spec[names(spec) != "name"]
    model <- do.call(fit_model, c(list(formula = formula, data = data), params))
    model_aic <- AIC(model)
    cat(sprintf("Model: %s, AIC: %f\n", spec$name, model_aic))
    if (model_aic < lowest_aic) {
      lowest_aic <- model_aic
      best_model <- model
      best_name <- spec$name
    }
  }, error = function(e) {
    cat(sprintf("Failed to fit %s model: %s\n", spec$name, e$message))
  })
}

  
  cat(sprintf("%s model selected at this node\n", best_name))
  return(best_model)
}

count_regression_tree <- function(data, original_data, nodes_list = list(), iteration = 1){
  print("Level")
  print(iteration)
  
  ####################################################################################
  #Choose a model - Poisson/Negative Binomial/Hurdle/Zero-Inflated
  #Fit a regression model
  control <- glm.control(maxit = 100000)  # Increase the max number of iterations
  #Identify columns that are factors or characters (potentially categorical)
  categorical_vars <- names(data)[sapply(data, function(x) is.factor(x) || is.character(x))]
  #Identify which categorical variables have more than one unique level
  valid_categorical_vars <- sapply(categorical_vars, function(var) length(unique(data[[var]])) > 1)
  #Prepare a model formula excluding single-level categorical variables
  formula <- paste("owners ~", paste(names(valid_categorical_vars[valid_categorical_vars]), collapse = " + "), 
                   "+ year + price + mileage + engine.size + seats + 
                   fuel.economy + max.speed + engine.power")
  #Fit the best model for the given node
  model <- select_model(data, formula)
  
  # Find the partial score residuals for each covariate
  model_residuals <- residuals(model, type = "response") # response = y_i - u(hat)_i
  #nb_model_residuals <- residuals(nb_model, type = "response") # response = y_i - u(hat)_i
  
  # Add their sign (+ or -) to the end of the dataframe
  data$ResidualSign <- ifelse(model_residuals >= 0, "Positive", "Negative")
  
  # ************** NEED TO IMPUTE DATA HERE ************************************
  
  ####################################################################################
  # Use the residual contingency tables to perform the 
  # Pearson Chi-Squared Test of Independence on the Residuals:
  # - Count the number of positive and negative residuals under the different quartiles/categories
  # - Perform test on these tables
  # - Calculate W1 and W2
  
  # Initialise...
  table_list <- list() #  The list for storing the tables
  w_1_vector <- c() #  The vector for storing W1 values
  # Remove count column (owners) and residual column
  covariates <- setdiff(colnames(data), c("owners", "ResidualSign"))
  original_colnames <- setdiff(colnames(data), "ResidualSign")
  
  for (i in covariates) # Loop through each covariate
  {
    # Order in terms of that covariate
    data <- data[order(data[,i]), ]
    # Initialise count_table
    count_table <- matrix()
    # For numerical columns - split into quartiles (the covariate column)
    if (is.numeric(data[,i]))
    {
      
      # # Calculate quartiles and add as a new column
      # data$quartiles <- cut(data[,i], breaks = c(-Inf, unique(quantile(data[,i], probs = c(0.25, 0.5, 0.75))), Inf),
      #                      include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
      
      # Calculate quantiles and ensure they are unique
      quantiles <- unique(quantile(data[, i], probs = c(0.25, 0.5, 0.75)))
      # Add -Inf and Inf to cover ranges outside the quantiles
      breaks <- c(-Inf, quantiles, Inf)
      # Generate appropriate labels based on the number of unique breaks
      labels <- switch(length(quantiles),
                       c("Q1", "Q2"),
                       c("Q1", "Q2", "Q3"),
                       c("Q1", "Q2", "Q3", "Q4"))
      # Use cut with these adjusted breaks
      data$quartiles <- cut(data[, i], breaks = breaks, include.lowest = TRUE, labels = labels)
      
      # Create empty contingency table
      count_table <- matrix(0, nrow = 2, ncol = length(labels), dimnames = list(c("Positive", "Negative"), labels))
      
      # Increment for each row in the dataframe
      for (j in 1:nrow(data)) {
        # Check the quartile and the sign of the row
        quartile <- data$quartiles[j]
        residual_sign <- data$ResidualSign[j]
        # Increment the table
        count_table[residual_sign, quartile] <- count_table[residual_sign, quartile] + 1
      }
      data <- subset(data, select = -c(quartiles)) # Remove quartiles column
      
    }else # For categorical columns - split into the different categories
    {
      # Create a vector of the different categorical options for the column
      categories <- unique(data[,i])
      # Create empty contingency table with these categories as columns
      count_table <- matrix(0, nrow = 2, ncol = length(categories), dimnames = list(c("Positive", "Negative"), categories))
      # Increment for each row in the dataframe
      for (j in 1:nrow(data)) {
        # Check the category and the sign of the row
        category <- data[j,i]
        residual_sign <- data$ResidualSign[j]
        # Increment the table
        count_table[residual_sign, category] <- count_table[residual_sign, category] + 1
      }
    }
    # Add table to list, naming the list element after the covariate
    table_list <- append(table_list, list(count_table))
    
    # Pearson Chi-Squared Test of Independence
    chi <- chisq.test(count_table)
    x_squared <- chi$statistic
    # Wilson-Hilferty Approximation
    v <- nrow(count_table) # v = number of rows in contingency table
    w_1 <- max(0,((7/9)+sqrt(v)*(((x_squared/v)^(1/3))-1+(2/(9*v))))^3)
    # Add to vector
    w_1_vector <- c(w_1_vector, w_1)
  }
  
  
  
  ##################################################################################################
  # Select split variable
  alpha_1 <- 0.05/length(covariates) # α1 value = 0.003846154
  alpha_2 <- 0.1/(length(covariates)*(length(covariates)-1)) # α2 value = 0.0006410256
  
  # Find the 	100(1-α)th percentile of the chi-squared distribution with one degree of freedom
  x_squared_alpha1 <- qchisq(1-alpha_1, df = 1)
  x_squared_alpha2 <- qchisq(1-alpha_2, df = 1)
  
  # Initialise...
  split_var <- NA # Split variable
  
  # Split Variable Choice
  if (max(w_1_vector, na.rm = TRUE)>x_squared_alpha1)
  {
    # Split variable is that which is associated with the largest value of x1
    split_var <- covariates[which.max(w_1_vector)]
    
  }else # If not, W_2 must be calculated - interaction effect detection
  {
    # Initialise...
    table_list_interaction <- list() #  The list for storing the tables
    w_2_vector <- c() #  The vector for storing W2 values
    # The matrix for storing the covariates used for each W2 value
    w_2_vector_covariates <- matrix(data = NA, nrow = 2,
                                    ncol = length(covariates)*(length(covariates))-1) 
    # The index for number of iterations
    index = 0
    for (i in covariates) # Loop through each covariate
    {
      # First, identify which groups each observation is in for the contingency table
      # Order in terms of that covariate
      data <- data[order(data[,i]), ]
      # Initialise count_table
      count_table <- matrix()
      # Initialise column names for count_table
      count_table_colnames <- c()
      # For numerical columns - split into halves (the covariate column)
      if (is.numeric(data[,i])) 
      {
        # Add to column names for contingency table
        count_table_colnames <- c("H1_i", "H2_i")
        # Calculate halves and add as a new column
        data$halves_i <- cut(data[,i], breaks = c(-Inf, unique(quantile(data[,i], probs = c(0.5))), Inf),
                             include.lowest = TRUE, labels = c("H1_i", "H2_i"))
        
        # # Calculate quantiles and ensure they are unique
        # quantiles <- unique(quantile(data[, i], probs = c(0.5)))
        # # Add -Inf and Inf to cover ranges outside the quantiles
        # breaks <- c(-Inf, quantiles, Inf)
        # # Generate appropriate labels based on the number of unique breaks
        # labels <- c("H1", "H2")
        # # Use cut with these adjusted breaks
        # data$quartiles <- cut(data[, i], breaks = breaks, include.lowest = TRUE, labels = labels)
        # 
        
      }else # For categorical columns - split into the different categories
      {
        # Create a vector of the different categorical options for the column
        count_table_colnames <- unique(data[,i])
      }
      
      # Remove i from covariates
      covariates_without_i <- setdiff(covariates, i)
      for (j in covariates_without_i) # Loop through again for second covariate
      {
        index = index+1 # Increment index
        # For numerical columns - split into halves (the covariate column)
        if (is.numeric(data[,j])) 
        {
          # Add to column names for contingency table - combine with other colnames
          grid <- expand.grid(count_table_colnames, c("H1_j", "H2_j"))
          combined_count_table_colnames <- paste(grid$Var1, grid$Var2, sep = "_")
          # Calculate halves and add as a new column
          data$halves_j <- cut(data[,j], breaks = c(-Inf, quantile(data[,j], probs = c(0.5)), Inf),
                               include.lowest = TRUE, labels = c("H1_j", "H2_j"))
        }else # For categorical columns - split into the different categories
        {
          # Add to column names for contingency table - combine with other colnames
          grid <- expand.grid(count_table_colnames, unique(data[,j]))
          combined_count_table_colnames <- paste(grid$Var1, grid$Var2, sep = "_")
        }
        # Create matrix with new column names
        count_table <- count_table <- matrix(0, nrow = 2, ncol = length(combined_count_table_colnames),
                                             dimnames = list(c("Positive", "Negative"), combined_count_table_colnames))
        # Iterate through the rows in the dataframe and increment the table
        for (h in 1:nrow(data)) 
        {
          category1 <- 0 # Initialise 'category1'
          category2 <- 0 # Initialise 'category2'
          # Check the column and the sign of the row
          # Column i
          if (is.numeric(data[,i])) 
          {
            # Select the half value for row h 
            category1 <- data$halves_i[h]
          }else
          {
            # Select the category value for row h
            category1 <- data[h,i]
          }
          # Column j
          if (is.numeric(data[,j])) 
          {
            # Select the half value for row h 
            category2 <- data$halves_j[h]
          }else
          {
            # Select the category value for row h
            category2 <- data[h,j]
          }
          
          residual_sign <- data$ResidualSign[h]
          # Increment the table
          count_table[residual_sign, paste(category1, category2, sep = "_")] <- 
            count_table[residual_sign, paste(category1, category2, sep = "_")] + 1
        }
        # Add table to list, naming the list element after the covariate
        table_name <- paste(i, j, sep = "_")
        table_list_interaction <- c(table_list_interaction, table_name = count_table)
        
        # Pearson Chi-Squared Test of Independence
        chi <- chisq.test(count_table)
        x_squared <- chi$statistic
        # Wilson-Hilferty Approximation
        v <- nrow(count_table) # v = number of rows in contingency table
        w_2 <- max(0,((7/9)+sqrt(v)*(((x_squared/v)^(1/3))-1+(2/(9*v))))^3)
        # Add to vector
        w_2_vector <- c(w_2_vector, w_2)
        # Add to matrix
        w_2_vector_covariates[1,index] <- i
        w_2_vector_covariates[2,index] <- j
      }
    }
    
    # Split Variable Choice
    if (max(w_2_vector, na.rm = TRUE)>x_squared_alpha1)
    {
      # Split variable is one of the two associated with the largest value of w2
      split_var_i <- w_2_vector_covariates[1,which.max(w_2_vector)]
      split_var_j <- w_2_vector_covariates[2,which.max(w_2_vector)]
      
      # Choose the split variable as the one of these two with the largest w1
      # Check for NA values
      if (is.na(w_1_vector[which(covariates==i)]))
      {
        if (is.na(w_1_vector[which(covariates==j)]))
        {
          split_var <- NA # Both have NA W1s
        }else
        {
          split_var <- split_var_j # i has an NA W1
        }
      }else
      {
        if(is.na(w_1_vector[which(covariates==j)]))
        {
          split_var <- split_var_i # j has an NA W1
        }else 
        {
          # Neither W1s are NA
          if (w_1_vector[which(covariates==i)] > w_1_vector[which(covariates==j)])
          {
            split_var <- split_var_i
          }else
          {
            split_var <- split_var_j
          }
        }
      }
    }
    
  }
  
  ##################################################################################################
  # Select split set
  print("split var")
  print(split_var)
  # Initialise...
  split_set_1_list <- list() # List of possible split sets (first)
  split_set_2_list <- list() # List of possible split sets (second)
  node_impurity_reduction <- c() # Vector of sizes of impurity reduction for each split 
  split_variable_is_zero <- FALSE # Boolean stating that the recursion stops here
  
  if (is.na(split_var)) # No split variable selected => no node split
  {
    # No node split => change boolean to true
    split_variable_is_zero <- TRUE
   
  }else
  {
    if (is.numeric(data[,split_var]))
    {
      # For numerical variables, try all possible split points 
      # (1 in group A, rest in B; 2 in A, rest in B; 3 in A,rest in B etc)
      split_points <- sort(unique(data[,split_var]))
      split_points <- split_points[-1]
      # Loop through each possible split of the data by numeric value
      for (i in split_points)
      {
        # Split the data frame into two (above and below i)
        split_set_1 <- data[data[, split_var] < i, ]
        split_set_2 <-data[data[, split_var] >= i, ]
        # Add these to the split set lists
        split_set_1_list <- append(split_set_1_list, list(split_set_1))
        split_set_2_list <- append(split_set_2_list, list(split_set_2))
      }
    }else
    {
      # For categorical variables, try all possible splits of categories
      # (Category 1 in A, rest in B; category 2 in A, rest in B; category 1 and 2 in A, rest in B etc)
      # Function to generate all non-empty subsets of a set
      subsets <- list()
      # Get unique categories from split variable column 
      unique_categories <- unique(data[, split_var])
      # Generate all possible non-empty subsets
      for (index in 1:(floor(length(unique_categories)/2))) # Only need to use combinations of half the size or less of the categories => other half is constituted by the other split set of values not selected as a category combination
      {
        if (index <= length(unique_categories))
        {
          # Create unique combinations
          combinations <- combn(unique_categories, index, simplify = FALSE)
          # Place in list form and append to subsets list
          subsets <- append(subsets, combinations)
        }
      }
      # Loop through each possible split of the data by category
      for (i in subsets)
      {
        # Split the data frame into two (one category subset)
        split_set_1 <- data[data[, split_var] %in% i, ]
        split_set_2 <-data[!data[, split_var] %in% i, ]
        # Add these to the lists
        split_set_1_list <- append(split_set_1_list, list(split_set_1))
        split_set_2_list <- append(split_set_2_list, list(split_set_2))
      }
    }
  
    # Test each split of the dataset to see which produces the maximum node impurity reduction
    for (i in 1:length(split_set_1_list)) 
    {
      # Re-initialise the selected split sets
      split_set_1 <- split_set_1_list[[i]]
      split_set_2 <- split_set_2_list[[i]]
      # Fit a model for each 
      #Initialise...
      split_set_1_model <- NA
      split_set_2_model <- NA
      
      if (nrow(split_set_1)>0)
      {
        # Split Set 1
        # Identify columns that are factors or characters (potentially categorical)
        categorical_vars <- names(split_set_1)[sapply(split_set_1, function(x) is.factor(x) || is.character(x))]
        # Identify which categorical variables have more than one unique level
        valid_categorical_vars <- sapply(categorical_vars, function(var) length(unique(split_set_1[[var]])) > 1)
        # Prepare a model formula excluding single-level categorical variables
        formula <- paste("owners ~", paste(names(valid_categorical_vars[valid_categorical_vars]), collapse = " + "), 
                         "+ year + price + mileage + engine.size + seats + 
                       fuel.economy + max.speed + engine.power")
        # Fit the Poisson model with the dynamic formula
        split_set_1_model <- glm(as.formula(formula), family=poisson(link="log"), data=split_set_1, control = control)
      }
      
      if (nrow(split_set_2)>0)
      {
        # Split Set 2
        # Identify columns that are factors or characters (potentially categorical)
        categorical_vars <- names(split_set_2)[sapply(split_set_2, function(x) is.factor(x) || is.character(x))]
        # Identify which categorical variables have more than one unique level
        valid_categorical_vars <- sapply(categorical_vars, function(var) length(unique(split_set_2[[var]])) > 1)
        # Prepare a model formula excluding single-level categorical variables
        formula <- paste("owners ~", paste(names(valid_categorical_vars[valid_categorical_vars]), collapse = " + "), 
                         "+ year + price + mileage + engine.size + seats + 
                       fuel.economy + max.speed + engine.power")
        split_set_2_model <- glm(as.formula(formula), family=poisson(link="log"), data=split_set_2, control = control)
      }
      # Find the reduction in node impurity
      # Node impurity: i(t) = n_t(mu_t - y.bar_t * log(mu_t))
      # Split set 1
      i_split_set_1 <- 0 # Initialise
      if(nrow(split_set_1)>0)
      {
        mu_hat_1 <- mean(predict(split_set_1_model, type = "response"))
        i_split_set_1 <- nrow(split_set_1)*(mu_hat_1-(mean(split_set_1$owners)*log(mu_hat_1)))
      }
      # Split set 2
      i_split_set_2 <- 0 # Initialise
      if(nrow(split_set_2)>0)
      {
        mu_hat_2 <- mean(predict(split_set_2_model, type = "response"))
        i_split_set_2 <- nrow(split_set_2)*(mu_hat_2-(mean(split_set_2$owners)*log(mu_hat_2)))
      }
      # Parent set
      mu_hat_data <- mean(predict(model, type = "response"))
      i_data <- nrow(data)*(mu_hat_data-(mean(data$owners)*log(mu_hat_data)))
      # Reduction in node impurity: i(t) - i(t_left) - i(t_right)
      i_reduction <- i_data - i_split_set_1 - i_split_set_2
      # Add to vector
      node_impurity_reduction <- c(node_impurity_reduction, i_reduction)
    }
  
    # Select the two split sets with the largest reduction in node impurity
    chosen_split_set_1 <- split_set_1_list[[which.max(node_impurity_reduction)]]
    chosen_split_set_2 <- split_set_2_list[[which.max(node_impurity_reduction)]]
    
    # Print the split point
    if (is.numeric(data[,split_var]))
    {
      print(split_points[[which.max(node_impurity_reduction)]])
    }else
    {
      print(subsets[[which.max(node_impurity_reduction)]])
    }
    
    # Print the number of rows in each split set
    print(nrow(chosen_split_set_1))
    print(nrow(chosen_split_set_2))
    
    # Recursive
    # Check if all covariates are homogenous
    # Does each covariate column contain one unique value
    is_column_unique <- function(column) {
      length(unique(column)) == 1
    }
    # Use function for all covariate columns (minus owners)
    all_columns_unique_1 <- all(sapply(chosen_split_set_1[,-1], is_column_unique))
    all_columns_unique_2 <- all(sapply(chosen_split_set_2[,-1], is_column_unique))
  }
  
  # If there is no split variable or one of the sets is empty, add the parent dataset as a leaf
  if (split_variable_is_zero || nrow(chosen_split_set_1)==0 || nrow(chosen_split_set_2)==0)
  {
    nodes_list <- append(nodes_list, list(data[, original_colnames]))
    print("Parent added as leaf")
  }else
  {
    # Split set 1
    # Continue with recursion if all covariates are homogenous
    if (!all_columns_unique_1)
    {
      # Continue with recursion if size of node is greater than 5% the size of total dataset
      if (nrow(chosen_split_set_1)>(0.05*nrow(original_data)))
      {
        nodes_list <- count_regression_tree(chosen_split_set_1[, original_colnames], original_data, nodes_list, iteration = iteration+1) # Recursive
        print("Full nodes list added - left")
      }else
      {
        nodes_list <- append(nodes_list, list(chosen_split_set_1[, original_colnames]))
        print("Added nodes smaller than0.05% of original data - left")
      }
    }else
    {
      nodes_list <- append(nodes_list, list(chosen_split_set_1[, original_colnames]))
      print("Added nodes that had all homogenous covariates - left")
    }
    
    
    # Split set 2
    # Continue with recursion if all covariates are homogenous
    if (!all_columns_unique_2)
    {
      # Continue with recursion if size of node is greater than 5% the size of total dataset
      if (nrow(chosen_split_set_2)>(0.05*nrow(original_data)))
      {
        nodes_list <- count_regression_tree(chosen_split_set_2[, original_colnames], original_data, nodes_list, iteration = iteration+1) # Recursive
        print("Full nodes list added - right")
      }else
      {
        nodes_list <- append(nodes_list, list(chosen_split_set_2[, original_colnames]))
        print("Added nodes smaller than0.05% of original data - right")
      }
    }else
    {
      nodes_list <- append(nodes_list, list(chosen_split_set_2[, original_colnames]))
      print("Added nodes that had all homogenous covariates - right")
    }
  }
  return(nodes_list)
}

# Run function
small_dataset_tree <- count_regression_tree(data_small, data_small)

#####################################################################################
# Finding the criteria for each node dataset (what covariate values are in each) 

# Initialise...
# These are the limits for each covariate for each leaf node dataset
# Change to a new leaf node every (number of covariates) iterations

find_covariate_limits <- function(tree)
{
  covariate_limits <- list() 
  # Calculate it for each node
  for (leaf_node in tree) 
  {
    # Loop through each covariate
    for (i in colnames(leaf_node)[-1]) 
    {
      if (is.numeric(leaf_node[,i])) # Two numbers (minimum and maximum for this covariate) if numeric
      {
        minimum <- min(leaf_node[,i]) # Minimum value for this covariate
        if (minimum == min(leaf_node[,i]))
        {
          # If the minimum value in this node is equal to the minimum value in the whole dataset
          # We assume all values lower than this min also enter this node
          minimum <- -Inf 
        }
        maximum <- max(leaf_node[,i])
        if (maximum == max(leaf_node[,i]))
        {
          # If the maximum value in this node is equal to the maximum value in the whole dataset
          # We assume all values higher than this min also enter this node
          maximum <- Inf 
        }
        # Add min and max to list as a vector
        covariate_limits <- append(covariate_limits, list(c(minimum, maximum)))
      }else # Categorical => multiple categories
      {
        # Add all unique categories in this leaf node to the list
        covariate_limits <- append(covariate_limits, list(unique(leaf_node[,i])))
      }
    }
  }
  return(covariate_limits)
}
# Run function
#covariate_limits <- find_covariate_limits(small_dataset_tree)

####################################################################################
# Classifying a new entry - return a list containing the branch the entry belongs to and the associated glm model
#observation <- data[2,]
# observation = new entry to 
# tree = list of leaf nodes of a CORE tree
# tree_limits = list of limits for each covariate of each leaf node
new_classification <- function(observation, tree, tree_limits)
{
  return_list <- list() # Initialise list for return
  # Loop through each potential class / leaf node
  for (i in 1:length(tree))
  {
    correct_classification <- TRUE # Boolean: leaf node i is correct leaf for this entry
    # Check for each covariate what group it should be in
    for (j in 1:(ncol(tree[[i]])-1)) 
    {
      if (correct_classification)
      {
        index <- (i-1)*(ncol(tree[[i]])-1)+(j) # Index of tree_limits for this covariate
        limits <- tree_limits[[index]] # Accepted limits for this covariate
        if (is.numeric(limits))
        {
          # Check if the value is in between the split points
          # j+1 used for column index as 'owners' column is avoided - not a covariate
          if ((observation[1,j+1]>=limits[1])&&(observation[1,j+1]<=limits[2])) 
          {
            # Continue the loop
          }else 
          {
            # Value is not in this leaf node
            correct_classification <- FALSE
          }
        }else
        {
          # Check if the value has only these covariates
          if (observation[1,j+1] %in% limits)
          {
            # Continue the loop
          }else 
          {
            # Value is not in this leaf node
            correct_classification <- FALSE
          }
        }
      }
    }
    if (correct_classification) # Select this leaf node
    {
      # Return this tree leaf
      classification_node <- tree[[i]]
      # Create model for selected leaf node
      categorical_covariates <- names(classification_node)[sapply(classification_node, function(x) is.factor(x) || is.character(x))]
      valid_categorical_covariates <- sapply(categorical_covariates, function(var) length(unique(classification_node[[var]])) > 1)
      covariates_for_modelling_formula <- paste("owners ~", paste(names(valid_categorical_covariates[valid_categorical_covariates]), collapse = " + "), 
                                                "+ year + price + mileage + engine.size + seats + 
                   fuel.economy + max.speed + engine.power")
      classification_node_poisson_model <- glm(as.formula(covariates_for_modelling_formula), 
                                               family = poisson(link = "log"), data = classification_node)
      # Return these objects
      return(list(classification_node, classification_node_poisson_model))
    }
  }
}

# Run function
#classification_list <- new_classification(observation, small_dataset_tree, covariate_limits)


#####################################################################################
# Run a model on this leaf node

# # Identify columns that are factors or characters (potentially categorical)
# categorical_covariates <- names(classification_list[[1]])[sapply(classification_list[[1]], function(x) is.factor(x) || is.character(x))]
# # Identify which categorical variables have more than one unique level
# valid_categorical_covariates <- sapply(categorical_covariates, function(var) length(unique((classification_list[[1]])[[var]])) > 1)
# # Prepare a model formula excluding single-level categorical variables
# covariates_for_modelling_formula <- paste("owners ~", paste(names(valid_categorical_covariates[valid_categorical_covariates]), collapse = " + "), 
#                  "+ year + price + mileage + engine.size + seats + 
#                    fuel.economy + max.speed + engine.power")
# # Fit the Poisson model with the dynamic formula
# classification_node_poisson_model <- glm(as.formula(covariates_for_modelling_formula), 
#                                          family = poisson(link = "log"), data = classification_list[[1]])
# 
# # Make a prediction
# new_observation <- data[1001,-1]
# predicted_value <- predict(classification_node_poisson_model, newdata = new_observation, type = "response")
# count_predicted_value <- round(predicted_value[[1]])

#####################################################################################
# Make a prediction using one Poisson model

# Training data (rows 1 to 1000)
control_training_data <- data[1:30000,]
# Create model
categorical_vars <- names(control_training_data)[sapply(control_training_data, function(x) is.factor(x) || is.character(x))]
valid_categorical_vars <- sapply(categorical_vars, function(var) length(unique(control_training_data[[var]])) > 1)
formula <- paste("owners ~", paste(names(valid_categorical_vars[valid_categorical_vars]), collapse = " + "), 
                 "+ year + price + mileage + engine.size + seats + 
                   fuel.economy + max.speed + engine.power")
control_poisson_model <- glm(as.formula(formula), family = poisson(link = "log"), data = control_training_data)

# Test data (rows 30001 to 40000)
control_test_data <- data[30001:40000,]
# Predict for all of these observations
control_predictions <- predict(control_poisson_model, newdata = control_test_data[,-1], type = "response")

#####################################################################################
# Make a prediction using CORE model

# Training data (rows 1 to 1000)
CORE_training_data <- data[1:30000,]
# Run CORE model and get split up datasets
CORE_tree <- count_regression_tree(CORE_training_data, CORE_training_data) # Takes really long
# Find the different limits for each covariate tree
CORE_limits <- find_covariate_limits(CORE_tree)
# Loop through each new observation to predict a value for each
control_test_data <- data[30001:40000,] # Test data
CORE_predictions <- c() # Initialise predictions vector
# Loop through each row
for (i in 1:nrow(control_test_data))
{
  # Find the correct leaf node for the data entry
  new_observation <- control_test_data[i,]
  new_leaf_node <- new_classification(new_observation, CORE_tree, CORE_limits)
  classification_node_poisson_model <- (new_classification(control_test_data[i,], CORE_tree, CORE_limits))[[2]]
  
  # Make a prediction
  predicted_value <- predict(classification_node_poisson_model, newdata = control_test_data[i,], type = "response")
  CORE_predictions <- c(CORE_predictions, predicted_value)
}
