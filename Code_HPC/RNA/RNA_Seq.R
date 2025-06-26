.libPaths( c("~/Rlibs/4.3", .libPaths()) )

library(dplyr)
library(randomForest)
library(caret)
library(DESeq2)
library(ggplot2)
library(tidyverse)

cnt <- read.delim("/work/users/w/x/wxueyao/Furey_Lab/RNA/plexus-RNA-Seq_ruvg_normalized_counts_k3.tsv", header=T, check.names = FALSE)
met <- read.delim("/work/users/w/x/wxueyao/Furey_Lab/RNA/coldata_for_plexus-RNA-Seq_ruvg_normalized_counts.tsv", header=T)

# Rename DEIDENTIFIED_MASTER_PATIENT_ID as ID and set it to chr variable
mdata <- met %>%
  select(ID, Group) %>%
  mutate(ID = as.character(ID))
mdata$Group <- as.character(gsub("[^[:alnum:]_.]", "_", mdata$Group))

# Adding row name
rownames(cnt) <- cnt[, 1]
cdata <- cnt[, -1]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: ATAC_Seq.R <seed> [<iteration_number>]")
}

seed_val <- as.integer(args[1])
if (length(args) >= 2) {
  iter_num <- as.integer(args[2])
} else {
  iter_num <- seed_val
}

set.seed(seed_val)

n <- nrow(mdata)

# Generate fold indices
folds <- cut(sample(1:n), breaks = 5, labels = FALSE)
folds

# Add the folds to metadata and count data
mdata$folds <- folds

mdata$Group <- as.factor(mdata$Group)

model_summary <- list()
misclassified <- list()
test_accuracy <- list()
test_sensitivity <- list()
test_specificity <- list()
importance_df <- list()
feature_count <- list()
  
# Perform 5-fold cross-validation
for (i in 1:5) {
  
  test_index <- which(mdata$folds == i, arr.ind = TRUE)  # Get the indices for the test set
  train_mdata <- mdata[-test_index, ]  # Training set: all rows except the test set
  test_mdata  <- mdata[test_index, ] # Test set: only rows in the current fold
  
  ID_all <- colnames(cdata)
  
  # Find the columns in second dataset that match values in ID column from the first dataset
  match <- ID_all[ID_all %in% train_mdata$ID]
  
  # Select columns based on 'match'
  train_cdata <- cdata %>%
    select(all_of(match))
  
  
  # Select columns NOT in 'match'
  test_cdata <- cdata %>%
    select(-all_of(match))
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData = train_cdata,
                                colData = train_mdata,
                                design = ~ Group)
  dds <- DESeq(dds)
  
  # View result
  res <- results(dds)
  
  # Extract the top 1000 features
  # Subset the dataframe to include only rows where baseMean >= 100
  res_filtered <- res[res$baseMean >= 100, ]
  
  # Check if the filtered dataframe is not empty
  if (nrow(res_filtered) > 0) {
    # Extract padj and calculate rank
    padj <- res_filtered[, 6]
    rankpadj <- rank(-padj, na.last = "keep", ties.method = "average")
    
    # Add RankPAdj to the dataframe
    res_filtered$RankPAdj <- rankpadj
    
    # Order by RankPAdj and select top 1000 features
    orderP <- res_filtered[order(-res_filtered$RankPAdj), ]
    fea_select <- head(orderP, 1000)
  } else {
    # Handle case where no rows meet the condition
    fea_select <- NULL
    message("No rows with baseMean >= 100")
  }
  
  # Train set reshape
  # Extract selected features
  fea <- rownames(fea_select)
  
  # Transpose the count data matrix
  trans_train_cdata <- data.frame(t(train_cdata))
  
  # Keep the features we want
  train_ana <- trans_train_cdata %>% 
    select(all_of(fea))
  train_ana$ID <- rownames(train_ana)
  
  # Left join with the metadata
  train_ana <- left_join(train_ana, mdata, by = "ID")
  rownames(train_ana) <- train_ana$ID
  train_ana <- train_ana %>% select(-ID)
  
  # Remove folds.x, folds.y and folds, in there are any
  if ("folds.x" %in% colnames(train_ana)) {
    train_ana <- train_ana %>% select(-folds.x)
  }
  if ("folds.y" %in% colnames(train_ana)) {
    train_ana <- train_ana %>% select(-folds.y)
  }
  if ("folds" %in% colnames(train_ana)) {
    train_ana <- train_ana %>% select(-folds)
  }
  
  # Make sure the format of the variables are correct
  cols_train <- setdiff(names(train_ana), "Group")
  train_ana[cols_train] <- lapply(train_ana[cols_train], function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }
    as.numeric(x)
  })
  
  # Test set reshape
  # Transpose the count data matrix
  trans_test_cdata <- data.frame(t(test_cdata))
  
  # Keep the features we want
  test_ana <- trans_test_cdata %>% 
    select(all_of(fea))
  test_ana$ID <- rownames(test_ana)
  
  # Left join with the metadata
  test_ana <- left_join(test_ana, mdata, by = "ID")
  rownames(test_ana) <- test_ana$ID
  test_ana <- test_ana %>% select(-ID)
  
  # Remove folds.x, folds.y and folds, in there are any
  if ("folds.x" %in% colnames(test_ana)) {
    test_ana <- test_ana %>% select(-folds.x)
  }
  if ("folds.y" %in% colnames(test_ana)) {
    test_ana <- test_ana %>% select(-folds.y)
  }
  if ("folds" %in% colnames(test_ana)) {
    test_ana <- test_ana %>% select(-folds)
  }
  
  # Make sure the format of the variables are correct
  cols_test <- setdiff(names(test_ana), "Group")
  test_ana[cols_test] <- lapply(test_ana[cols_test], function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }
    as.numeric(x)
  })
  
  # Run Random Forest model
  rf <- randomForest(as.factor(Group) ~ ., data = train_ana, importance = TRUE)
  
  # Predictions
  # predictions_train <- predict(rf, newdata = train_ana)
  
  # Model Summary
  cat(paste0("\n\nVariable Ranking result for Loop ", i, "\n"))
  print(rf)
  
  # Misclassified Individuals
  #misclassified <- test_mdata[test_ana$Group != predictions, ]
  #cat("\nMisclassified Individuals:\n")
  #print(misclassified)
  
  # Accuracy for Train Dataset
  #confusion_matrix_train <- table(predictions_train, train_ana$Group)
  #accuracy_train <- sum(diag(confusion_matrix_train)) / sum(confusion_matrix_train)
  #cat(paste("Accuracy with Train Dataset:", round(accuracy_train * 100, 4), "%\n"))
  
  # Accuracy for Test Dataset
  #confusion_matrix <- table(predictions, test_ana$Group)
  #accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  #cat(paste("Accuracy with Test Dataset:", round(accuracy * 100, 4), "%\n"))
  
  # Variable Importance
  importance_df <- data.frame(importance(rf))
  importance_df$RankAccuracy <- rank(importance_df$MeanDecreaseAccuracy, na.last = F, ties.method = "average")
  importance_df$RankGini <- rank(importance_df$MeanDecreaseGini, na.last = F, ties.method = "average")
  importance_df$Score <- importance_df$RankAccuracy + importance_df$RankGini
  importance_df <- importance_df[order(-importance_df$Score), ]
  importance_df$Feature <- rownames(importance_df)
  
  # Plot Variable Importance
  #ggplot(importance_df, aes(x = seq_along(Score), y = Score)) +
  #geom_point(color = "blue") +
  #geom_line(color = "blue") +
  #labs(title = paste("Importance Score Plot (Fold", i, ")"),
  #x = "Feature Index",
  #y = "Importance Score") +
  #theme_minimal()
  
  # Plot Score
  #plot(importance_df$Score, main = paste("Importance Score Plot (Fold", i, ")"), ylab = "Score", xlab = "Features", type = "b")
  
  # Store results in list for further analysis if needed
  #results[[i]] <- list(
  #model_summary = rf,
  #misclassified = misclassified,
  #accuracy_train = accuracy_train,
  #accuracy_test = accuracy,
  #importance_df = importance_df
  #)
  
  ### Step 2: Select optimal number of features
  
  fea_score <- list()
  fea_err <- data.frame(numfea = numeric(), error = numeric())
  
  for (j in 1:100) {
    
    model_compare <- list()
    
    importance_mod <- importance_df
    train_mod <- train_ana
    num_fea <- 1000
    
    for (a in 1:10) {
      # Determine features used
      numfea <- ceiling(nrow(importance_mod) / 2)
      var_select <- importance_mod$Feature[1:numfea]
      
      # Preparing training count data for RF
      cdata_fil <- train_cdata[rownames(train_cdata) %in% var_select, ]
      
      # Feature Selection for train set
      # Keep the features we want
      train_mod <- data.frame(t(cdata_fil))
      train_mod$ID <- rownames(train_mod)
      
      # Left join with the metadata
      train_mod <- left_join(train_mod, mdata, by = "ID")
      train_mod <- train_mod %>% select(-ID)
      
      # Remove folds.x, folds.y and folds, in there are any
      if ("folds.x" %in% colnames(train_mod)) {
        train_mod <- train_mod %>% select(-folds.x)
      }
      if ("folds.y" %in% colnames(train_mod)) {
        train_mod <- train_mod %>% select(-folds.y)
      }
      if ("folds" %in% colnames(train_mod)) {
        train_mod <- train_mod %>% select(-folds)
      }
      
      # Make sure the format of the variables are correct
      cols_train_mod <- setdiff(names(train_mod), "Group")
      train_mod[cols_train_mod] <- lapply(train_mod[cols_train_mod], function(x) {
        if (is.factor(x)) {
          x <- as.character(x)
        }
        as.numeric(x)
      })
      
      # Run Random Forest model
      rf <- randomForest(as.factor(Group) ~ ., data = train_mod, importance = TRUE)
      
      # Variable Importance
      importance_mod <- data.frame(importance(rf))
      importance_mod$RankAccuracy <- rank(importance_mod$MeanDecreaseAccuracy, na.last = F, ties.method = "average")
      importance_mod$RankGini <- rank(importance_mod$MeanDecreaseGini, na.last = F, ties.method = "average")
      importance_mod$Score <- importance_mod$RankAccuracy + importance_mod$RankGini
      importance_mod <- importance_mod[order(-importance_mod$Score), ]
      importance_mod$Feature <- rownames(importance_mod)
      
      error <- rf$err.rate[nrow(rf$err.rate), "OOB"]
      
      importance_mod$numfea <- numfea
      importance_mod$error <- error
      
      model_compare[[a]] <- list(
        numfea = numfea,
        error = error,
        model_summary = rf,
        importance = importance_mod
      )
    }
    
    #for (d in 1:10) {
    #print(model_compare[[d]][1])
    #print(model_compare[[d]][2])
    #}
    
    for (b in 1:10) {
      # Extract values from model_compare
      numfea_value <- model_compare[[b]][1]
      error_value <- model_compare[[b]][2]
      
      # Assign values row-wise
      fea_err[(j-1)*10+b, ] <- c(numfea_value, error_value)
    }
    
    fea_score[[j]] <- list(
      importance_500 = model_compare[[1]][4],
      importance_250 = model_compare[[2]][4],
      importance_125 = model_compare[[3]][4],
      importance_63 = model_compare[[4]][4],
      importance_32 = model_compare[[5]][4],
      importance_16 = model_compare[[6]][4],
      importance_8 = model_compare[[7]][4],
      importance_4 = model_compare[[8]][4],
      importance_2 = model_compare[[9]][4],
      importance_1 = model_compare[[10]][4]
    )
  }
  
  mean_error <- fea_err %>%
    group_by(numfea) %>%
    summarise(mean_error = mean(error, na.rm = TRUE)) %>%
    arrange(desc(numfea)) %>%
    ungroup()
  
  # Create scatter plot with a smooth trend line
  ggplot(mean_error, aes(x = numfea, y = mean_error)) +
    geom_point(color = "blue", size = 3) +  # Scatter points
    geom_line(color = "blue", linetype = "dashed") +  # Connect points
    labs(title = "Number of Features vs Mean OOB Error Rate",
         x = "Number of Features",
         y = "Mean OOB Error Rate") +
    theme_minimal()
  
  # Check the final data frame
  min_err <- mean_error[which.min(mean_error$mean_error), ]
  min <- min(min_err$numfea)
  
  # Initialize an empty vector to store all row names from importance_16
  all_rownames <- c()
  
  # Loop through all j elements in fea_score and extract row names from importance_16
  for (e in seq_along(fea_score)) {
    col_name <- paste0("importance_", min)  # Correct column name
    if (col_name %in% names(fea_score[[e]])) {  
      all_rownames <- c(all_rownames, rownames(data.frame(fea_score[[e]][[col_name]])))
    }
  }
  
  # Create a frequency table and convert it into a data frame
  rownames_count <- as.data.frame(table(all_rownames), stringsAsFactors = FALSE)
  
  # Ensure correct column names
  colnames(rownames_count) <- c("Features", "Count")
  
  # Sort in descending order of count
  rownames_count <- rownames_count[order(-rownames_count$Count), ]
  
  # Print result
  rownames_count
  
  # Print number of features
  min_err
  
  # Try a full model
  
  # Create train and test datasets based on the features selected
  train_final_full <- train_ana %>% 
    select(any_of(rownames_count$Features), Group)
  
  test_final_full <- test_ana %>% 
    select(any_of(rownames_count$Features), Group)
  
  rf_final_full <- randomForest(as.factor(Group) ~ ., data = train_final_full, importance = TRUE)
  
  # Model Summary
  print(rf_final_full)
  
  names(train_final_full) <- make.names(names(train_final_full))
  names(test_final_full) <- make.names(names(test_final_full))
  
  # Variable Importance
  importance_final_full <- importance(rf_final_full)
  
  # Get variable importance from the model
  var_imp <- importance(rf_final_full)
  
  # Convert to a dataframe
  var_imp_df <- data.frame(
    Feature = rownames(var_imp),
    MDA = var_imp[, "MeanDecreaseAccuracy"],
    MDG = var_imp[, "MeanDecreaseGini"]
  )
  
  # Plot Mean Decrease in Accuracy (MDA)
  #ggplot(var_imp_df, aes(x = reorder(Feature, MDA), y = MDA)) +
  #geom_bar(stat = "identity", fill = "steelblue") +
  #coord_flip() +
  #labs(title = "Variable Importance: Mean Decrease in Accuracy (MDA)",
  #x = "Features", y = "MDA") +
  #theme_minimal()
  
  # Plot Mean Decrease in Gini (MDG)
  #ggplot(var_imp_df, aes(x = reorder(Feature, MDG), y = MDG)) +
  #geom_bar(stat = "identity", fill = "darkred") +
  #coord_flip() +
  #labs(title = "Variable Importance: Mean Decrease in Gini (MDG)",
  #x = "Features", y = "MDG") +
  #theme_minimal()
  
  # Predictions on the training set
  train_predictions <- predict(rf_final_full, train_final_full)
  
  # Predictions on the test set
  test_predictions <- predict(rf_final_full, test_final_full)
  
  # Add the predictions as a new column to the test set
  test_results_full <- test_final_full %>% 
    mutate(Predicted = test_predictions)
  
   # Filter for misclassified individuals (where the true label differs from the prediction)
  mis_i <- test_results_full %>% 
    filter(Group != Predicted)
  
  # Compute confusion matrices
  train_conf_matrix <- confusionMatrix(train_predictions, train_final_full$Group)
  test_conf_matrix <- confusionMatrix(test_predictions, test_final_full$Group)
  
  # Extract accuracy
  train_accuracy <- train_conf_matrix$overall['Accuracy']
  test_accuracy <- test_conf_matrix$overall['Accuracy']
  
  # Print results
  cat("Train Accuracy for Loop ",i, ": ", formatC(train_accuracy, format = "f", digits = 4), "\n")
  cat("Test Accuracy for Loop ",i, ": ", formatC(test_accuracy, format = "f", digits = 4), "\n")
  
  # Accuracy, Sensitivity, and Specificity
  train_sensitivity <- train_conf_matrix$byClass['Sensitivity']
  test_sensitivity <- test_conf_matrix$byClass['Sensitivity']
  train_specificity <- train_conf_matrix$byClass['Specificity']
  test_specificity <- test_conf_matrix$byClass['Specificity']
  
  # Print Metrics
  cat("Train Sensitivity for Loop ",i, ": ", train_sensitivity, "\n")
  cat("Test Sensitivity for Loop ",i, ": ", test_sensitivity, "\n")
  cat("Train Specificity for Loop ",i, ": ", train_specificity, "\n")
  cat("Test Specificity for Loop ",i, ": ", test_specificity, "\n")
  
  # Store results in list for further analysis if needed
  model_summary[[i]] <- list(Model = rf_final_full)
  misclassified[[i]] <- list(Misclassified = mis_i)
  test_accuracy[[i]] <- list(Test_Accuracy = test_accuracy)
  test_sensitivity[[i]] <- list(Test_Sensitivity = test_sensitivity)
  test_specificity[[i]] <- list(Test_Specificity = test_specificity)
  importance_df[[i]] <- list(Importance = var_imp_df)
  feature_count[[i]] <- list(Feature_count = rownames_count)
}

# Extract all row names across the list and unlist
rownames_list <- lapply(misclassified, function(x) {
  if (length(x) >= 1 && is.data.frame(x[[1]]) && nrow(x[[1]]) > 0) {
    rownames(x[[1]])
  } else {
    character(0)
  }
})

mis_rows <- unlist(rownames_list, use.names = FALSE)

# Count how many times each appears
mis_counts <- as.data.frame(table(mis_rows))

# Optionally, sort in descending order
mis_counts <- mis_counts[order(-mis_counts$Freq), ]

# Rename columns for clarity
colnames(mis_counts) <- c("Features", "Count")

# Convert the lists to numeric vectors
accuracy_vector    <- unlist(lapply(test_accuracy, function(x) test_accuracy))
sensitivity_vector <- unlist(lapply(test_sensitivity, function(x) test_sensitivity))
specificity_vector <- unlist(lapply(test_specificity, function(x) test_specificity))

# Calculate averages for each metric
average_accuracy    <- mean(accuracy_vector, na.rm = TRUE)
average_sensitivity <- mean(sensitivity_vector, na.rm = TRUE)
average_specificity <- mean(specificity_vector, na.rm = TRUE)

# Combine averages into a final dataframe
final_metrics <- data.frame(
  Average_Accuracy    = average_accuracy,
  Average_Sensitivity = average_sensitivity,
  Average_Specificity = average_specificity
)

# Extract and combine the dataframes
all_features <- do.call(rbind, lapply(feature_count, function(x) x$Feature_count))

# Use aggregate to group by Feature and sum Count
combined_features <- aggregate(Count ~ Features, data = all_features, FUN = sum)

# Optionally, rename the column for clarity
names(combined_features)[2] <- "Total_Count"

# -------------------------------
# Output the datasets
# -------------------------------

# Write combined_features to a CSV file with an iteration-specific name
output_combined <- paste0("combined_features_", iter_num, ".csv")
write.csv(combined_features, file = output_combined, row.names = FALSE)

# Write final_metrics to a CSV file with an iteration-specific name
output_metrics <- paste0("final_metrics_", iter_num, ".csv")
write.csv(final_metrics, file = output_metrics, row.names = FALSE)

# Write mis_counts to a CSV file with an iteration-specific name
output_mis <- paste0("mis_counts_", iter_num, ".csv")
write.csv(mis_counts, file = output_mis, row.names = FALSE)
