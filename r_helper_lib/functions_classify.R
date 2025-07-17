# a binary thresholdclassifier
threshold_classifier <- function(predictor, threshold) {
  # Check if input predictor and threshold are valid
  if (!is.numeric(predictor) || !is.numeric(threshold)) {
    stop("Both predictor and threshold must be numeric values.")
  }
  # Apply the threshold rule to classify
  eps = 1e-15 # small value to avoid mismatch due to rounding
  classification <- ifelse(predictor >= threshold - eps, 1, 0)
  return(classification)
}

# the classifier function that takes *precomputed* NB values
MALT_classifier_lite <- function(arg_NB_val_low_H2rise, arg_threshold_low_H2rise, arg_NB_val_high_CH4rise, arg_threshold_high_CH4rise) {
    #Check if input predictor and threshold are valid
    #if (!is.numeric(predictor) || !is.numeric(threshold)) {
    #  stop("Both predictor and threshold must be numeric values.")
    #}
    if(arg_NB_val_high_CH4rise >= arg_threshold_high_CH4rise) {
        return(0)
    } else if(arg_NB_val_low_H2rise < arg_threshold_low_H2rise) {
        return(0)
    } else {
        # low CH4 and low H2 rise implies MALT
        return(1)
    }  
}

evaluate_classification = function(arg_ground_truth, arg_prediction, print_confusion_matrix = TRUE) {    
    arg_ground_truth = factor(arg_ground_truth, levels = c(0, 1))
    arg_prediction = factor(arg_prediction, levels = c(0, 1))
    confusion_matrix <- table(arg_ground_truth, arg_prediction)
    if(print_confusion_matrix) {
        print(confusion_matrix)
    }

    # Extract True Positives, True Negatives, False Positives, False Negatives
    TP <- confusion_matrix[2, 2]  # True Positives (A_actual = 1 and A_predicted = 1)
    TN <- confusion_matrix[1, 1]  # True Negatives (A_actual = 0 and A_predicted = 0)
    FP <- confusion_matrix[1, 2]  # False Positives (A_actual = 0 and A_predicted = 1)
    FN <- confusion_matrix[2, 1]  # False Negatives (A_actual = 1 and A_predicted = 0)

    # Calculate Sensitivity and Specificity
    recall = TP / (TP + FN) # aka sensitivity
    specificity = TN / (TN + FP)
    precision = TP / (TP + FP)
    f1_score = ifelse(precision == 0 | recall == 0, 0, # define F1 = 0 when precision or recall = 0
        2 * (precision * recall) / (precision + recall))
    g_mean = sqrt(specificity * recall)

    #return(lapply(list(recall = recall, specificity = specificity, precision = precision, f1_score = f1_score, g_mean = g_mean), function(x) {round(x, digits = 3)}))
    return(list(recall = recall, specificity = specificity, precision = precision, f1_score = f1_score, g_mean = g_mean))
}

do_auprc = function(arg_nb_vals, arg_classif_factor, output_text = TRUE, make_plot = TRUE) {
    prc_res = pr.curve(scores.class0 = arg_nb_vals, weights.class0 = arg_classif_factor, 
        curve = TRUE, rand.compute = TRUE, max.compute = TRUE, min.compute = TRUE)
    #print(prc_res)
    
    # Extract precision, recall, and threshold values from pr_results
    recall <- prc_res$curve[, 1]
    precision <- prc_res$curve[, 2]
    thresholds <- prc_res$curve[, 3]

    # better for imbalanced classes than integral method
    auprc = prc_res$auc.davis.goadrich
    auprc_random = prc_res$rand$auc.davis.goadrich        
    auprc_ratio = auprc / auprc_random

    # Calculate the F1 score for each threshold    
    f1_scores = ifelse(precision == 0 | recall == 0, 0,
        2 * (precision * recall) / (precision + recall))

    # Find the threshold that maximizes the F1 score
    optimal_index <- which.max(f1_scores)
    optimal_threshold_classif_factor <- median(thresholds[optimal_index])

    if(output_text) {
        # Print the optimal threshold and corresponding F1 score
        cat("\nAUPRC (Davis & Goadrich): ", auprc, "\n")
        cat("Random AUPRC (Davis & Goadrich): ", auprc_random, "\n")
        cat("AUPRC Ratio: ", auprc_ratio, "\n")
        cat("Optimal Threshold (arg max F1 score): ", optimal_threshold_classif_factor, "\n")
        cat("F1 Score at Optimal Threshold: ", f1_scores[optimal_index], "\n")
        cat("Recall at Optimal Threshold: ", recall[optimal_index], "\n")
        cat("Precision at Optimal Threshold: ", precision[optimal_index], "\n")    
    }

    # plot PRC curve
    if(make_plot) {
        p.dims(6, 5)
        plot(prc_res, rand.plot = TRUE, fill.area = TRUE, fill.color = rgb(0.8,1,0.8), maxminrand.col = "blue", auc.type = "davis.goadrich") 
    }

    c(  auprc = auprc, 
        auprc_random = auprc_random,
        auprc_ratio = auprc_ratio,
        opt_threshold = optimal_threshold_classif_factor,
        f1_opt = f1_scores[optimal_index],
        recall_opt = recall[optimal_index], 
        precision_opt = precision[optimal_index])
}


# Partial AUROC / Early Retrieval ! See details.
do_auroc = function(arg_nb_vals, arg_classif_factor, output_text = TRUE, make_plot = TRUE, metric_to_optimize = "g_mean") {
    # compute the number of positive cases    
    freq_pos = length(which(arg_classif_factor == 1))/length(arg_classif_factor)
    #print(freq_pos)
        
    roc_res = tryCatch({
        suppressWarnings({
            roc(arg_classif_factor, arg_nb_vals, quiet = TRUE,
            direction = "<",
            partial.auc = c(0.9, 1.0), # measure AUC for values [0.9; 1]
            partial.auc.focus = "spec", # ...of specificity -  that is [0.1; 0] for FPR
            partial.auc.correct = TRUE, 
            #allow.invalid.partial.auc.correct = FALSE # default way
            allow.invalid.partial.auc.correct = TRUE # tweak to allow invalid values instead of NA's to avoid introducing a bias against low pAUCs under small test sizes
        )
        })

    }, error = function(e) {
        message("Caught error: ", e$message)
        return(NULL)  # Return NULL instead of stopping execution
    })
    
    # Extract 
    recall = roc_res$sensitivities    
    spec = roc_res$specificities    
    thresholds = roc_res$thresholds

    # compute precision
    precision = recall * freq_pos / (recall * freq_pos + (1 - spec) * (1 - freq_pos))
    #print(precision)

    # Calculate the G-means for each threshold
    g_means = sqrt(spec * recall)
    #print(g_means)

    # Calculate the F1 score for each threshold, position wise    
    f1_scores = ifelse(is.na(precision) | is.na(recall) | precision == 0 | recall == 0, 0,
        2 * (precision * recall) / (precision + recall))    
    #print(f1_scores)

    # Find the threshold that maximizes the:
    # - G means
    if(metric_to_optimize == "g_mean") {
        optimal_index = which.max(g_means)
        #print(g_means[optimal_index])
    } else if(metric_to_optimize == "f1_score") {
        # - or F1 score
        optimal_index = which.max(f1_scores)
        #print(f1_scores[optimal_index])
    } else {
        stop("Invalid metric_to_optimize. Choose 'g_mean' or 'f1_score'.")
    }
    #print(optimal_index)

    optimal_threshold_classif_factor = median(thresholds[optimal_index])

    if(output_text) {
        # Print the optimal threshold and corresponding F1 score
        cat("\nAUROC: ", roc_res$auc, "\n")
        cat(paste0("Optimal Threshold (arg max " , metric_to_optimize, "): ", optimal_threshold_classif_factor, "\n"))
        cat("G-mean at Optimal Threshold: ", g_means[optimal_index], "\n")
        cat("F1 score at Optimal Threshold: ", f1_scores[optimal_index], "\n")        
        cat("Recall at Optimal Threshold: ", recall[optimal_index], "\n")
        cat("Specificity at Optimal Threshold: ", spec[optimal_index], "\n")    
    }

    # plot PRC curve
    if(make_plot) {
        p.dims(6, 5)
        #plot(roc_res, print.thres = "best", print.auc = TRUE, print.thres.best.method = "closest.topleft", xlab = "1 - Specificity")
        plot(roc_res, print.auc = TRUE, xlab = "1 - Specificity")
    }

    c(  auroc = as.numeric(roc_res$auc), 
        opt_threshold = optimal_threshold_classif_factor,
        f1_score_opt = f1_scores[optimal_index],
        g_mean_opt = g_means[optimal_index],
        recall_opt = recall[optimal_index], 
        spec_opt = spec[optimal_index])
}

# Ordinary AUROC
do_auroc_ordinary = function(arg_nb_vals, arg_classif_factor, output_text = TRUE, make_plot = TRUE, metric_to_optimize = "g_mean") {
    # compute the number of positive cases    
    freq_pos = length(which(arg_classif_factor == 1))/length(arg_classif_factor)
    #print(freq_pos)

    roc_res = tryCatch({
        suppressWarnings({
            roc(arg_classif_factor, arg_nb_vals, quiet = TRUE,
            direction = "<"
            )
        })

    }, error = function(e) {
        message("Caught error: ", e$message)
        return(NULL)  # Return NULL instead of stopping execution
    })

    # Extract 
    recall = roc_res$sensitivities    
    spec = roc_res$specificities    
    thresholds = roc_res$thresholds

    # compute precision
    precision = recall * freq_pos / (recall * freq_pos + (1 - spec) * (1 - freq_pos))
    #print(precision)

    # Calculate the G-means for each threshold
    g_means = sqrt(spec * recall)
    #print(g_means)

    # Calculate the F1 score for each threshold, position wise    
    f1_scores = ifelse(is.na(precision) | is.na(recall) | precision == 0 | recall == 0, 0,
        2 * (precision * recall) / (precision + recall))    
    #print(f1_scores)

    # Find the threshold that maximizes the:
    # - G means
    if(metric_to_optimize == "g_mean") {
        optimal_index = which.max(g_means)
        #print(g_means[optimal_index])
    } else if(metric_to_optimize == "f1_score") {
        # - or F1 score
        optimal_index = which.max(f1_scores)
        #print(f1_scores[optimal_index])
    } else {
        stop("Invalid metric_to_optimize. Choose 'g_mean' or 'f1_score'.")
    }
    #print(optimal_index)

    optimal_threshold_classif_factor = median(thresholds[optimal_index])

    if(output_text) {
        # Print the optimal threshold and corresponding F1 score
        cat("\nAUROC: ", roc_res$auc, "\n")
        cat(paste0("Optimal Threshold (arg max " , metric_to_optimize, "): ", optimal_threshold_classif_factor, "\n"))
        cat("G-mean at Optimal Threshold: ", g_means[optimal_index], "\n")
        cat("F1 score at Optimal Threshold: ", f1_scores[optimal_index], "\n")        
        cat("Recall at Optimal Threshold: ", recall[optimal_index], "\n")
        cat("Specificity at Optimal Threshold: ", spec[optimal_index], "\n")    
    }

    # plot PRC curve
    if(make_plot) {
        p.dims(6, 5)
        #plot(roc_res, print.thres = "best", print.auc = TRUE, print.thres.best.method = "closest.topleft", xlab = "1 - Specificity")
        plot(roc_res, print.auc = TRUE, xlab = "1 - Specificity")
    }

    c(  auroc = as.numeric(roc_res$auc), 
        opt_threshold = optimal_threshold_classif_factor,
        f1_score_opt = f1_scores[optimal_index],
        g_mean_opt = g_means[optimal_index],
        recall_opt = recall[optimal_index], 
        spec_opt = spec[optimal_index])
}