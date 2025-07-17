# wrapper function to run 1 Nearest Balance analysis
make_nb_formula_to_1factor = function(arg_train, arg_test, arg_meta_train, arg_formula, arg_factor_coef)
{
	form = reformulate(arg_formula, response = "as.matrix(arg_train)")
	lm_res = lm(form, arg_meta_train)	
	coef_a = coefficients(lm_res)
	# convert rownames to alphanumeric
	rownames(coef_a) = gsub("[^[:alnum:]]", "_", rownames(coef_a))
	coef <- coef_a[arg_factor_coef,]
	print(coef)
	nb_curr <- find_nearest_balance_clr(coef)
	nb_curr
}


# inspired by groupKFold - but sampling without replace
leaveBatchOut = function (group, k = length(unique(group))) 
{
    g_unique <- unique(group)
    m <- length(g_unique)
    if (k > m) {
        stop("`k` should be less than ", m)
    }
    g_folds <- sample(k, size = m, replace = FALSE) #TRUE)
    out <- split(seq_along(group), g_folds[match(group, g_unique)])
    names(out) <- paste0("Fold", gsub(" ", "0", format(seq_along(out))))
    lapply(out, function(z) seq_along(group)[-z])
}

# leave N batches 
leaveNBatchesIn = function (arg_group, arg_n_batches_in, arg_n_splits) 
{	
	ub = arg_group %>% unique()	
	ub_samples = list()
	for(i in 1:arg_n_splits) {
		batches_subset <- sample(ub, size = arg_n_batches_in, replace = FALSE)
		cur_b = which(arg_group %in% batches_subset)	
		ub_samples[[paste0("Fold", i)]] = cur_b
	}
	ub_samples
}

#####

# run an LM model over the CLR abundance matrix and get the coefficients for the factor of choice
do_clr_lm_get_fac_coeffs = function(arg_lm_nb_formula, arg_clr, arg_meta_df, arg_sel_factor_coef)
{
	# LM over CLR data for our factor of interest
	lm_res = lm(reformulate(arg_lm_nb_formula, response = "as.matrix(arg_clr)"), arg_meta_df)
	#summary(lm_res)

	coef = coefficients(lm_res) # get all coefficients of LM
	# convert rownames to alphanumeric
	rownames(coef) = gsub("[^[:alnum:]]", "_", rownames(coef))
	print("Names of coefficients: ")
	print(rownames(coef))

	# get coefficients for the factor of choice
	data.frame(t(coef)) %>% 
		rownames_to_column(var = "name") %>%
		select(name, rlang::as_name(arg_sel_factor_coef)) %>%
		rename(lm_coef = rlang::as_name(arg_sel_factor_coef))
}


parallel_run_nb = function(arg_lm_nb_formula, arg_sel_factor_coef, arg_clr, arg_meta_df, arg_splits, arg_n_sim, arg_num_rparallel_cores) {
	# Initialize a parallel cluster with the specified number of cores
	cl = makeCluster(arg_num_rparallel_cores)
	registerDoParallel(cl)
	# Run the loop in parallel
	nb_res = foreach(i = 1:arg_n_sim, .packages = "NearestBalance", .export = "make_nb_formula_to_1factor") %dopar% {
		split = arg_splits[[i]]
		train_data = arg_clr[split,]
		test_data = arg_clr[-split,]
		meta_train = arg_meta_df[split,,drop=F]
		
		#make_nb_formula_to_1factor(train_data, test_data, meta_train, arg_lm_nb_formula, arg_sel_factor_coef)

		tryCatch({
			make_nb_formula_to_1factor(train_data, test_data, meta_train, arg_lm_nb_formula, arg_sel_factor_coef)
		}, error = function(e) NULL	)
	}

	# Stop the parallel cluster
	stopCluster(cl)

	# if there are nb_res entries, produce a warning and remove them		
	i_failed <- which(sapply(nb_res, is.null))	
	if(length(i_failed) > 0) {		
		nb_res = nb_res[-i_failed]
		print(paste0(length(i_failed), "/", arg_n_sim, " NB iterations results failed. Removed them."))
	}

	nb_res
}

flip_balance = function(arg_nb) {
	res = arg_nb
	res$b1$num = arg_nb$b1$den
	res$b1$den = arg_nb$b1$num
	res$sbp = -arg_nb$sbp        
	res
}

compute_balance = function(arg_abundance, arg_nb) {
	tibble(Sample = rownames(arg_abundance), NB_Value = as.vector(balance.fromSBP(arg_abundance, arg_nb$sbp)))
}

get_sbp_consensus_with_reprod = function(arg_sbp_iters) {
	arg_sbp_iters %>% select(taxName, b1 = b1_consensus, reprod)
}

# get consensus balance from multiple NB across iterations
aggregate_balance_iterations = function(arg_nb_res, arg_reproducibility_threshold) {
	# merge SBP from all nb_res entries via purrr
	sbp_iters = purrr::list_cbind(purrr::map(arg_nb_res, ~ .x$sbp), name_repair = "unique_quiet") %>% 
		dplyr::rename_all(~ gsub("b1...", "b1_", .))	
	rownames_taxo = sbp_iters %>% rownames_to_column(var = "taxName") %>% select(taxName)
	# for each row of sbp_iters, count the occurrence of "-1" and "1" values.
	# and set b1_consensus to 1 if freq_num > reproducibility_threshold, to -1 if freq_den > reproducibility_threshold, to 0 otherwise
	sbp_iters = sbp_iters %>% mutate(freq_den = rowSums(sbp_iters == -1)/n_sim,
		freq_num = rowSums(sbp_iters == 1)/n_sim) %>%
		mutate(b1_consensus = ifelse(freq_num >= arg_reproducibility_threshold, 1, ifelse(freq_den >= arg_reproducibility_threshold, -1, 0)),
		reprod = pmax(freq_num, freq_den))

	#sbp_iters = tidytable(cbind(rownames_taxo, sbp_iters))
	sbp_iters = tibble(cbind(rownames_taxo, sbp_iters))

	# and get a data frame version alike to  NB sbp format 	
	sbp_consensus = get_sbp_consensus_with_reprod(sbp_iters) %>% select(-reprod) %>% column_to_rownames("taxName")

	list(sbp_iters = sbp_iters, sbp_consensus = sbp_consensus)
}

# bring the consensus NB into a format similar to those of the produced by find_nearest_balance_clr
format_consensus_balance = function(arg_sbp_consensus) {
	res = list()
	res$b1 = list()	
	res$b1$num = rownames(arg_sbp_consensus)[which(arg_sbp_consensus$b1 == 1)]	
	res$b1$den = rownames(arg_sbp_consensus)[which(arg_sbp_consensus$b1 == -1)]	
	res$sbp = arg_sbp_consensus
	res
}

# join sbp_consensus_reprod with coef_oa_bacs
join_coefs_and_sbp_reprod = function(arg_sbp_consensus_reprod, arg_coef_oa_bacs) {
	arg_sbp_consensus_reprod %>% inner_join(arg_coef_oa_bacs, by = c("taxName" = "name")) %>%
		filter(b1 %in% c(1, -1)) #%>%
		#mutate(sign_lm_coef = factor(b1))
}

# balance size
balance_size = function(arg_nb) {
	list(n_num = arg_nb$b1$num %>% length(), n_den = arg_nb$b1$den %>% length())
}

# Perform the inverse clr transformation
inverse_clr = function(arg_clr) {
  exp_clr = exp(arg_clr)
  res_abundance = exp_clr / rowSums(exp_clr)  # Normalize rows
  return(res_abundance)
}

# get the mean abundance per group (i.e. Country)
get_clr_means_per_group = function(arg_clr_my, arg_meta, arg_factor, arg_join_meta_by = c('Sample' = 'ParticipantID')) {
    arg_clr_my %>% rownames_to_column('Sample') %>% 
        inner_join(arg_meta %>% select(ParticipantID, !!sym(arg_factor)), by = arg_join_meta_by) %>%
        pivot_longer(cols = -c(Sample, !!sym(arg_factor)), names_to = "Taxon", values_to = "Abundance") %>% 
        group_by(!!sym(arg_factor), Taxon) %>%
        mutate(mean_per_group = mean(Abundance))  %>% 
        ungroup() %>% 
        select(!!sym(arg_factor), Taxon, mean_per_group) %>% 
        distinct() %>% 
        arrange(!!sym(arg_factor), Taxon)
}

# subtract the precomputed mean abundance per group (i.e. Country) from the clr matrix
subtract_clr_means_per_group = function(arg_clr_my, arg_clr_means, arg_meta, arg_factor, arg_join_meta_by = c('Sample' = 'ParticipantID')) {
    res = arg_clr_my %>% rownames_to_column('Sample') %>% 
        inner_join(arg_meta %>% select(ParticipantID, !!sym(arg_factor)), by = arg_join_meta_by) %>%
        pivot_longer(cols = -c(Sample, !!sym(arg_factor)), names_to = "Taxon", values_to = "Abundance") %>% 
        inner_join(arg_clr_means, by = c(arg_factor, "Taxon")) %>%        
        mutate(Abundance = Abundance - mean_per_group) %>% 
        select(-!!sym(arg_factor), -mean_per_group) %>% 
        pivot_wider(names_from = Taxon, values_from = Abundance) %>% 		
        column_to_rownames('Sample') %>%     
        as.data.frame
	res = res[, sort(colnames(res))]	
	res
}

# exclude a set of taxa from a balance
exclude_taxa_from_balance = function(arg_nb, arg_taxa_list) {    
    res_nb = arg_nb
    if(length(arg_taxa_list) > 0) {
        res_nb$b1$num = setdiff(arg_nb$b1$num, arg_taxa_list)
        res_nb$b1$den = setdiff(arg_nb$b1$den, arg_taxa_list)
        res_nb$sbp = arg_nb$sbp[!rownames(arg_nb$sbp) %in% arg_taxa_list, , drop=FALSE]
    }
    res_nb
}