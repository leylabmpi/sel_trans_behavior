# collapse Bracken species-matrices at a higher taxonomy level (genus, family, order)
collapse_bracken_taxon = function(arg_br_species, arg_tax_level) {
    if(arg_tax_level == "Genus") {
        #res = arg_br_species %>% group_by(Genus, Sample)
        res = arg_br_species %>% group_by(Domain, Phylum, Class, Order, Family, Genus, Sample)        
    } else if(arg_tax_level == "Family") {
        #res = arg_br_species %>% group_by(Family, Sample)
        res = arg_br_species %>% group_by(Domain, Phylum, Class, Order, Family, Sample)        
    } else if(arg_tax_level == "Order") {
        #res = arg_br_species %>% group_by(Order, Sample)
        res = arg_br_species %>% group_by(Domain, Phylum, Class, Order, Sample)        
    } else {
        stop("Invalid taxonomic level. Please choose from 'Genus', 'Family', 'Order'.")
    }

    res = res %>% 
        mutate(Abundance = sum(Abundance)) %>%
        ungroup() %>%
        mutate(name = !!sym(arg_tax_level)) #%>%
        #arrange(Abundance) 
        #arrange(name, Sample, Abundance)

    if(arg_tax_level == "Genus") {
        res = res %>% select(-Species)
    } else if(arg_tax_level == "Family") {
        res = res %>% select(-Species, -Genus)
    } else if(arg_tax_level == "Order") {
        res = res %>% select(-Species, -Genus, -Family)
    }

    res %>% distinct()
}

# do collapse_bracken_taxon for all taxonomic levels
collapse_bracken_taxon_multiple = function(arg_br_species, arg_levels) {   
    res = map(arg_levels, 
        function(taxon) {
            collapse_bracken_taxon(arg_br_species, taxon)        
        }
    )
    res = purrr::set_names(res, arg_levels)
}

collapse_bracken_taxon_multiple_GFO = function(arg_br_species) {
    res_genus = arg_br_species %>%     
    #mutate(Abundance = sum(Abundance), .by = c(Sample, Genus)) %>% 
    mutate(Abundance = sum(Abundance), .by = c(Sample, Domain, Phylum, Class, Order, Family, Genus)) %>% 
    mutate(name = Genus) %>%    
    select(-Species) %>% 
    #arrange(name, Sample, Abundance)  %>% 
    distinct()

    res_family = res_genus %>%
    #mutate(Abundance = sum(Abundance), .by = c(Sample, Family)) %>% 
    mutate(Abundance = sum(Abundance), .by = c(Sample, Domain, Phylum, Class, Order, Family)) %>%     
    mutate(name = Family) %>%    
    select(-Genus) %>% 
    #arrange(name, Sample, Abundance)  %>% 
    distinct()

    res_order = res_family %>%
    #mutate(Abundance = sum(Abundance), .by = c(Sample, Order)) %>% 
    mutate(Abundance = sum(Abundance), .by = c(Sample, Domain, Phylum, Class, Order)) %>% 
    mutate(name = Order) %>%    
    select(-Family) %>% 
    #arrange(name, Sample, Abundance)  %>% 
    distinct()

    list(Genus = res_genus, Family = res_family, Order = res_order)

}

# filtering low prevalence at fixed abundance
filter_low_prev <- function(arg_bracken, arg_abund_cutoff_perc, arg_prev_cutoff) {        
    tmp_n_samples = arg_bracken %>% select(Sample) %>% distinct() %>% nrow()    
    res = arg_bracken %>%        
        mutate(prev = sum(Abundance > arg_abund_cutoff_perc) / tmp_n_samples * 100, 
                .by=c(name)) %>%
        filter(prev >= arg_prev_cutoff) %>%
        select(-prev)
    #print(res %>% unique_n('taxa', name))    
	return(res)
}

# check the long table is full (includes the same numbers of zeros as if it were wide) before computing the prevalence
check_full = function(x, arg_feature_col) {
    tmp_n_samples = x %>% select(Sample) %>% distinct() %>% nrow()
    tmp_n_feats = x %>% select(arg_feature_col) %>% distinct() %>% nrow()
    test_that("Check the table is full before computing the prevalence", {    
        expect_equal(tmp_n_samples * tmp_n_feats,  x %>% nrow())
    })
}

adjust_via_residuals = function(arg_meta_samples, arg_features_long, arg_factor_to_adj_for, arg_feature_col = "feature", arg_value_col = "value") {
    check_full(arg_features_long, arg_feature_col)

    feat_df = 
        arg_features_long %>%          
        select(!!sym(arg_feature_col), Sample, !!sym(arg_value_col)) %>% 
        pivot_wider(names_from = !!sym(arg_feature_col), values_from = !!sym(arg_value_col), values_fill = 0) %>% 
        column_to_rownames("Sample") %>% 
        as.data.frame

    meta_df_tmp = data.frame(arg_meta_samples %>% column_to_rownames("Sample"))
    meta_df_tmp = meta_df_tmp[rownames(feat_df),]

    test_that("feat_df and arg_meta_samples have the same names", {
        expect_equal(rownames(feat_df), rownames(meta_df_tmp))
    })

    # collect the residuals
    lm_res_ = lm(reformulate(arg_factor_to_adj_for, response = "as.matrix(feat_df)"), meta_df_tmp)

    feat_df_adj = data.frame(resid(lm_res_))
    # non-alphanum characters appear to be replaced with dots by lm, so restore the original names
    colnames(feat_df_adj) = colnames(feat_df)

    tmp = feat_df_adj %>% rownames_to_column("Sample") %>%
        pivot_longer(cols = -Sample, names_to = arg_feature_col, values_to = arg_value_col) 

    tmp = tmp %>% inner_join(arg_features_long %>% select(-arg_value_col), by = c("Sample", arg_feature_col))   
}

do_lm_tidy = function(arg_lm_in, arg_formula, arg_response_col = "value", arg_feature_col = "feature") {
    arg_lm_in %>%
        nest(data = -!!sym(arg_feature_col)) %>%    
        mutate(fit = map(data, function(x) {        
            lm1 = lm(reformulate(arg_formula, response = arg_response_col), data = x)
            return(lm1)
        }),
        tidied = map(fit, tidy)) %>%    
        select(feature, tidied) %>% 
        unnest(tidied) %>%
        filter(term != "(Intercept)") %>% 
        #filter(effect != "ran_pars") %>%            
        arrange(desc(term), p.value) #%>%
        #mutate(p.adj = p.adjust(p.value, method = "fdr"))   
}

do_lmer_tidy = function(arg_lm_in, arg_formula, arg_response_col = "value", arg_feature_col = "feature") {
    arg_lm_in %>%
        nest(data = -!!sym(arg_feature_col)) %>%    
        mutate(fit = map(data, function(x) {        
            lm1 = lmer(reformulate(arg_formula, response = arg_response_col), data = x)
            return(lm1)
        }),
        tidied = map(fit, tidy)) %>%    
        select(feature, tidied) %>% 
        unnest(tidied) %>%
        filter(term != "(Intercept)") %>% 
        filter(effect != "ran_pars") %>%            
        arrange(desc(term), p.value) #%>%
        #mutate(p.adj = p.adjust(p.value, method = "fdr"))   
}