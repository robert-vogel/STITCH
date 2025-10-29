
mk_model_data <- function() {}


mk_bgen_model <- function(
    sampleRange,
    B_bit_prob,
    allAlphaBetaBlocks,
    allPhasing,
    alphaMatCurrentLocal_tc,
    eHapsCurrent_tc,
    transMatRateLocal_tc_H,
    transMatRateLocal_tc_D,    
    first_grid_in_region,
    last_grid_in_region,    
    i_output_block,
    read_starts_and_ends,
    first_snp_in_region,
    last_snp_in_region,
    nSNPsInOutputBlock,
    output_format,
    method,
    grid,
    highCovInLow,
    allSampleReads,
    outputdir,
    niterations,
    maxEmissionMatrixDifference,
    maxDifferenceBetweenReads,
    Jmax,
    useTempdirWhileWriting,
    output_haplotype_dosages,
    do_phasing) {

    bundledSampleReads <- NULL
    bundledSampleProbs <- NULL
    bundledAlphaBetaBlocks <- NULL
    bundledPhasing <- NULL
    
    ## load sample
    ## load pRgivenH1
    K <- dim(eHapsCurrent_tc)[1]        
    nSNPs <- dim(eHapsCurrent_tc)[2]
    S <- dim(eHapsCurrent_tc)[3]
    who_to_run <- sampleRange[1]:sampleRange[2]
    N_core <- sampleRange[2] - sampleRange[1] + 1 ## number in this core
    hweCount <- array(0, c(nSNPsInOutputBlock, 3))
    infoCount <- array(0, c(nSNPsInOutputBlock, 2))
    afCount <- array(0, nSNPsInOutputBlock)

    gp_raw_t <- matrix(
            data = raw(0),
            nrow = N_core * 2 * (B_bit_prob / 8),
            ncol = nSNPsInOutputBlock
        )
        vcf_matrix_to_out <- NULL
    
    pRgivenH1_m<- NULL
    pRgivenH2_m<- NULL


    
    # Set model  
    if (first_grid_in_region == last_grid_in_region) {
        impute_vals <- make_fbsoL_for_single_grid
    else
        impute_vals <- run_forward_backwards


    ## Per sample
    return(function(sample_name) {
        ## get these from list?
        out <- get_sampleReads_from_dir_for_sample(
            dir = tempdir,
            regionName = regionName,
            iSample = sample_name,
            bundling_info = bundling_info,
            bundledSampleReads = bundledSampleReads,
            allSampleReads = allSampleReads
        )
        sampleReads <- out$sampleReads
        bundledSampleReads <- out$bundledSampleReads
        
        out <- get_alphaBetaBlocks_from_dir_for_sample(
            dir = tempdir,
            regionName = regionName,
            iSample = sample_name,
            bundling_info = bundling_info,
            bundledAlphaBetaBlocks = bundledAlphaBetaBlocks,
            allAlphaBetaBlocks = allAlphaBetaBlocks
        )
        alphaBetaBlocks <- out$alphaBetaBlocks
        bundledAlphaBetaBlocks <- out$bundledAlphaBetaBlocks
        
        ## note, inefficient in current form
        if (do_phasing) {
            out <- get_phasing_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = bundling_info,
                bundledPhasing = bundledPhasing,
                allPhasing = allPhasing
            )
            phasing <- out$phasing
            bundledPhasing <- out$bundledPhasing
        }
        
        if (method == "pseudoHaploid") {
            out <- get_sampleProbs_from_dir_for_sample(
                dir = tempdir,
                regionName = regionName,
                iSample = iSample,
                bundling_info = bundling_info,
                bundledSampleProbs = bundledSampleProbs
            )
            pRgivenH1_m <- out$pRgivenH1_m
            pRgivenH2_m <- out$pRgivenH2_m
            srp <- out$srp
            bundledSampleProbs <- out$bundledSampleProbs
        }
        
        ## run forward backwards
        s <- read_starts_and_ends[iSample, i_output_block, 1]
        e <- read_starts_and_ends[iSample, i_output_block, 2]
        which_reads <- s:e
        if (is.na(s) | is.na(e))
            which_reads <- NULL
    

        fbsoL <- impute_genotypes()

        return(list(gp_raw_t = gp_raw_t,
                vcf_matrix_to_out = vcf_matrix_to_out,
                infoCount = infoCount,
                hweCount = hweCount,
                afCount = afCount))
    })
}



mk_bgvcf_model <- function() {}


mk_bcf_model <- function() {
    vcf_matrix_to_out <- data.frame(matrix(
        data = NA,
        nrow = nSNPsInOutputBlock,
        ncol = N_core
    ))
}





impute_genotypes <- function(sample_name) {

    if (first_grid_in_region == last_grid_in_region) {
    
        fbsoL <- make_fbsoL_for_single_grid(
            alphaBetaBlocks = alphaBetaBlocks,
            nSNPsInOutputBlock = nSNPsInOutputBlock,
            S = S,
            K = K,
            eHapsCurrent_tc = eHapsCurrent_tc,
            grid = grid,
            first_snp_in_region = first_snp_in_region,
            last_snp_in_region = last_snp_in_region,
            first_grid_in_region = first_grid_in_region,
            i_output_block = i_output_block,
            method = method
        )
        
    } else {    
    
        fbsoL <- run_forward_backwards(
            sampleReads = sampleReads[which_reads],
            pRgivenH1_m = pRgivenH1_m[which_reads, , drop = FALSE],
            pRgivenH2_m = pRgivenH2_m[which_reads, , drop = FALSE],
            method = method,
            priorCurrent_m = array(-1, c(K, S)), ## irrelevant here
            alphaMatCurrent_tc = alphaMatCurrentLocal_tc,
            eHapsCurrent_tc = eHapsCurrent_tc,
            transMatRate_tc_H = transMatRateLocal_tc_H,
            transMatRate_tc_D = transMatRateLocal_tc_D, 
            list_of_alphaBetaBlocks = alphaBetaBlocks, ## yes this is the right input
            i_snp_block_for_alpha_beta = i_output_block,
            run_fb_grid_offset = first_grid_in_region,
            run_fb_subset = TRUE,
            Jmax = Jmax,
            maxDifferenceBetweenReads = maxDifferenceBetweenReads,
            maxEmissionMatrixDifference = maxEmissionMatrixDifference,
            niterations = niterations,
            iteration = niterations,
            return_genProbs = TRUE,
            grid = grid,
            snp_start_1_based = first_snp_in_region,
            snp_end_1_based = last_snp_in_region,
            output_haplotype_dosages = output_haplotype_dosages ## whether to return states
        )
        
    }
    
    
    gp_t <- calculate_gp_t_from_fbsoL(
        fbsoL = fbsoL,
        method = method
    )
    
    if (iSample %in% highCovInLow) {
        save(gp_t, file = file_dosages(tempdir, iSample, regionName, "piece.gp_t"))
    }
    
    ## 
    eij <- round(gp_t[2, ] + 2 * gp_t[3, ], 3) ## prevent weird rounding issues
    fij <- round(gp_t[2, ] + 4 * gp_t[3, ], 3) ##
    
    infoCount[, 1] <- infoCount[, 1, drop = FALSE] + eij
    infoCount[, 2] <- infoCount[, 2, drop = FALSE] + (fij - eij**2)
    ## this returns un-transposed results
    max_gen <- get_max_gen_rapid(gp_t)
    ## hweCount is NOT transposed!
    hweCount[max_gen] <- hweCount[max_gen] + 1 ## hmmmmm not ideal
    afCount <- afCount + (eij) / 2
    
    ##
    if (output_format == "bgvcf") {
        if (output_haplotype_dosages) {
            if (method == "pseudoHaploid") {
                q_t <- fbsoL[[1]][["gammaEK_t"]] + fbsoL[[2]][["gammaEK_t"]]
            } else {
                ## do this here I suppose? 
                q_t <- 2 * fbsoL[[1]][["gammaEK_t"]]
            }
        } else {
            q_t <- matrix()
        }
        vcf_matrix_to_out[, iiSample] <- rcpp_make_column_of_vcf(
            gp_t = gp_t,
            use_read_proportions = FALSE,
            use_state_probabilities = output_haplotype_dosages,
            add_x_2_cols = FALSE,
            read_proportions = matrix(),
            q_t = q_t,
            x_t = matrix()
        )
        if (do_phasing) {
    
            ##
            w <- first_snp_in_region:last_snp_in_region
            hap1 <- rcpp_int_expand(phasing[, 1], nSNPs)[w]
            hap2 <- rcpp_int_expand(phasing[, 2], nSNPs)[w]
    
            vcf_matrix_to_out[, iiSample] <-
                paste0(
                    hap1, "|", hap2, 
                    substring(vcf_matrix_to_out[, iiSample], first = 4, last = 100L)
                )
    
        }
        
    } else if (output_format == "bgen") {
        rrbgen::rcpp_place_gp_t_into_output(
            gp_t,
            gp_raw_t, ## storage matrix
            iiSample, ## relative position
            nSNPs = ncol(gp_t), ## ncol(gp_t) = nSNPsInRegion here
            B_bit_prob
        )
    }
    }
    
