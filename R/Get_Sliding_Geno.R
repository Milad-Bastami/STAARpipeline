favor_names <- c("VarInfo","vid","variant_vcf","variant_annovar","chromosome","start_position","end_position","ref_annovar","alt_annovar","position","ref_vcf","alt_vcf","aloft_value","aloft_description","apc_conservation","apc_conservation_v2","apc_epigenetics_active","apc_epigenetics","apc_epigenetics_repressed","apc_epigenetics_transcription","apc_local_nucleotide_diversity","apc_local_nucleotide_diversity_v2","apc_local_nucleotide_diversity_v3","apc_mappability","apc_micro_rna","apc_mutation_density","apc_protein_function","apc_protein_function_v2","apc_protein_function_v3","apc_proximity_to_coding","apc_proximity_to_coding_v2","apc_proximity_to_tsstes","apc_transcription_factor","bravo_an","bravo_af","filter_status","cage_enhancer","cage_promoter","cage_tc","clnsig","clnsigincl","clndn","clndnincl","clnrevstat","origin","clndisdb","clndisdbincl","geneinfo","polyphen2_hdiv_score","polyphen2_hvar_score","mutation_taster_score","mutation_assessor_score","metasvm_pred","fathmm_xf","funseq_value","funseq_description","genecode_comprehensive_category","genecode_comprehensive_info","genecode_comprehensive_exonic_category","genecode_comprehensive_exonic_info","genehancer","af_total","af_asj_female","af_eas_female","af_afr_male","af_female","af_fin_male","af_oth_female","af_ami","af_oth","af_male","af_ami_female","af_afr","af_eas_male","af_sas","af_nfe_female","af_asj_male","af_raw","af_oth_male","af_nfe_male","af_asj","af_amr_male","af_amr_female","af_sas_female","af_fin","af_afr_female","af_sas_male","af_amr","af_nfe","af_eas","af_ami_male","af_fin_female","linsight","gc","cpg","min_dist_tss","min_dist_tse","sift_cat","sift_val","polyphen_cat","polyphen_val","priphcons","mamphcons","verphcons","priphylop","mamphylop","verphylop","bstatistic","chmm_e1","chmm_e2","chmm_e3","chmm_e4","chmm_e5","chmm_e6","chmm_e7","chmm_e8","chmm_e9","chmm_e10","chmm_e11","chmm_e12","chmm_e13","chmm_e14","chmm_e15","chmm_e16","chmm_e17","chmm_e18","chmm_e19","chmm_e20","chmm_e21","chmm_e22","chmm_e23","chmm_e24","chmm_e25","gerp_rs","gerp_rs_pval","gerp_n","gerp_s","encodeh3k4me1_sum","encodeh3k4me2_sum","encodeh3k4me3_sum","encodeh3k9ac_sum","encodeh3k9me3_sum","encodeh3k27ac_sum","encodeh3k27me3_sum","encodeh3k36me3_sum","encodeh3k79me2_sum","encodeh4k20me1_sum","encodeh2afz_sum","encode_dnase_sum","encodetotal_rna_sum","grantham","freq100bp","rare100bp","sngl100bp","freq1000bp","rare1000bp","sngl1000bp","freq10000bp","rare10000bp","sngl10000bp","remap_overlap_tf","remap_overlap_cl","cadd_rawscore","cadd_phred","k24_bismap","k24_umap","k36_bismap","k36_umap","k50_bismap","k50_umap","k100_bismap","k100_umap","nucdiv","rdhs","recombination_rate","refseq_category","refseq_info","refseq_exonic_category","refseq_exonic_info","super_enhancer","tg_afr","tg_all","tg_amr","tg_eas","tg_eur","tg_sas","ucsc_category","ucsc_info","ucsc_exonic_category","ucsc_exonic_info")
favor_colnames <- paste0('annotation/info/FunctionalAnnotation/FunctionalAnnotation/', favor_names)
# based on Sliding_Window_Multiple
# save variants and genotypes for sliding window results
# start_loc and end_loc are coordinates of the final
# sliding window retutned by staarpipeline

Get_Sliding_Geno <-
  function(chr,
           start_loc,
           end_loc,
           #sliding_window_length = 2000,
           genofile,
           obj_nullmodel,
           rare_maf_cutoff = 0.01,
           rv_num_cutoff = 2,
           QC_label = "annotation/filter",
           variant_type = c("SNV", "Indel", "variant"),
           geno_missing_imputation = c("mean", "minor"),
           Annotation_dir = "annotation/info/FunctionalAnnotation",
           Annotation_name_catalog,
           Use_annotation_weights = c(TRUE, FALSE),
           Annotation_name = NULL,
           SPA_p_filter = FALSE,
           p_filter_cutoff = 0.05,
           silent = FALSE,
           vaf = 0.35,
           dp = 10,
           gq = 80,
           save_df=TRUE,
           out_dir,
           case_idx=1:1644,
           control_idx = 1645:4849,
           doSKAT=FALSE,
           skatNULL=NULL) {

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- as.character(obj_nullmodel$id_include)
  n_pheno <- obj_nullmodel$n.pheno

  ## get SNV id
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant")
  {
    SNVlist <- filter == "PASS"
  }

  if(variant_type=="SNV")
  {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }

  if(variant_type=="Indel")
  {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  variant.id <- seqGetData(genofile, "variant.id")

  ## Position
  position <- as.numeric(seqGetData(genofile, "position"))
  position_SNV <- position[SNVlist]

  # first filtering
  is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
  seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

  ## genotype id
  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index

  if(sum(is.in)>=2)
  {
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage_alt")
    Geno <- Geno[id.genotype.match,,drop=FALSE]

    ## impute missing
    if(!is.null(dim(Geno)))
    {
      if(dim(Geno)[2]>0)
      {
        flip = flip_alleles(Geno, snps = 'cols')
        Geno.flip <- flip$dosmat
        flipped_loci <- flip$flipped_loci
        ## filtering for VAF,GQ,DP
        VAF <- seqGetData(genofile, 'annotation/format/VAF')
        VAF <- VAF[id.genotype.match, , drop = FALSE]
        # for flipped loci, ref is minor allele
        # need to correct VAF accordingly
        if(!is.null(flipped_loci)){
          VAF[, flipped_loci] <- 1 - VAF[, flipped_loci]
        }
        GQ <- seqGetData(genofile, 'annotation/format/GQ')
        GQ <- GQ[id.genotype.match, , drop = FALSE]
        DP <- seqGetData(genofile, 'annotation/format/DP')
        DP <- DP[id.genotype.match, , drop = FALSE]
        QC <- (VAF >= vaf) & GQ >= gq & DP >= dp
        isgt0 <- Geno.flip > 0
        # non-reference GTs (i.e. 1,2) are set to zero
        # if they fail defined QC; NA: because it wont
        # affect MAF calculation by staar.
        Geno.adj <-
          ifelse(QC == FALSE & isgt0 == TRUE, 0 , Geno.flip)
        Geno.adj <- ifelse(is.na(Geno.adj), 0, Geno.adj)

        if (geno_missing_imputation == "mean") {
          imputed <- matrix_flip_mean(Geno.adj)
          Geno <- imputed$Geno
          MAF <- imputed$MAF
        }
        # second filtering:
        # filtering variants with maf above threshold
        is.in.rare <- (MAF<rare_maf_cutoff)&(MAF>0)
        Geno <- Geno[, is.in.rare]
        Geno.adj <- Geno.adj[, is.in.rare]
        MAF <- MAF[is.in.rare]
        variant.id.rare <- seqGetData(genofile, 'variant.id')
        variant.id.rare <- variant.id.rare[is.in.rare]
        seqSetFilter(genofile, variant.id = variant.id.rare)

        ## perform SKAT test
        if(doSKAT) {
          skat = SKAT::SKATBinary(as.matrix(Geno.adj),
                                  get(skatNULL),
                                  method = 'Burden')
          skat_res <- cbind.data.frame(
            window = paste0("chr", chr, "_", start_loc, "_", end_loc),
            pval_SKATbin_burden = skat$p.value,
            SKATbin_burden_n.test = skat$param$n.marker.test,
            SKATbin_burden_n.marker = skat$param$n.marker,
            MAC = skat$MAC
          )
          fname = paste0(out_dir, "/chr",chr,"_",start_loc,"_",end_loc,"_skat.txt")
          write.table(
            skat_res,
            file = fname,
            sep = '\t',
            row.names = F,
            col.names = T,
            quote = F
          )
          rm(list = c('fname', 'skat_res', 'skat'))
          gc()

        }
      }
    }

    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL

    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }

            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }

        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }

    Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[is.in.rare, ]
    ################################
    #
    ################################
    if (save_df) {
      ## df_var & df_sam
      nrow_df = sum(is.in.rare)
      df_var = data.frame(
        id = as.character(seqGetData(genofile, 'annotation/id')),
        chr = as.numeric(seqGetData(genofile, "chromosome")),
        pos = as.numeric(seqGetData(genofile, "position")),
        ref = as.character(seqGetData(genofile, "$ref")),
        alt = as.character(seqGetData(genofile, "$alt")),
        acM = apply(Geno.adj, 2, function(x)
          (2 * sum(x == 0L, na.rm = T)) + sum(x == 1L, na.rm = T)),
        acm = apply(Geno.adj, 2, function(x)
          sum(x, na.rm = T)),
        maf = apply(Geno.adj, 2, function(x)
          (sum(x, na.rm = T) / 2) / sum(!is.na(x))),
        nMM = apply(Geno.adj, 2, function(x)
          sum(x == 0L, na.rm = T)),
        nMm = apply(Geno.adj, 2, function(x)
          sum(x == 1L, na.rm = T)),
        nmm = apply(Geno.adj, 2, function(x)
          sum(x == 2L, na.rm = T)),
        nMissing = apply(Geno.adj, 2, function(x)
          sum(is.na(x))),
        fMissing_info = as.numeric(seqGetData(genofile, 'annotation/info/F_MISSING')),
        type = as.character(isSNV(genofile)),
        an_info = as.numeric(seqGetData(genofile, 'annotation/info/AN')),
        ac_info = as.numeric(seqGetData(genofile, 'annotation/info/AC')),
        ac_hom = as.numeric(seqGetData(genofile, 'annotation/info/AC_Hom')),
        ac_het = as.numeric(seqGetData(genofile, 'annotation/info/AC_Het')),
        maf_info = as.numeric(seqGetData(genofile, 'annotation/info/MAF')),
        maf_staar = MAF
      )
      # minor allele
      df_var$minor_allele = 'alt'
      if (!is.null(flipped_loci))
      {
        df_var[flipped_loci, 'minor_allele'] = 'ref'
      }
      rm(list ='flipped_loci')
      gc()
      # ISKS variant freq
      Geno.adj.isks <- Geno.adj[case_idx, , drop = FALSE]
      df_var$isks_acM = apply(Geno.adj.isks, 2, function(x)
        (2 * sum(x == 0L, na.rm = T)) + sum(x == 1L, na.rm = T))
      df_var$isks_acm = apply(Geno.adj.isks, 2, function(x)
        sum(x, na.rm = T))
      df_var$isks_maf = df_var$isks_acm / df_var$isks_acM
      df_var$isks_nMissing = apply(Geno.adj.isks, 2, function(x)
        sum(is.na(x)))
      df_var$isks_nMM = apply(Geno.adj.isks, 2, function(x)
        sum(x == 0L, na.rm = T))
      df_var$isks_nMm = apply(Geno.adj.isks, 2, function(x)
        sum(x == 1L, na.rm = T))
      df_var$isks_nmm = apply(Geno.adj.isks, 2, function(x)
        sum(x == 2L, na.rm = T))

      ## MGRB variant freq
      Geno.adj.mgrb <- Geno.adj[control_idx, , drop = FALSE]
      df_var$mgrb_acM = apply(Geno.adj.mgrb, 2, function(x)
        (2 * sum(x == 0L, na.rm = T)) + sum(x == 1L, na.rm = T))
      df_var$mgrb_acm = apply(Geno.adj.mgrb, 2, function(x)
        sum(x, na.rm = T))
      df_var$mgrb_maf = df_var$mgrb_acm / df_var$mgrb_acM
      df_var$mgrb_nMissing = apply(Geno.adj.mgrb, 2, function(x)
        sum(is.na(x)))
      df_var$mgrb_nMM = apply(Geno.adj.mgrb, 2, function(x)
        sum(x == 0L, na.rm = T))
      df_var$mgrb_nMm = apply(Geno.adj.mgrb, 2, function(x)
        sum(x == 1L, na.rm = T))
      df_var$mgrb_nmm = apply(Geno.adj.mgrb, 2, function(x)
        sum(x == 2L, na.rm = T))
      ## adding favor annotations to varDF
      favor_anno <- seqGetData(genofile, var.name = favor_colnames)
      favor_anno <- as.data.frame(favor_anno)
      colnames(favor_anno) <- favor_names
      gnocchi <- seqGetData(genofile, 'annotation/info/Annovar/gnocchi')
      ensemblRB <- seqGetData(genofile, 'annotation/info/Annovar/ensemblRB')
      cCRE <- seqGetData(genofile, 'annotation/info/Annovar/cCRE')
      gnomad_MAF_grpmax <- seqGetData(genofile, 'annotation/info/Annovar/gnomad40_genome_MAF_grpmax')
      gnomad_MAF_nfe <- seqGetData(genofile, 'annotation/info/Annovar/gnomad40_genome_MAF_nfe')
      favor_anno[, 'gnocchi'] <- gnocchi
      favor_anno[, 'ensemblRB'] <- ensemblRB
      favor_anno[, 'cCRE'] <- cCRE
      favor_anno[, 'gnomad_MAF_grpmax'] <- gnomad_MAF_grpmax
      favor_anno[, 'gnomad_MAF_nfe'] <- gnomad_MAF_nfe
      df_var <- cbind.data.frame(df_var, favor_anno)
      ## sample table
      df_sam <- data.frame(
        sample.id = seqGetData(genofile, 'sample.id'),
        hit_AD = apply(Geno.adj, 1, function(x)
          sum(x, na.rm = T) >= 1),
        hit_AR = apply(Geno.adj, 1, function(x)
          sum(x, na.rm = T) >= 2),
        total_RVs = rep(ncol(Geno.adj), nrow(Geno.adj)),
        nMM = apply(Geno.adj, 1, function(x)
          sum(x == 0L, na.rm = T)),
        nMm = apply(Geno.adj, 1, function(x)
          sum(x == 1L, na.rm = T)),
        nmm = apply(Geno.adj, 1, function(x)
          sum(x == 2L, na.rm = T)),
        nMissing = apply(Geno.adj, 1, function(x)
          sum(is.na(x)))
      )
      df_sam$fMissing = df_sam$nMissing / df_sam$total_RVs
      ## appennding genotype gene-level
      df_geno <- t(Geno.adj)
      colnames(df_geno) <- phenotype.id
      df_geno <- cbind(varID = df_var$id, df_geno)
      ###### savings ####
      ## add PHRED weights
      if (!is.null(Anno.Int.PHRED.sub)) {
        if (nrow(Anno.Int.PHRED.sub) != ncol(Geno.adj)) {
          message('different number of rows (Anno & Geno)')
        } else {
          df_geno <- cbind.data.frame(df_geno, Anno.Int.PHRED.sub)
          fname <- paste0(out_dir,'/chr',chr,'_',start_loc,'_', end_loc, '_Geno_wt.txt')
          write.table(
            df_geno,
            file = fname,
            row.names = F,
            quote = F,
            sep = '\t',
            col.names = TRUE)
          rm(fname)
          gc()

          df_var <- cbind.data.frame(df_var, Anno.Int.PHRED.sub)
          fname <- paste0(out_dir,'/chr',chr,'_',start_loc,'_', end_loc, '_vDF_wt.txt')
          write.table(
            df_var,
            file = fname,
            row.names = F,
            quote = F,
            sep = '\t',
            col.names = TRUE)
          rm(fname)
          gc()
        }

      } else {
        fname <- paste0(out_dir,'/chr',chr,'_',start_loc,'_', end_loc, '_Geno.txt')
        write.table(
          df_geno,
          file = fname,
          row.names = F,
          quote = F,
          sep = '\t',
          col.names = TRUE)
        rm(fname)
        gc()

        fname <- paste0(out_dir,'/chr',chr,'_',start_loc,'_', end_loc, '_vDF.txt')
        write.table(
          df_var,
          file = fname,
          row.names = F,
          quote = F,
          sep = '\t',
          col.names = TRUE)
        rm(fname)
        gc()
      }
      fname <- paste0(out_dir,'/chr',chr,'_',start_loc,'_', end_loc, '_sDF.txt')
      write.table(
        df_sam,
        file = fname,
        row.names = F,
        quote = F,
        sep = '\t',
        col.names = TRUE,
      )
      rm(fname)
      gc()

      rm(list = c(
        'df_sam',
        'df_var',
        'df_geno',
        'Geno.adj.mgrb',
        'Geno.adj.isks'
      ))
      gc()
    }
    seqResetFilter(genofile)
  }
  }

