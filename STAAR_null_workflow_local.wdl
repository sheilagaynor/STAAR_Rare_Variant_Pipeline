workflow STAAR_null {
    File pheno_file
    String null_file
    String sample_id
    String outcome
    String? outcome_type = "continuous"
    String? covariates = "NA"
    File? kinship_file
    String? group_id = "NA"

    call run_null_model {
            input:
                pheno_file = pheno_file,
                null_file = null_file,
                sample_id = sample_id,
                outcome = outcome,
                outcome_type = outcome_type,
                covariates = covariates,
                kinship_file = kinship_file,
                group_id = group_id
    }

    output {
        File null_model = run_null_model.null_model
    }
}

task run_null_model {
    File pheno_file
    String null_file
    String sample_id
    String outcome
    String outcome_type
    String covariates
    File kinship_file
    String group_id

    command {
        Rscript ~/Desktop/STAAR_Rare_Variant_Pipeline/STAAR_null_model.R ${pheno_file} ${null_file} ${sample_id} ${outcome} ${outcome_type} ${covariates} ${kinship_file} ${group_id}
    }

    output {
        File null_model = select_first(glob("${null_file}*"))
    }
}
