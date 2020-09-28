workflow STAAR_genomewide {

    # run_null_model inputs
    File? null_file_precompute
    File? pheno_file
    String? null_file_name
    String? sample_id
    String? outcome
    String? outcome_type = "continuous"
    String? covariates = "NA"
    File? kinship_file
    String? group_id = "NA"
    Int? null_memory = 25
    Int? null_disk = 50

    # run_genomewide inputs
    Array[File] geno_files
    Array[File]? annot_files
    File? agg_file
    String results_file
    String? agds_file = "None"
    String? agds_annot_channels = "None"
    String? maf_thres = "0.05"
    String? mac_thres = "1"
    Int? window_length = 2000
    Int? step_length = 1000
    # Compute inputs
    Int? num_cores = 3
    Int? num_chunk_div = 3
    Int? test_memory = 25
    Int? test_disk = 50

    if (!defined(null_file_precompute)) {
        call run_null_model {
            input:
                pheno_file = pheno_file,
                null_file_name = null_file_name,
                sample_id = sample_id,
                outcome = outcome,
                outcome_type = outcome_type,
                covariates = covariates,
                kinship_file = kinship_file,
                group_id = group_id,
                null_memory = null_memory,
                null_disk = null_disk
        }
    }

    File null_file = select_first([null_file_precompute, run_null_model.null_model])

    Array[File] annot_avail = select_first([ annot_files, []])

    if (length(annot_avail) == length(geno_files)) {
        Array[Pair[File,File]] geno_annot_pairs = zip(geno_files, annot_avail)
        scatter (geno_annot_set in geno_annot_pairs) {
        File geno_in = geno_annot_set.left
        File annot_in = geno_annot_set.right
            call run_genomewide {
                input:
                    null_file = null_file,
                    geno_file = geno_in,
                    annot_file = annot_in,
                    results_file = results_file,
                    agg_file = agg_file,
                    maf_thres = maf_thres,
                    mac_thres = mac_thres,
                    window_length = window_length,
                    step_length = step_length,
                    num_cores = num_cores,
                    num_chunk_div = num_chunk_div,
                    test_memory = test_memory,
                    test_disk = test_disk
            }
        }
    }

    if (!defined(annot_files)) {
        scatter (geno_in in geno_files) {
            call run_genomewide as run_genomewide_annotfree {
                input:
                    null_file = null_file,
                    geno_file = geno_in,
                    results_file = results_file,
                    agds_file = agds_file,
                    agds_annot_channels = agds_annot_channels,
                    agg_file = agg_file,
                    maf_thres = maf_thres,
                    mac_thres = mac_thres,
                    window_length = window_length,
                    step_length = step_length,
                    num_cores = num_cores,
                    num_chunk_div = num_chunk_div,
                    test_memory = test_memory,
                    test_disk = test_disk
            }
        }
    }

    output {
        File null_model = null_file
        Array[File]? result_genomewide = run_genomewide.results
        Array[File]? result_genomewide_annotfree = run_genomewide_annotfree.results
    }
}



task run_null_model {
    File pheno_file
    String null_file_name
    String sample_id
    String outcome
    String outcome_type
    String covariates
    File kinship_file
    String group_id
    Int null_memory
    Int null_disk

    command {
        Rscript /STAAR_null_model.R ${pheno_file} ${null_file_name} ${sample_id} ${outcome} ${outcome_type} ${covariates} ${kinship_file} ${group_id}
    }

    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "${null_memory} GB"
        disks: "local-disk ${null_disk} HDD"
    }

    output {
        File null_model = select_first(glob("${null_file_name}*"))
    }
}

task run_genomewide {
    File null_file
    File geno_file
    File? annot_file
    String results_file
    String? agds_file
    String? agds_annot_channels
    File? agg_file
    String maf_thres
    String mac_thres
    Int window_length
    Int step_length
    Int num_cores
    Int num_chunk_div
    Int test_memory
    Int test_disk

    command {
        Rscript /STAAR_genomewide.R ${null_file} ${geno_file} ${default="None" annot_file} ${results_file} ${default="None" agds_file} ${default="None" agds_annot_channels} ${default="None" agg_file} ${maf_thres} ${mac_thres} ${window_length} ${step_length} ${num_cores} ${num_chunk_div}
    }
    runtime {
        docker: "quay.io/sheilagaynor/staar_slim"
        memory: "${test_memory} GB"
        disks: "local-disk ${test_disk} HDD"
    }
    output {
        File results = select_first(glob("*.gz"))
    }
}
