scripts_path=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts
separate_callers_dir=$validation_pipeline_dir/fusion_pipeline_output_files/smc_rna_sim45-separate
renamed_dir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/TEST/TEST_rename_cff_file_genes/output
file=
python $scripts_path/rename_cff_file_genes.py $separate_callers_dir/$file > $renamed_dir/$file_renamed
