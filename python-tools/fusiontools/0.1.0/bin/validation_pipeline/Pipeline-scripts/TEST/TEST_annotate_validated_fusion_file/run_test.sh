#!/bin/bash
source setup.sh
pipeline_codedir=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts
#in_file=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/stjude_validation/process_validated_fusions_file/105-SJ-Valid-fusions_COMPARED_St-Jude-fusions-removed_duplicates-removed_brokensamples.valfile
#out_file=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/stjude_validation/process_validated_fusions_file/105-SJ-Valid-fusions_COMPARED_St-Jude-fusions-removed_duplicates-removed_brokensamples-bedfile_intersected.valfile
in_file=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/kumar_2016_validation/kumar-validated_fusions_grch37_infofile.tsv
out_file=/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/kumar_2016_validation/kumar-validated_fusions_grch37_infofile.bedfile_intersected.valfile
python $pipeline_codedir/annotate_validated_fusion_file.py $in_file stjude_5callers $out_file
