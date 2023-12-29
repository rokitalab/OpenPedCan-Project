#!/bin/bash
# C. Bethell and C. Savonen for CCDL 2019, J. Rokita for D3b 2023
# Run focal-cn-file-preparation module
#
# Usage: bash run-prepare-cn.sh

set -e
set -o pipefail

# Run original files - will not by default
RUN_ORIGINAL=${RUN_ORIGINAL:-0}

# Run testing files for circle CI - will not by default
IS_CI=${OPENPBTA_TESTING:-0}

RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

scratch_dir=../../scratch
data_dir=../../data
results_dir=../../analyses/focal-cn-file-preparation/results
gtf_file=${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz

if [[ "$RUN_FOR_SUBTYPING" -eq "1" ]]
then
  histologies_file=${data_dir}/histologies-base.tsv
else
  histologies_file=${data_dir}/histologies.tsv
fi

# Prep the consensus SEG file data
Rscript --vanilla -e "rmarkdown::render('02-add-ploidy-consensus.Rmd', clean = TRUE)"

# Run annotation step for consensus file
Rscript --vanilla 04-prepare-cn-file.R \
--cnv_file ${results_dir}/consensus_seg_with_status.tsv \
--gtf_file $gtf_file \
--metadata $histologies_file \
--filename_lead "consensus_seg_annotated_cn" \
--seg


# if we want to process the CNV data from the original callers
# (e.g., CNVkit, ControlFreeC)
if [ "$RUN_ORIGINAL" -gt "0" ]
then

# Prep the CNVkit data
Rscript --vanilla -e "rmarkdown::render('01-add-ploidy-cnvkit.Rmd', clean = TRUE)"

# Run annotation step for CNVkit WXS only
Rscript --vanilla 04-prepare-cn-file.R \
--cnv_file ${results_dir}/cnvkit_with_status.tsv \
--gtf_file $gtf_file \
--metadata $histologies_file \
--filename_lead "cnvkit_annotated_cn" \
--seg \
--runWXSonly

# Run annotation step for ControlFreeC tumor only
Rscript --vanilla 04-prepare-cn-file.R \
--cnv_file ${data_dir}/cnv-controlfreec-tumor-only.tsv.gz \
--gtf_file $gtf_file \
--metadata $histologies_file \
--filename_lead "freec-tumor-only_annotated_cn" \
--controlfreec

# Run merging for all annotated files 
Rscript --vanilla 07-consensus-annotated-merge.R \
--cnvkit_auto ${results_dir}/cnvkit_annotated_cn_wxs_autosomes.tsv.gz \
--cnvkit_x_and_y ${results_dir}/cnvkit_annotated_cn_wxs_x_and_y.tsv.gz \
--consensus_auto ${results_dir}/consensus_seg_annotated_cn_autosomes.tsv.gz \
--consensus_x_and_y ${results_dir}/consensus_seg_annotated_cn_x_and_y.tsv.gz \
--cnv_tumor_auto ${results_dir}/freec-tumor-only_annotated_cn_autosomes.tsv.gz \
--cnv_tumor_x_and_y ${results_dir}/freec-tumor-only_annotated_cn_x_and_y.tsv.gz \
--outdir ${results_dir}

fi
