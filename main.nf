nextflow.enable.dsl=2

params.input = null
params.outdir = 'results'
params.min_genes = 200
params.max_pct_mt = 20
params.n_top_genes = 2000
params.resolution = 0.5

workflow {
    if (!params.input) {
        error "Please provide --input pointing to a .h5ad file or 10x mtx directory"
    }

    Channel
        .fromPath(params.input, checkIfExists: true)
        .set { input_data }

    RUN_SCANPY(input_data)
}

process RUN_SCANPY {
    tag "single-cell-qc-cluster"
    publishDir params.outdir, mode: 'copy'

    input:
    path input_data

    output:
    path "*.h5ad"
    path "figures/*.png"
    path "tables/*.csv"
    path "reports/*.md"

    script:
    """
    mkdir -p figures tables reports

    python /app/bin/run_scanpy_pipeline.py \
        --input ${input_data} \
        --out_h5ad ad_processed.h5ad \
        --figures figures \
        --tables tables \
        --report reports/summary.md \
        --min_genes ${params.min_genes} \
        --max_pct_mt ${params.max_pct_mt} \
        --n_top_genes ${params.n_top_genes} \
        --resolution ${params.resolution}
    """
}
