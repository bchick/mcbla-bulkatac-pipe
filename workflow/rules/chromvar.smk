# TF activity: port of 05_atac_temporal_clustering/07_chromvar.Rmd.
#   counts over chromvar.peaks (min fragments filter) -> addGCBias
#   -> JASPAR2020 CORE vertebrate PWMs matched on 200 bp peak-centred windows
#      (motifmatchr) -> getBackgroundPeaks(niterations = 200, w = 0.1)
#   -> computeDeviations -> deviation Z-scores + variability.
# set.seed(2025) before GC bias and background sampling, as in the source.
# The genome is either a BSgenome package (chromvar.bsgenome) or the FASTA.

CV = config["chromvar"]

if RUN_CHROMVAR:

    rule chromvar:
        input:
            counts=counts_matrix(CV_PEAKS),
            fasta="results/reference/genome.fa",
            fai="results/reference/genome.fa.fai",
        output:
            z="results/chromvar/deviation_zscores.tsv",
            dev="results/chromvar/deviations.tsv",
            var="results/chromvar/variability.tsv",
            heatmap="results/chromvar/top_variable_heatmap.pdf",
            rds="results/chromvar/chromvar_deviations.rds",
        log:
            "logs/chromvar/chromvar.log",
        conda:
            CV.get("conda_env") or "../envs/r.yaml"
        threads: threads("chromvar", 8)
        resources:
            mem_mb=32000,
            runtime=480,
        params:
            helpers=workflow.source_path("../scripts/common.R"),
            libs=STAT_LIBS,
            conditions=[LIBRARIES.loc[l, "condition"] for l in STAT_LIBS],
            bsgenome=CV.get("bsgenome") or "",
            collection=CV["jaspar_collection"],
            tax_group=CV["jaspar_tax_group"],
            peak_width=CV["motif_window"],
            min_fragments=CV["min_peak_fragments"],
            iterations=CV["background_iterations"],
            seed=CV["seed"],
            n_top=CV["n_top_motifs"],
        script:
            "../scripts/chromvar.R"
