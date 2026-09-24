# Normalization sensitivity check: port of 08b_atac_normalization_comparison
# and the decision rule of 22_normalization_decision.
#   csaw     : dba.normalize(method = DBA_DESEQ2, normalize = DBA_NORM_NATIVE,
#              background = TRUE)  (csaw-style 15 kb background bins)
#   quantile : log2(count + 0.5) -> preprocessCore quantile normalization
#              -> limma lmFit / contrasts.fit / eBayes          (optional)
# The same contrasts as the diff module are re-run on the same counted DBA.
# A contrast whose gained or lost count moves by more than
# normcheck.sensitivity_threshold (default 20%) versus the default
# normalization is flagged "normalization-sensitive".

NC = config["normcheck"]
NC_METHODS = [m for m in NC.get("methods", ["csaw"]) if m in ("csaw", "quantile")]

if RUN_NORMCHECK:

    rule normcheck_csaw:
        input:
            dba="results/counts/dba_counted.rds",
            contrasts=config["contrasts"],
        output:
            rds="results/normcheck/csaw/dba_csaw.rds",
            summary="results/normcheck/csaw/summary.tsv",
            sizefactors="results/normcheck/csaw/size_factors.tsv",
            tables=expand(
                "results/normcheck/csaw/tables/{label}_all.tsv", label=CONTRASTS
            ),
        log:
            "logs/normcheck/csaw.log",
        conda:
            "../envs/r.yaml"
        threads: threads("diffbind", 8)
        resources:
            mem_mb=64000,
            runtime=720,
        params:
            helpers=workflow.source_path("../scripts/common.R"),
            fdr=DIFF["fdr"],
            lfc=DIFF["lfc"],
            tabdir=lambda wildcards, output: os.path.dirname(output.tables[0]),
        script:
            "../scripts/normcheck_csaw.R"

    rule normcheck_quantile:
        input:
            dba="results/counts/dba_counted.rds",
            contrasts=config["contrasts"],
        output:
            summary="results/normcheck/quantile/summary.tsv",
            diagnostics="results/normcheck/quantile/quantile_diagnostics.tsv",
            tables=expand(
                "results/normcheck/quantile/tables/{label}_all.tsv", label=CONTRASTS
            ),
        log:
            "logs/normcheck/quantile.log",
        conda:
            "../envs/r.yaml"
        resources:
            mem_mb=32000,
            runtime=240,
        params:
            helpers=workflow.source_path("../scripts/common.R"),
            fdr=DIFF["fdr"],
            lfc=DIFF["lfc"],
            tabdir=lambda wildcards, output: os.path.dirname(output.tables[0]),
        script:
            "../scripts/normcheck_quantile.R"

    rule normcheck_compare:
        input:
            depth_summary="results/diff/depth/summary.tsv",
            depth_tables=expand(
                "results/diff/depth/tables/{label}_all.tsv", label=CONTRASTS
            ),
            summaries=expand("results/normcheck/{m}/summary.tsv", m=NC_METHODS),
            tables=expand(
                "results/normcheck/{m}/tables/{label}_all.tsv",
                m=NC_METHODS,
                label=CONTRASTS,
            ),
        output:
            table="results/normcheck/norm_comparison.tsv",
            verdict="results/normcheck/norm_verdict.tsv",
            barplot="results/normcheck/norm_comparison_barplot.pdf",
        log:
            "logs/normcheck/compare.log",
        conda:
            "../envs/r.yaml"
        params:
            helpers=workflow.source_path("../scripts/common.R"),
            methods=NC_METHODS,
            labels=CONTRASTS,
            fdr=DIFF["fdr"],
            threshold=NC["sensitivity_threshold"],
            min_abs=NC["min_abs_change"],
        script:
            "../scripts/normcheck_compare.R"
