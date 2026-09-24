# Differential accessibility: port of 08_atac_differential/01_atac_diffbind.Rmd.
#   counted DBA -> dba.normalize() (DiffBind default: library size, DESeq2)
#   -> contrasts by Condition masks (config/contrasts.tsv)
#   -> dba.analyze (DESeq2) -> dba.report(th = 1) per contrast
#   -> full + significant tables, summary, MA plots, PCA, correlation heatmap,
#      analysed DBA saved as RDS.

if MODULES.get("diff", True) and CONTRASTS:

    rule diffbind_analyze:
        input:
            dba="results/counts/dba_counted.rds",
            contrasts=config["contrasts"],
        output:
            rds="results/diff/depth/dba_analyzed.rds",
            summary="results/diff/depth/summary.tsv",
            tables=expand("results/diff/depth/tables/{label}_all.tsv", label=CONTRASTS),
            sig=expand("results/diff/depth/tables/{label}_sig.tsv", label=CONTRASTS),
            ma="results/diff/depth/plots/ma_plots.pdf",
            pca="results/diff/depth/plots/pca.pdf",
            heatmap="results/diff/depth/plots/correlation_heatmap.pdf",
            barplot="results/diff/depth/plots/summary_barplot.pdf",
        log:
            "logs/diff/diffbind_analyze.log",
        conda:
            "../envs/r.yaml"
        threads: threads("diffbind", 8)
        resources:
            mem_mb=32000,
            runtime=480,
        params:
            helpers=workflow.source_path("../scripts/common.R"),
            fdr=DIFF["fdr"],
            lfc=DIFF["lfc"],
            tabdir=lambda wildcards, output: os.path.dirname(output.tables[0]),
        script:
            "../scripts/diffbind_analyze.R"
