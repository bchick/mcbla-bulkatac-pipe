# Temporal programs: port of 05_atac_temporal_clustering (00_preprocessing +
# 01_degpatterns_clustering). Runs only when timecourse.run is true (needs
# complete `treatment` and `time` columns and timecourse.peaks).
#
# Per treatment series (baseline_treatments, e.g. unstimulated time 0, are
# shared by every series):
#   counts over timecourse.peaks -> keep rowMeans >= 10 -> VST (blind) on all samples
#   -> DESeq2 LRT (~time vs ~1) on the series -> padj < 0.01
#   -> range of per-time VST means >= 0.5 -> DEGreport::degPatterns
#      (minc 50, cutoff 0.5, set.seed(42))

if RUN_TIMECOURSE:

    rule timecourse_series:
        input:
            counts=counts_matrix(TC_PEAKS),
        output:
            lrt="results/timecourse/{series}/lrt_results.tsv",
            clusters="results/timecourse/{series}/degpatterns_clusters.tsv",
            profiles="results/timecourse/{series}/degpatterns_profiles.tsv",
            plot="results/timecourse/{series}/degpatterns_clusters.pdf",
            rds="results/timecourse/{series}/timecourse.rds",
        log:
            "logs/timecourse/{series}.log",
        conda:
            "../envs/r.yaml"
        threads: threads("deseq2", 4)
        resources:
            mem_mb=32000,
            runtime=480,
        params:
            helpers=workflow.source_path("../scripts/common.R"),
            series=lambda wildcards: wildcards.series,
            libs=STAT_LIBS,
            treatments=[LIBRARIES.loc[l, "treatment"] for l in STAT_LIBS],
            times=[LIBRARIES.loc[l, "time"] for l in STAT_LIBS],
            conditions=[LIBRARIES.loc[l, "condition"] for l in STAT_LIBS],
            baseline=TC.get("baseline_treatments") or [],
            min_mean_counts=TC["min_mean_counts"],
            lrt_padj=TC["lrt_padj"],
            vst_range=TC["vst_range"],
            minc=TC["minc"],
            cutoff=TC["cutoff"],
            seed=TC["seed"],
        script:
            "../scripts/timecourse.R"
