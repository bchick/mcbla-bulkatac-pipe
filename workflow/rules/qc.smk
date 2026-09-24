# QC: port of 2.0_atac_qc.sh (deepTools) plus MultiQC.
#   RPGC bigWigs (bin 10), TSS profile +/-2 kb, bamPEFragmentSize,
#   multiBamSummary 500 bp bins -> Spearman heatmap + PCA, fingerprint.

DT = config["qc"]


rule bigwig:
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
    output:
        "results/bigwig/{lib}.bw",
    log:
        "logs/qc/{lib}.bamCoverage.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        bin=DT["bigwig_bin_size"],
        gsize=REF["effective_genome_size"],
        extra=DT.get("bamcoverage_extra", ""),
    shell:
        """
        bamCoverage --bam {input.bam} --outFileName {output} --outFileFormat bigwig \
            --binSize {params.bin} --normalizeUsing RPGC \
            --effectiveGenomeSize {params.gsize} {params.extra} \
            --numberOfProcessors {threads} > {log} 2>&1
        """


rule tss_matrix:
    input:
        bw=expand("results/bigwig/{lib}.bw", lib=LIBS),
        gtf=REF["gtf"],
    output:
        matrix="results/qc/deeptools/tss_matrix.gz",
    log:
        "logs/qc/computeMatrix.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=16000,
        runtime=240,
    params:
        flank=DT["tss_flank"],
    shell:
        """
        computeMatrix reference-point -S {input.bw} -R {input.gtf} \
            --referencePoint TSS -a {params.flank} -b {params.flank} --skipZeros \
            -p {threads} -o {output.matrix} > {log} 2>&1
        """


rule tss_profile:
    input:
        "results/qc/deeptools/tss_matrix.gz",
    output:
        png="results/qc/deeptools/tss_enrichment_profile.png",
        data="results/qc/deeptools/tss_enrichment_profile.tab",
    log:
        "logs/qc/plotProfile.log",
    conda:
        "../envs/deeptools.yaml"
    params:
        labels=LIBS,
    shell:
        """
        plotProfile -m {input} --plotTitle "TSS Enrichment" \
            --samplesLabel {params.labels} --outFileNameData {output.data} \
            -o {output.png} > {log} 2>&1
        """


rule fragment_size_plot:
    input:
        bam=all_lib_bams(),
        bai=[f"{b}.bai" for b in all_lib_bams()],
    output:
        png="results/qc/deeptools/fragment_size_distribution.png",
        table="results/qc/deeptools/fragment_size_table.tsv",
        raw="results/qc/deeptools/fragment_size_raw.tsv",
    log:
        "logs/qc/bamPEFragmentSize.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        labels=LIBS,
    shell:
        """
        bamPEFragmentSize --bamfiles {input.bam} --samplesLabel {params.labels} \
            --numberOfProcessors {threads} \
            --plotTitle "Fragment Size Distribution (ATAC-seq)" \
            --table {output.table} --outRawFragmentLengths {output.raw} \
            -o {output.png} > {log} 2>&1
        """


rule multibamsummary:
    input:
        bam=all_lib_bams(),
        bai=[f"{b}.bai" for b in all_lib_bams()],
    output:
        "results/qc/deeptools/multiBamSummary.npz",
    log:
        "logs/qc/multiBamSummary.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=16000,
        runtime=480,
    params:
        labels=LIBS,
        bin=DT["summary_bin_size"],
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam} --labels {params.labels} \
            --binSize {params.bin} --numberOfProcessors {threads} \
            --outFileName {output} > {log} 2>&1
        """


rule plot_correlation:
    input:
        "results/qc/deeptools/multiBamSummary.npz",
    output:
        png="results/qc/deeptools/correlation_spearman.png",
        tsv="results/qc/deeptools/correlation_spearman.tsv",
    log:
        "logs/qc/plotCorrelation.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotCorrelation --corData {input} --corMethod spearman --whatToPlot heatmap \
            --plotNumbers --skipZeros --plotTitle "Spearman Correlation (ATAC-seq)" \
            --outFileCorMatrix {output.tsv} --plotFile {output.png} > {log} 2>&1
        """


rule plot_pca:
    input:
        "results/qc/deeptools/multiBamSummary.npz",
    output:
        png="results/qc/deeptools/pca_plot.png",
        tsv="results/qc/deeptools/pca_data.tsv",
    log:
        "logs/qc/plotPCA.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotPCA --corData {input} --plotTitle "PCA of ATAC-seq Samples" \
            --outFileNameData {output.tsv} --plotFile {output.png} > {log} 2>&1
        """


rule fingerprint:
    input:
        bam=all_lib_bams(),
        bai=[f"{b}.bai" for b in all_lib_bams()],
    output:
        png="results/qc/deeptools/fingerprint.png",
        metrics="results/qc/deeptools/fingerprint_metrics.tsv",
        counts="results/qc/deeptools/fingerprint_counts.tsv",
    log:
        "logs/qc/plotFingerprint.log",
    conda:
        "../envs/deeptools.yaml"
    threads: threads("deeptools", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        labels=LIBS,
    shell:
        """
        plotFingerprint --bamfiles {input.bam} --labels {params.labels} \
            --numberOfProcessors {threads} --plotTitle "Fingerprint (ATAC-seq)" \
            --outQualityMetrics {output.metrics} --outRawCounts {output.counts} \
            --plotFile {output.png} > {log} 2>&1
        """


rule multiqc:
    input:
        multiqc_inputs,
    output:
        "results/qc/multiqc/multiqc_report.html",
    log:
        "logs/qc/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        config=workflow.source_path("../resources/multiqc_config.yaml"),
    shell:
        """
        (
        stage=$(mktemp -d)
        for f in {input}; do
            base=$(basename "$f")
            case "$f" in
              *alignment_qc_report.tsv) dest=alignment_qc_mqc.tsv ;;
              *individual/peak_summary.tsv) dest=peaks_individual_mqc.tsv ;;
              *merged_stringent/peak_summary.tsv) dest=peaks_merged_stringent_mqc.tsv ;;
              *idr_summary.tsv) dest=idr_summary_mqc.tsv ;;
              *) dest="$base" ;;
            esac
            ln -sf "$(realpath "$f")" "$stage/$dest"
        done
        multiqc --force --config {params.config} -o {params.outdir} "$stage"
        rm -rf "$stage"
        ) > {log} 2>&1
        """
