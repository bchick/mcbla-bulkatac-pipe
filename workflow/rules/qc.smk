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


# ---------------------------------------------------------------------------
# QC metrics with ENCODE ATAC-seq thresholds -> results/qc/qc_summary.tsv
# ---------------------------------------------------------------------------
if INPUT_MODE == "fastq":

    rule fastqc:
        """FastQC on the raw reads (runs of one library concatenated)."""
        input:
            unpack(trim_inputs),
        output:
            r1="results/qc/fastqc/{lib}_R1_fastqc.zip",
            r2="results/qc/fastqc/{lib}_R2_fastqc.zip",
        log:
            "logs/qc/{lib}.fastqc.log",
        conda:
            "../envs/fastqc.yaml"
        threads: 2
        resources:
            mem_mb=2000,
            runtime=240,
        params:
            outdir=lambda wildcards, output: os.path.dirname(output.r1),
        shell:
            """
            (
            set -euo pipefail
            tmp=$(mktemp -d)
            # name the reads after the library so MultiQC shows library names
            ln -s "$(realpath {input.r1})" "$tmp/{wildcards.lib}_R1.fastq.gz"
            ln -s "$(realpath {input.r2})" "$tmp/{wildcards.lib}_R2.fastq.gz"
            fastqc --threads {threads} --outdir "$tmp" \
                "$tmp/{wildcards.lib}_R1.fastq.gz" "$tmp/{wildcards.lib}_R2.fastq.gz"
            mv "$tmp/{wildcards.lib}_R1_fastqc.zip" {output.r1}
            mv "$tmp/{wildcards.lib}_R2_fastqc.zip" {output.r2}
            rm -rf "$tmp"
            ) > {log} 2>&1
            """


rule tss_sites:
    """Unique transcript TSSs from the GTF (gene TSSs if it has no transcripts)."""
    input:
        gtf=REF["gtf"],
    output:
        "results/reference/tss.bed",
    log:
        "logs/qc/tss_sites.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        set -euo pipefail
        for feature in transcript gene; do
            zcat -f {input.gtf} \
              | awk -F'\t' -v f=$feature 'BEGIN{{OFS="\t"}} $3 == f {{
                    p = ($7 == "-") ? $5 - 1 : $4 - 1; print $1, p, p + 1, ".", ".", $7}}' \
              | LC_ALL=C sort -u -k1,1 -k2,2n -k6,6 > {output}
            [ -s {output} ] && break
        done
        echo "TSSs: $(wc -l < {output})"
        ) > {log} 2>&1
        """


rule tss_enrichment:
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
        tss="results/reference/tss.bed",
    output:
        tsv="results/qc/tss/{lib}.tss_enrichment.tsv",
        profile="results/qc/tss/{lib}.tss_profile.tsv",
    log:
        "logs/qc/{lib}.tss_enrichment.log",
    conda:
        "../envs/deeptools.yaml"
    resources:
        mem_mb=4000,
        runtime=240,
    params:
        flank=DT["tss_flank"],
        edge=DT.get("tss_edge", 100),
        smooth=DT.get("tss_smooth", 20),
        mito=MITO,
    script:
        "../scripts/tss_enrichment.py"


rule frip_library:
    """FRiP of one library over its condition's merged stringent peaks."""
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
        peaks=lambda wildcards: (
            "results/peaks/merged_stringent/"
            f"{LIBRARIES.loc[wildcards.lib, 'condition']}_peaks.narrowPeak"
        ),
    output:
        "results/qc/frip/{lib}.frip.tsv",
    log:
        "logs/qc/{lib}.frip.log",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 4)
    shell:
        """
        (
        set -euo pipefail
        total=$(samtools view -@ {threads} -c -F 2304 {input.bam})
        inpk=0
        if [ -s {input.peaks} ]; then
            inpk=$(samtools view -@ {threads} -c -F 2304 -L {input.peaks} {input.bam})
        fi
        frip=$(awk -v a="$inpk" -v b="$total" 'BEGIN{{ if (b > 0) printf "%.4f", a / b; else print "NA" }}')
        printf "Sample\tReads\tReads_In_Peaks\tFRiP\n%s\t%s\t%s\t%s\n" \
            "{wildcards.lib}" "$total" "$inpk" "$frip" > {output}
        ) > {log} 2>&1
        """


rule qc_summary:
    input:
        unpack(qc_summary_inputs),
    output:
        "results/qc/qc_summary.tsv",
    log:
        "logs/qc/qc_summary.log",
    conda:
        "../envs/python.yaml"
    params:
        libs=LIBS,
        conditions=[LIBRARIES.loc[l, "condition"] for l in LIBS],
        thresholds=QC_THRESHOLDS,
    script:
        "../scripts/qc_summary.py"


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
              *qc_summary.tsv) dest=qc_summary_mqc.tsv ;;
              *) dest="$base" ;;
            esac
            ln -sf "$(realpath "$f")" "$stage/$dest"
        done
        multiqc --force --config {params.config} -o {params.outdir} "$stage"
        rm -rf "$stage"
        ) > {log} 2>&1
        """
