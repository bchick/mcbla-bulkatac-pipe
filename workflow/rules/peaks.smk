# Peak calling: port of 2.1_atac_peaks.sh (four phases).
#   1. merge replicate BAMs per condition
#   2. per-replicate MACS2, relaxed (-p 0.01)            -> IDR input
#   3. per-condition merged MACS2, relaxed (-p 0.01)     -> IDR oracle
#   4. per-condition merged MACS2, stringent (-q 0.05)   -> production
# All narrowPeak files are blacklist-filtered; FRiP + median width per file.

PK = config["peaks"]


rule merge_condition_bam:
    input:
        cond_bam_inputs,
    output:
        bam="results/bam/merged/{condition}.merged.bam",
        bai="results/bam/merged/{condition}.merged.bam.bai",
    log:
        "logs/peaks/{condition}.merge.log",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 8)
    resources:
        mem_mb=4000,
        runtime=240,
    params:
        mode=lambda wildcards: (
            "link" if nfcore_mrp_bam(wildcards.condition) else "merge"
        ),
    shell:
        """
        (
        if [ "{params.mode}" = link ]; then
            ln -sf "$(realpath {input})" {output.bam}
        else
            samtools merge -f -@ {threads} {output.bam} {input}
        fi
        samtools index -@ {threads} {output.bam}
        ) > {log} 2>&1
        """


rule macs2_individual:
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
        blacklist=blacklist_input(),
    output:
        peaks="results/peaks/individual/{lib}_peaks.narrowPeak",
        xls="results/peaks/individual/{lib}_peaks.xls",
        summits="results/peaks/individual/{lib}_summits.bed",
    log:
        "logs/peaks/individual/{lib}.macs2.log",
    conda:
        "../envs/macs2.yaml"
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.peaks),
        name=lambda wildcards: wildcards.lib,
        threshold=f"-p {PK['relaxed_pvalue']}",
        gsize=REF["macs2_gsize"],
        extra=PK.get("macs2_extra", ""),
    shell:
        MACS2_SHELL


rule macs2_merged:
    input:
        bam="results/bam/merged/{condition}.merged.bam",
        bai="results/bam/merged/{condition}.merged.bam.bai",
        blacklist=blacklist_input(),
    output:
        peaks="results/peaks/{mergedset}/{condition}_peaks.narrowPeak",
        xls="results/peaks/{mergedset}/{condition}_peaks.xls",
        summits="results/peaks/{mergedset}/{condition}_summits.bed",
    log:
        "logs/peaks/{mergedset}/{condition}.macs2.log",
    wildcard_constraints:
        mergedset="merged_relaxed|merged_stringent",
    conda:
        "../envs/macs2.yaml"
    resources:
        mem_mb=16000,
        runtime=480,
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.peaks),
        name=lambda wildcards: wildcards.condition,
        threshold=lambda wildcards: (
            f"-p {PK['relaxed_pvalue']}"
            if wildcards.mergedset == "merged_relaxed"
            else f"-q {PK['stringent_qvalue']}"
        ),
        gsize=REF["macs2_gsize"],
        extra=PK.get("macs2_extra", ""),
    shell:
        MACS2_SHELL


rule peak_stats:
    """Peaks, FRiP (reads in peaks / reads, -F 2304) and median width."""
    input:
        peaks="results/peaks/{peakset}/{name}_peaks.narrowPeak",
        bam=peak_stats_bam,
    output:
        temp("results/peaks/{peakset}/{name}.stats.tsv"),
    log:
        "logs/peaks/{peakset}/{name}.stats.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        n=$(wc -l < {input.peaks})
        frip=NA; median=NA
        if [ "$n" -gt 0 ]; then
            total=$(samtools view -c -F 2304 {input.bam})
            inpk=$(samtools view -c -F 2304 -L {input.peaks} {input.bam})
            frip=$(awk -v a="$inpk" -v b="$total" 'BEGIN{{ if (b>0) printf "%.4f", a/b; else print "NA" }}')
            median=$(awk '{{print $3-$2}}' {input.peaks} | sort -n \
              | awk '{{a[n++]=$1}} END{{ if (n%2==1) print a[int(n/2)]; else print (a[n/2-1]+a[n/2])/2 }}')
        fi
        printf "%s\\t%s\\t%s\\t%s\\n" "{wildcards.name}" "$n" "$frip" "$median" > {output}
        ) > {log} 2>&1
        """


rule peak_summary:
    input:
        peakset_members,
    output:
        "results/peaks/{peakset}/peak_summary.tsv",
    log:
        "logs/peaks/{peakset}/peak_summary.log",
    conda:
        "../envs/align.yaml"
    params:
        label=lambda wildcards: (
            "Sample" if wildcards.peakset == "individual" else "Condition"
        ),
    shell:
        """
        (printf "{params.label}\\tPeaks\\tFRiP\\tMedian_Width\\n"; cat {input}) > {output} 2> {log}
        """


rule union_stringent:
    """Union of merged-stringent peaks over all conditions (bedtools merge)."""
    input:
        expand("results/peaks/merged_stringent/{c}_peaks.narrowPeak", c=CONDITIONS),
    output:
        "results/peaks/consensus/union_stringent.bed",
    log:
        "logs/peaks/union_stringent.log",
    conda:
        "../envs/align.yaml"
    shell:
        "(cat {input} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > {output}) 2> {log}"
