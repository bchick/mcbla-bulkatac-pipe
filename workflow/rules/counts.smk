# Counting (built only for the analyses that are switched on).
#   * featureCounts (fragments, -p --countReadPairs) over the peak set named in
#     timecourse.peaks / chromvar.peaks (as in the lab's
#     05_atac_temporal_clustering/00_preprocessing.Rmd), one matrix per set:
#     results/counts/<peak set>/.
#   * DiffBind sample sheet + dba.count(minOverlap): the counted DBA object is
#     shared by the diff and normcheck modules, so the same count matrix is
#     reused under every normalization (as in 08b, which reused 08's DBA).

CNT = config["counts"]
DIFF = config["diff"]


rule peaks_saf:
    input:
        lambda wildcards: COUNTSETS[wildcards.countset],
    output:
        "results/counts/{countset}/peaks.saf",
    log:
        "logs/counts/{countset}.saf.log",
    conda:
        "../envs/counts.yaml"
    params:
        min_width=CNT["min_width"],
        max_width=CNT["max_width"],
    shell:
        """
        (printf "GeneID\\tChr\\tStart\\tEnd\\tStrand\\n"
         awk -v lo={params.min_width} -v hi={params.max_width} 'BEGIN{{OFS="\\t"}}
             /^(#|track|browser)/ {{next}}
             {{w=$3-$2}} w>=lo && w<=hi {{print $1":"$2+1"-"$3, $1, $2+1, $3, "."}}' {input}
        ) > {output} 2> {log}
        """


rule featurecounts:
    input:
        saf="results/counts/{countset}/peaks.saf",
        bams=all_lib_bams(),
    output:
        raw="results/counts/{countset}/featurecounts.txt",
        summary="results/counts/{countset}/featurecounts.txt.summary",
        matrix="results/counts/{countset}/peak_counts.tsv",
    log:
        "logs/counts/{countset}.featurecounts.log",
    conda:
        "../envs/counts.yaml"
    threads: threads("featurecounts", 8)
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        libs=LIBS,
    shell:
        """
        (
        featureCounts -F SAF -a {input.saf} -o {output.raw} -p --countReadPairs \
            -T {threads} {input.bams}
        {{ printf "peak_id"; for l in {params.libs}; do printf "\\t%s" "$l"; done; echo;
          tail -n +3 {output.raw} | cut -f1,7-; }} > {output.matrix}
        ) > {log} 2>&1
        """


rule diffbind_samplesheet:
    input:
        bams=all_lib_bams(STAT_LIBS),
        peaks=expand("results/peaks/individual/{lib}_peaks.narrowPeak", lib=STAT_LIBS),
    output:
        "results/counts/diffbind_samplesheet.csv",
    log:
        "logs/counts/diffbind_samplesheet.log",
    conda:
        "../envs/python.yaml"
    params:
        libs=STAT_LIBS,
        conditions=[LIBRARIES.loc[l, "condition"] for l in STAT_LIBS],
        treatments=[
            LIBRARIES.loc[l, "treatment"] or LIBRARIES.loc[l, "condition"]
            for l in STAT_LIBS
        ],
        replicates=[LIBRARIES.loc[l, "replicate"] for l in STAT_LIBS],
        tissue=DIFF.get("tissue", "NA"),
    script:
        "../scripts/diffbind_samplesheet.py"


rule diffbind_count:
    input:
        sheet="results/counts/diffbind_samplesheet.csv",
        bams=all_lib_bams(STAT_LIBS),
        consensus=DIFF_PEAKS[1] if DIFF_PEAKS and DIFF_PEAKS[1] else [],
    output:
        "results/counts/dba_counted.rds",
    log:
        "logs/counts/diffbind_count.log",
    conda:
        "../envs/r.yaml"
    threads: threads("diffbind", 8)
    resources:
        mem_mb=64000,
        runtime=720,
    params:
        min_overlap=DIFF["min_overlap"],
        summits=DIFF.get("summits"),
        peakset=DIFF_PEAKS[0] if DIFF_PEAKS else "",
    script:
        "../scripts/diffbind_count.R"
