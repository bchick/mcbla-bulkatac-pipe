# TF footprinting: port of 10_tf_footprinting/scripts/01-03.
#   TOBIAS ATACorrect per condition on the merged BAM over footprint.peaks
#   -> TOBIAS ScoreBigwig -> one TOBIAS BINDetect across footprint.conditions
#   (empty = every condition).
# Motifs: footprint.motifs (JASPAR format). If empty, JASPAR2020 CORE
# vertebrates (the chromVAR motif set) is exported with TFBSTools.

FP = config["footprint"]

if RUN_FOOTPRINT:

    if not FP.get("motifs"):

        rule export_jaspar_motifs:
            output:
                "results/reference/jaspar_motifs.jaspar",
            log:
                "logs/footprint/export_jaspar_motifs.log",
            conda:
                "../envs/r.yaml"
            params:
                collection=CV["jaspar_collection"],
                tax_group=CV["jaspar_tax_group"],
            script:
                "../scripts/export_jaspar.R"

    rule tobias_atacorrect:
        input:
            bam="results/bam/merged/{condition}.merged.bam",
            bai="results/bam/merged/{condition}.merged.bam.bai",
            fasta="results/reference/genome.fa",
            fai="results/reference/genome.fa.fai",
            peaks=FP_PEAKS[1],
            blacklist=blacklist_input(),
        output:
            corrected="results/footprint/atacorrect/{condition}_corrected.bw",
            bias="results/footprint/atacorrect/{condition}_bias.bw",
            expected="results/footprint/atacorrect/{condition}_expected.bw",
            uncorrected="results/footprint/atacorrect/{condition}_uncorrected.bw",
        log:
            "logs/footprint/{condition}.atacorrect.log",
        conda:
            "../envs/tobias.yaml"
        threads: threads("tobias", 8)
        resources:
            mem_mb=32000,
            runtime=720,
        params:
            outdir=lambda wildcards, output: os.path.dirname(output.corrected),
            blacklist=lambda wildcards, input: (
                f"--blacklist {input.blacklist}" if input.blacklist else ""
            ),
            extra=FP.get("atacorrect_extra", ""),
        shell:
            """
            TOBIAS ATACorrect --bam {input.bam} --genome {input.fasta} \
                --peaks {input.peaks} {params.blacklist} --outdir {params.outdir} \
                --cores {threads} --prefix {wildcards.condition} {params.extra} \
                > {log} 2>&1
            """

    rule tobias_scorebigwig:
        input:
            signal="results/footprint/atacorrect/{condition}_corrected.bw",
            peaks=FP_PEAKS[1],
        output:
            "results/footprint/footprintscores/{condition}_footprints.bw",
        log:
            "logs/footprint/{condition}.scorebigwig.log",
        conda:
            "../envs/tobias.yaml"
        threads: threads("tobias", 8)
        resources:
            mem_mb=16000,
            runtime=480,
        params:
            extra=FP.get("scorebigwig_extra", ""),
        shell:
            """
            TOBIAS ScoreBigwig --signal {input.signal} --regions {input.peaks} \
                --output {output} --cores {threads} {params.extra} > {log} 2>&1
            """

    rule tobias_bindetect:
        input:
            signals=expand(
                "results/footprint/footprintscores/{c}_footprints.bw", c=FP_CONDITIONS
            ),
            motifs=tobias_motifs(),
            fasta="results/reference/genome.fa",
            fai="results/reference/genome.fa.fai",
            peaks=FP_PEAKS[1],
        output:
            "results/footprint/bindetect/bindetect_results.txt",
        log:
            "logs/footprint/bindetect.log",
        conda:
            "../envs/tobias.yaml"
        threads: threads("tobias_bindetect", 32)
        resources:
            mem_mb=64000,
            runtime=1440,
        params:
            outdir=lambda wildcards, output: os.path.dirname(output[0]),
            conditions=FP_CONDITIONS,
            extra=FP.get("bindetect_extra", ""),
        shell:
            """
            TOBIAS BINDetect --signals {input.signals} --cond-names {params.conditions} \
                --motifs {input.motifs} --genome {input.fasta} --peaks {input.peaks} \
                --outdir {params.outdir} --cores {threads} {params.extra} > {log} 2>&1
            """
