# Alignment and filtering (input_mode: fastq), or import of filtered
# merged-library BAMs from an nf-core/atacseq outdir (input_mode: nfcore).
#
# FASTQ mode is a port of the lab's 1.1_atac_align_cc.sh:
#   cutadapt (Nextera, -q 20,20 -m 30 --pair-filter=any)
#   -> bowtie2 --maxins 2000 --no-mixed --no-discordant
#   -> samtools view -q 30 -f 2 -F 2828, drop chrM
#   -> bedtools intersect -v blacklist (optional)
#   -> sort -n -> fixmate -r -m -> drop orphaned mates (-f 2) -> sort
#   -> markdup -r -s -> index, flagstat,
#      fragment-size table, alignment QC table.


# ---------------------------------------------------------------------------
# Reference preparation (always available; cheap)
# ---------------------------------------------------------------------------
rule prepare_fasta:
    input:
        fasta=REF["fasta"],
    output:
        fasta="results/reference/genome.fa",
        fai="results/reference/genome.fa.fai",
        sizes="results/reference/chrom.sizes",
    log:
        "logs/reference/prepare_fasta.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        case "{input.fasta}" in
            *.gz) pigz -dc "{input.fasta}" > {output.fasta} ;;
            *)    ln -sf "$(realpath {input.fasta})" {output.fasta} ;;
        esac
        samtools faidx {output.fasta}
        cut -f1,2 {output.fai} > {output.sizes}
        ) > {log} 2>&1
        """


if not REF.get("bowtie2_index") and INPUT_MODE == "fastq":

    rule bowtie2_build:
        input:
            fasta="results/reference/genome.fa",
        output:
            multiext(
                "results/reference/bowtie2/genome",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        log:
            "logs/reference/bowtie2_build.log",
        conda:
            "../envs/align.yaml"
        threads: threads("bowtie2_build", 8)
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            prefix=lambda wildcards, output: output[0][: -len(".1.bt2")],
        shell:
            "bowtie2-build --threads {threads} {input.fasta} {params.prefix} > {log} 2>&1"


# ---------------------------------------------------------------------------
# FASTQ mode
# ---------------------------------------------------------------------------
if INPUT_MODE == "fastq":

    rule merge_fastq:
        """Concatenate sequencing runs of one library (nf-core semantics)."""
        input:
            merge_fastq_inputs,
        output:
            temp("results/fastq/{lib}_R{read}.merged.fastq.gz"),
        log:
            "logs/align/{lib}_R{read}.merge_fastq.log",
        wildcard_constraints:
            read="1|2",
        conda:
            "../envs/align.yaml"
        shell:
            "cat {input} > {output} 2> {log}"

    rule cutadapt:
        input:
            unpack(trim_inputs),
        output:
            r1=temp("results/fastq/{lib}_R1.trimmed.fastq.gz"),
            r2=temp("results/fastq/{lib}_R2.trimmed.fastq.gz"),
        log:
            "logs/align/{lib}.cutadapt.log",
        conda:
            "../envs/align.yaml"
        threads: threads("cutadapt", 8)
        resources:
            mem_mb=4000,
            runtime=240,
        params:
            a1=config["trimming"]["adapter_r1"],
            a2=config["trimming"]["adapter_r2"],
            q=config["trimming"]["quality"],
            m=config["trimming"]["min_length"],
            extra=config["trimming"].get("extra", ""),
        shell:
            """
            cutadapt -a {params.a1} -A {params.a2} \
                -q {params.q} -m {params.m} --pair-filter=any \
                -j {threads} {params.extra} \
                -o {output.r1} -p {output.r2} \
                {input.r1} {input.r2} > {log} 2>&1
            """

    rule bowtie2_align:
        input:
            r1="results/fastq/{lib}_R1.trimmed.fastq.gz",
            r2="results/fastq/{lib}_R2.trimmed.fastq.gz",
            idx=bowtie2_index_files(),
        output:
            bam=temp("results/bam/tmp/{lib}.aligned.bam"),
        log:
            "logs/align/{lib}.bowtie2.log",
        conda:
            "../envs/align.yaml"
        threads: threads("bowtie2", 16)
        resources:
            mem_mb=16000,
            runtime=720,
        params:
            index=lambda wildcards, input: re.sub(r"\.1\.bt2l?$", "", input.idx[0]),
            maxins=config["align"]["maxins"],
            extra=config["align"].get("bowtie2_extra", ""),
        shell:
            """
            bowtie2 -p {threads} -x {params.index} \
                --maxins {params.maxins} --no-mixed --no-discordant {params.extra} \
                -1 {input.r1} -2 {input.r2} 2> {log} \
            | samtools view -b -@ 2 -o {output.bam} -
            """

    rule filter_bam:
        """MAPQ / proper-pair / flag filter, drop mito reads, optional blacklist."""
        input:
            bam="results/bam/tmp/{lib}.aligned.bam",
            blacklist=blacklist_input(),
        output:
            bam=temp("results/bam/tmp/{lib}.filtered.bam"),
            stats="results/qc/filter_stats/{lib}.filter_stats.tsv",
        log:
            "logs/align/{lib}.filter_bam.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 8)
        resources:
            mem_mb=8000,
            runtime=240,
        params:
            mapq=config["align"]["min_mapq"],
            flag_exclude=config["align"]["flag_exclude"],
            mito=MITO,
        shell:
            """
            (
            set -euo pipefail
            tmp={output.bam}.nomito.bam
            # total records and mito records in the raw alignment
            samtools view -@ {threads} {input.bam} \
              | awk -v m={params.mito} 'BEGIN{{n=0;c=0}} {{n++}} $3==m{{c++}} END{{print n"\\t"c}}' \
              > {output.stats}.counts
            read total mito < {output.stats}.counts
            rm -f {output.stats}.counts

            samtools view -h -@ {threads} -q {params.mapq} -f 2 -F {params.flag_exclude} {input.bam} \
              | awk -v m={params.mito} '$1 ~ /^@/ || $3 != m' \
              | samtools view -b -@ 2 -o $tmp -
            post_mito=$(samtools view -c $tmp)

            if [ -n "{input.blacklist}" ]; then
                bedtools intersect -v -abam $tmp -b {input.blacklist} > {output.bam}
                rm -f $tmp
            else
                mv $tmp {output.bam}
            fi
            post_bl=$(samtools view -c {output.bam})

            printf "total_aligned\\t%s\\nmito_reads\\t%s\\npost_mito_filter\\t%s\\npost_blacklist\\t%s\\nblacklist_removed\\t%s\\n" \
                "$total" "$mito" "$post_mito" "$post_bl" "$((post_mito - post_bl))" > {output.stats}
            ) > {log} 2>&1
            """

    rule library_complexity:
        """ENCODE NRF / PBC1 / PBC2 on the orphan-free BAM before deduplication.

        One fragment per pair (the leftmost mate, TLEN > 0), keyed by
        chrom, start, TLEN and strand:
          NRF  = distinct fragments / all fragments
          PBC1 = fragments seen exactly once / distinct fragments
          PBC2 = fragments seen exactly once / fragments seen exactly twice
        """
        input:
            "results/bam/tmp/{lib}.paired.bam",
        output:
            "results/qc/complexity/{lib}.complexity.tsv",
        log:
            "logs/align/{lib}.complexity.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 4)
        resources:
            mem_mb=8000,
            runtime=240,
        params:
            tmp=lambda wildcards: f"results/bam/tmp/{wildcards.lib}.complexity",
        shell:
            """
            (
            set -euo pipefail
            mkdir -p {params.tmp}
            samtools view -@ {threads} -F 2308 {input} \
              | awk 'BEGIN{{OFS="\\t"}} $9 > 0 {{print $3, $4, $9, int($2 / 16) % 2}}' \
              | LC_ALL=C sort -S 4G --parallel={threads} -T {params.tmp} \
              | uniq -c \
              | awk 'BEGIN{{OFS="\\t"}}
                     {{t += $1; d++; if ($1 == 1) m1++; if ($1 == 2) m2++}}
                     END{{
                       nrf  = (t > 0) ? sprintf("%.4f", d / t) : "NA"
                       pbc1 = (d > 0) ? sprintf("%.4f", m1 / d) : "NA"
                       pbc2 = (m2 > 0) ? sprintf("%.4f", m1 / m2) : "NA"
                       print "Total_Fragments", "Distinct_Fragments", "One_Read", "Two_Reads", "NRF", "PBC1", "PBC2"
                       print t + 0, d + 0, m1 + 0, m2 + 0, nrf, pbc1, pbc2
                     }}' > {output}
            rm -rf {params.tmp}
            ) > {log} 2>&1
            """

    rule fix_pairs:
        """sort -n -> fixmate -r -m -> drop orphans (-f 2) -> sort.

        The per-read MAPQ filter and the blacklist intersect remove reads one at
        a time, so they can leave a mate without its partner. fixmate clears the
        proper-pair flag of such orphans and the -f 2 filter drops them (as the
        ENCODE ATAC pipeline does), so every downstream step sees whole pairs.
        """
        input:
            "results/bam/tmp/{lib}.filtered.bam",
        output:
            bam=temp("results/bam/tmp/{lib}.paired.bam"),
            stats="results/qc/filter_stats/{lib}.pairs.tsv",
        log:
            "logs/align/{lib}.fix_pairs.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 8)
        resources:
            mem_mb=16000,
            runtime=240,
        params:
            tmp=lambda wildcards: f"results/bam/tmp/{wildcards.lib}",
            mem=config["align"].get("sort_mem_per_thread", "768M"),
            flag_exclude=config["align"]["flag_exclude"],
        shell:
            """
            (
            set -euo pipefail
            samtools sort -@ {threads} -m {params.mem} -n -T {params.tmp}.nsort \
                -o {params.tmp}.nsorted.bam {input}
            samtools fixmate -@ {threads} -r -m {params.tmp}.nsorted.bam - \
              | samtools view -b -@ 2 -f 2 -F {params.flag_exclude} \
                  -o {params.tmp}.fixmate.bam -
            before=$(samtools view -c -@ {threads} {params.tmp}.nsorted.bam)
            rm -f {params.tmp}.nsorted.bam
            samtools sort -@ {threads} -m {params.mem} -T {params.tmp}.csort \
                -o {output.bam} {params.tmp}.fixmate.bam
            rm -f {params.tmp}.fixmate.bam
            after=$(samtools view -c -@ {threads} {output.bam})
            printf "pre_pair_filter\t%s\npost_pair_filter\t%s\norphans_removed\t%s\n" \
                "$before" "$after" "$((before - after))" > {output.stats}
            ) > {log} 2>&1
            """

    rule markdup:
        """markdup -r -s on the orphan-free BAM (duplicates removed)."""
        input:
            "results/bam/tmp/{lib}.paired.bam",
        output:
            bam="results/bam/{lib}.final.bam",
            stats="results/qc/markdup/{lib}.markdup.txt",
        log:
            "logs/align/{lib}.markdup.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 8)
        resources:
            mem_mb=8000,
            runtime=240,
        shell:
            """
            samtools markdup -@ {threads} -r -s -f {output.stats} \
                {input} {output.bam} > {log} 2>&1
            """

    rule index_bam:
        input:
            "results/bam/{lib}.final.bam",
        output:
            "results/bam/{lib}.final.bam.bai",
        log:
            "logs/align/{lib}.index.log",
        conda:
            "../envs/align.yaml"
        threads: threads("samtools", 4)
        shell:
            "samtools index -@ {threads} {input} > {log} 2>&1"

    rule fragment_sizes:
        input:
            "results/bam/{lib}.final.bam",
        output:
            "results/qc/fragment_sizes/{lib}_fragment_sizes.tsv",
        log:
            "logs/align/{lib}.fragment_sizes.log",
        conda:
            "../envs/align.yaml"
        shell:
            """
            (samtools view {input} \
              | awk '$9 > 0 && $9 < 1000 {{print $9}}' \
              | sort -n | uniq -c | awk '{{print $2"\\t"$1}}' > {output}) 2> {log}
            """

    rule alignment_qc_report:
        input:
            cutadapt=expand("logs/align/{lib}.cutadapt.log", lib=LIBS),
            bowtie2=expand("logs/align/{lib}.bowtie2.log", lib=LIBS),
            filt=expand("results/qc/filter_stats/{lib}.filter_stats.tsv", lib=LIBS),
            pairs=expand("results/qc/filter_stats/{lib}.pairs.tsv", lib=LIBS),
            markdup=expand("results/qc/markdup/{lib}.markdup.txt", lib=LIBS),
            complexity=expand("results/qc/complexity/{lib}.complexity.tsv", lib=LIBS),
            flagstat=expand("results/qc/flagstat/{lib}.flagstat.txt", lib=LIBS),
            frag=expand("results/qc/fragment_sizes/{lib}_fragment_sizes.tsv", lib=LIBS),
        output:
            "results/qc/alignment_qc_report.tsv",
        log:
            "logs/align/alignment_qc_report.log",
        conda:
            "../envs/python.yaml"
        params:
            libs=LIBS,
        script:
            "../scripts/alignment_qc_report.py"


# ---------------------------------------------------------------------------
# nf-core/atacseq mode: take the filtered merged-library BAMs
# ---------------------------------------------------------------------------
if INPUT_MODE == "nfcore":

    rule import_nfcore_bam:
        """Symlink <aligner>/merged_library/<lib>.mLb.clN.sorted.bam (+ .bai)."""
        input:
            bam=nfcore_lib_bam,
        output:
            bam="results/bam/{lib}.final.bam",
            bai="results/bam/{lib}.final.bam.bai",
        log:
            "logs/import/{lib}.import.log",
        conda:
            "../envs/align.yaml"
        shell:
            """
            (
            ln -sf "$(realpath {input.bam})" {output.bam}
            if [ -s "{input.bam}.bai" ]; then
                ln -sf "$(realpath {input.bam}.bai)" {output.bai}
            else
                samtools index {output.bam}
            fi
            ) > {log} 2>&1
            """


# ---------------------------------------------------------------------------
# Both modes
# ---------------------------------------------------------------------------
rule flagstat:
    input:
        bam="results/bam/{lib}.final.bam",
        bai="results/bam/{lib}.final.bam.bai",
    output:
        "results/qc/flagstat/{lib}.flagstat.txt",
    log:
        "logs/align/{lib}.flagstat.log",
    conda:
        "../envs/align.yaml"
    threads: threads("samtools", 4)
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output} 2> {log}"
