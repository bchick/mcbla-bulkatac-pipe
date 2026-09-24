# IDR + consensus: port of 2.2_atac_idr.sh, generalised to any number of
# replicates.
#
#   sort -k8,8rn -> idr --rank p.value --peak-list <merged relaxed oracle>
#   --idr-threshold 0.05 --plot -> keep column 5 >= 540 (scaled IDR 0.05)
#   -> consensus = cat | cut -f1-3 | sort | bedtools merge
#
# Replicates:
#   2 reps  : rep1 vs rep2 (identical to the source).
#   >2 reps : multi_rep_strategy "encode_max_pair" (default) runs IDR on every
#             pair of true replicates and keeps the pair with the most
#             reproducible peaks (the ENCODE ATAC "conservative" set);
#             "first_two" uses rep1 vs rep2 and warns.
#   1 rep   : IDR skipped; single_rep_fallback "stringent" uses the merged
#             stringent (-q 0.05) peaks for the consensus, "skip" leaves the
#             condition out.

if MODULES.get("idr", True):

    rule sort_rep_peaks:
        input:
            "results/peaks/individual/{lib}_peaks.narrowPeak",
        output:
            temp("results/peaks/idr/sorted/{lib}.sorted.narrowPeak"),
        log:
            "logs/idr/{lib}.sort.log",
        conda:
            "../envs/idr.yaml"
        shell:
            "sort -k8,8rn {input} > {output} 2> {log}"

    rule sort_oracle_peaks:
        input:
            "results/peaks/merged_relaxed/{condition}_peaks.narrowPeak",
        output:
            temp("results/peaks/idr/sorted/{condition}.oracle.sorted.narrowPeak"),
        log:
            "logs/idr/{condition}.sort_oracle.log",
        conda:
            "../envs/idr.yaml"
        shell:
            "sort -k8,8rn {input} > {output} 2> {log}"

    rule idr_pair:
        input:
            a="results/peaks/idr/sorted/{a}.sorted.narrowPeak",
            b="results/peaks/idr/sorted/{b}.sorted.narrowPeak",
            oracle="results/peaks/idr/sorted/{condition}.oracle.sorted.narrowPeak",
        output:
            all="results/peaks/idr/pairs/{condition}/{a}__{b}.idr_all.narrowPeak",
            filt="results/peaks/idr/pairs/{condition}/{a}__{b}.idr.narrowPeak",
            status="results/peaks/idr/pairs/{condition}/{a}__{b}.status",
        log:
            "logs/idr/{condition}/{a}__{b}.idr.log",
        wildcard_constraints:
            a=_alt(LIBS),
            b=_alt(LIBS),
        conda:
            "../envs/idr.yaml"
        resources:
            mem_mb=8000,
            runtime=120,
        params:
            thr=IDR["threshold"],
            scaled=idr_scaled_threshold(),
            rank=IDR.get("rank", "p.value"),
            allow_failure="1" if IDR.get("allow_failure", True) else "0",
        shell:
            """
            if idr --samples {input.a} {input.b} --peak-list {input.oracle} \
                   --input-file-type narrowPeak --rank {params.rank} \
                   --output-file {output.all} --idr-threshold {params.thr} \
                   --plot > {log} 2>&1; then
                awk -v t={params.scaled} 'BEGIN{{OFS="\\t"}} $5 >= t' {output.all} \
                    | cut -f1-10 > {output.filt}
                echo OK > {output.status}
            elif [ "{params.allow_failure}" = 1 ]; then
                echo "WARNING: IDR failed for {wildcards.a} vs {wildcards.b}; recorded as FAILED" >> {log}
                : > {output.all}; : > {output.filt}
                echo FAILED > {output.status}
            else
                exit 1
            fi
            """

    rule idr_select:
        """Pick the reproducible set per condition (max pair when >2 reps)."""
        input:
            pairs=idr_pair_files,
            oracle="results/peaks/merged_relaxed/{condition}_peaks.narrowPeak",
            reps=lambda wildcards: expand(
                "results/peaks/individual/{lib}_peaks.narrowPeak",
                lib=LIBS_BY_COND[wildcards.condition],
            ),
        output:
            peaks="results/peaks/idr/{condition}_idr.narrowPeak",
            stats=temp("results/peaks/idr/{condition}.idr_stats.tsv"),
        log:
            "logs/idr/{condition}.select.log",
        conda:
            "../envs/python.yaml"
        params:
            condition=lambda wildcards: wildcards.condition,
            libs=lambda wildcards: LIBS_BY_COND[wildcards.condition],
            strategy=IDR.get("multi_rep_strategy", "encode_max_pair"),
        script:
            "../scripts/idr_select.py"

    rule idr_summary:
        input:
            stats=expand("results/peaks/idr/{c}.idr_stats.tsv", c=IDR_CONDITIONS),
            single=expand(
                "results/peaks/merged_stringent/{c}_peaks.narrowPeak",
                c=SINGLE_REP_CONDITIONS,
            ),
        output:
            "results/peaks/idr/idr_summary.tsv",
        log:
            "logs/idr/idr_summary.log",
        conda:
            "../envs/python.yaml"
        params:
            single=SINGLE_REP_CONDITIONS,
            fallback=IDR.get("single_rep_fallback", "stringent"),
        script:
            "../scripts/idr_summary.py"

    rule consensus_idr:
        input:
            consensus_inputs,
        output:
            bed="results/peaks/consensus/consensus_idr.bed",
            summary="results/peaks/consensus/consensus_idr.summary.tsv",
        log:
            "logs/idr/consensus.log",
        conda:
            "../envs/idr.yaml"
        shell:
            """
            (
            cat {input} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > {output.bed}
            n=$(wc -l < {output.bed})
            med=$(awk '{{print $3-$2}}' {output.bed} | sort -n \
              | awk '{{a[n++]=$1}} END{{ if (n==0) print "NA"; else if (n%2==1) print a[int(n/2)]; else print (a[n/2-1]+a[n/2])/2 }}')
            printf "Peaks\\tMedian_Width\\n%s\\t%s\\n" "$n" "$med" > {output.summary}
            ) > {log} 2>&1
            """
