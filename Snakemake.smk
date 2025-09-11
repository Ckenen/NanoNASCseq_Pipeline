#!/usr/bin/env runsnakemake

RUNS = ["20220729_K562_1"]
RUN_CELLS = []
for run in RUNS:
    for line in open("data/%s.tsv" % run):
        cell = line.strip("\n").split("\t")[0]
        RUN_CELLS.append("%s/%s" % (run, cell))

GENOME_FASTA = "data/reference/GRCh38.canonical.genome.fa"
SPLICE_MMI = "data/reference/GRCh38.canonical.mm2.splice.mmi"
TRANSCRIPT_BED = "data/reference/gencode.v39.transcripts.bed"
TRANSCRIPT_RICH_NAME_BED = "data/reference/gencode.v39.transcripts.rich_name.bed.gz"
SNP_BED = "data/snp151.3.lite.bed.gz"
OUTDIR = "results"

rule all:
    input:
        # expand(OUTDIR + "/1_fbilr/{run}.tsv.gz", run=RUNS),
        # expand(OUTDIR + "/2_splitted/{run}", run=RUNS),
        # expand(OUTDIR + "/3_trimmed/{run_cell}", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/4_minimap2/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/5_filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/6_extract_umi/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/7_stat_clip/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/8_mark_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/9_remove_duplicate/{run_cell}.bam", run_cell=RUN_CELLS), # (optional)
        # expand(OUTDIR + "/10_marked_events/{run_cell}.bam", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/11_ratio_all/{run_cell}.tsv", run_cell=RUN_CELLS), # (optional)
        # expand(OUTDIR + "/12_ratio_rmdup/{run_cell}.tsv", run_cell=RUN_CELLS), # (optional)
        # expand(OUTDIR + "/13_ratio_corrected/{run_cell}.tsv", run_cell=RUN_CELLS),
        # expand(OUTDIR + "/14_marked_new/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/15_quant_tgs/{run_cell}", run_cell=RUN_CELLS),
        


rule fbilr:
    input:
        fq = "data/{run}.fastq.gz",
        fa = "data/{run}.fa"
    output:
        txt = OUTDIR + "/1_fbilr/{run}.tsv.gz"
    log:
        OUTDIR + "/1_fbilr/{run}.log"
    threads:
        48
    shell:
        """
        ( fbilr -t {threads} -r 10000 -m PE -b {input.fa} {input.fq} | gzip -c > {output.txt} ) &> {log}
        """

rule split_reads:
    input:
        fq = "data/{run}.fastq.gz",
        mtx = rules.fbilr.output.txt,
        txt = "data/{run}.tsv"
    output:
        out = directory(OUTDIR + "/2_splitted/{run}")
    log:
        OUTDIR + "/2_splitted/{run}.log"
    threads:
        12
    shell:
        """
        ./scripts/split_reads.py {input} {output} &> {log}
        pigz -p {threads} {output}/*/*.fastq
        """

rule trim_reads:
    input:
        fqs = rules.split_reads.output.out
    output:
        out = directory(OUTDIR + "/3_trimmed/{run}/{cell}")
    log:
        OUTDIR + "/3_trimmed/{run}/{cell}.log"
    params:
        fq = rules.split_reads.output.out + "/succeed/{cell}.fastq.gz"
    shell:
        """
        ./scripts/trim_reads.py {params.fq} {output.out} &> {log}
        gzip {output.out}/trimmed.fastq
        """

rule minimap2:
    input:
        fqd = rules.trim_reads.output.out,
        mmi = SPLICE_MMI,
        bed = TRANSCRIPT_BED,
    output:
        bam = OUTDIR + "/4_minimap2/{run}/{cell}.bam",
        bai = OUTDIR + "/4_minimap2/{run}/{cell}.bam.bai",
        txt = OUTDIR + "/4_minimap2/{run}/{cell}.flagstat",
    log:
        OUTDIR + "/4_minimap2/{run}/{cell}.log"
    params:
        rg = '@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}'
    threads:
        12
    shell:
        """(
        minimap2 -ax splice -u f -Y --MD -R '{params.rg}' -t {threads} \
            --junc-bed {input.bed} {input.mmi} {input.fqd}/trimmed.fastq.gz \
            | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {output.bam} -o {output.bam} - 
        samtools index -@ {threads} {output.bam} 
        samtools flagstat -@ {threads} {output.bam} > {output.txt} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/5_filtered/{run}/{cell}.bam",
        bai = OUTDIR + "/5_filtered/{run}/{cell}.bam.bai",
        txt = OUTDIR + "/5_filtered/{run}/{cell}.flagstat",
    log:
        OUTDIR + "/5_filtered/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} --expr 'rname =~ "^chr([0-9]+|[XY])$"' \
            -q 30 -m 200 -F 2308 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} 
        samtools flagstat -@ {threads} {output.bam} > {output.txt} ) &> {log}
        """

rule extract_umi:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = OUTDIR + "/6_extract_umi/{run}/{cell}.bam",
        bai = OUTDIR + "/6_extract_umi/{run}/{cell}.bam.bai",
    log:
        OUTDIR + "/6_extract_umi/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        ./scripts/extract_umi.py {input.bam} {output.bam} 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule stat_clip:
    input:
        bam = rules.extract_umi.output.bam
    output:
        bam = OUTDIR + "/7_stat_clip/{run}/{cell}.bam",
        bai = OUTDIR + "/7_stat_clip/{run}/{cell}.bam.bai",
        txt = OUTDIR + "/7_stat_clip/{run}/{cell}.flagstat",
        tsv = OUTDIR + "/7_stat_clip/{run}/{cell}.tsv",
    log:
        OUTDIR + "/7_stat_clip/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        nasctools StatClip -c 20 -s {output.tsv} -i {input.bam} -o {output.bam}
        samtools index -@ {threads} {output.bam} 
        samtools flagstat -@ {threads} {output.bam} > {output.txt} ) &> {log}
        """

rule mark_duplicate: 
    input:
        bam = rules.stat_clip.output.bam
    output:
        bam = OUTDIR + "/8_mark_duplicate/{run}/{cell}.bam",
        bai = OUTDIR + "/8_mark_duplicate/{run}/{cell}.bam.bai",
        txt = OUTDIR + "/8_mark_duplicate/{run}/{cell}.flagstat",
        tsv = OUTDIR + "/8_mark_duplicate/{run}/{cell}.tsv"
    log:
        OUTDIR + "/8_mark_duplicate/{run}/{cell}.log"
    shell:
        """(
        nasctools MarkDuplicate -s {output.tsv} -i {input.bam} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.txt} ) &> {log}
        """

rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = OUTDIR + "/9_remove_duplicate/{run}/{cell}.bam",
        bai = OUTDIR + "/9_remove_duplicate/{run}/{cell}.bam.bai",
        txt = OUTDIR + "/9_remove_duplicate/{run}/{cell}.flagstat",
    log:
        OUTDIR + "/9_remove_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -F 1024 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.txt} ) &> {log}
        """

rule mark_event:
    input:
        bam = rules.mark_duplicate.output.bam,
        bed = SNP_BED,
        fa = GENOME_FASTA,
    output:
        bam = OUTDIR + "/10_marked_events/{run}/{cell}.bam",
        bai = OUTDIR + "/10_marked_events/{run}/{cell}.bam.bai"
    log:
        OUTDIR + "/10_marked_events/{run}/{cell}.log"
    threads:
        8
    shell:
        """(
        nasctools MarkEvent -c -t {threads} -f {input.fa} \
            -s {input.bed} -i {input.bam} -o {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule report_mismatch_all:
    input:
        bam = rules.mark_event.output.bam
    output:
        tsv = OUTDIR + "/11_ratio_all/{run}/{cell}.tsv"
    log:
        OUTDIR + "/11_ratio_all/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools ReportMismatch -t {threads} \
            -i {input.bam} -o {output.tsv} &> {log}
        """

rule report_mismatch_rmdup:
    input:
        bam = rules.mark_event.output.bam
    output:
        tsv = OUTDIR + "/12_ratio_rmdup/{run}/{cell}.tsv"
    log:
        OUTDIR + "/12_ratio_rmdup/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        nasctools ReportMismatch --discard-duplicates -t {threads} \
            -i {input.bam} -o {output.tsv} &> {log}
        """

rule report_mismatch_consensus:
    input:
        bam = rules.mark_event.output.bam
    output:
        tsv = OUTDIR + "/13_ratio_corrected/{run}/{cell}.tsv"
    log:
        OUTDIR + "/13_ratio_corrected/{run}/{cell}.log"
    shell:
        """
        nasctools ReportMismatch --event-tag HE --reads 2 --discard-duplicates \
            -t {threads} -i {input.bam} -o {output.tsv} &> {log}
        """

rule mark_new:
    input:
        bam = rules.mark_event.output.bam
    output:
        bam = OUTDIR + "/14_marked_new/{run}/{cell}.bam",
        bai = OUTDIR + "/14_marked_new/{run}/{cell}.bam.bai",
    log:
        OUTDIR + "/14_marked_new/{run}/{cell}.log"
    shell:
        """
        nasctools MarkNew --event-tag HE --tc 2 --quality 0 \
            -i {input.bam} -o {output.bam} &> {log}
        samtools index {output.bam}
        """

rule quant_tgs:
    input:
        bam = rules.mark_new.output.bam,
        bed = TRANSCRIPT_RICH_NAME_BED,
    output:
        out = directory(OUTDIR + "/15_quant_tgs/{run}/{cell}"),
    log:
        OUTDIR + "/15_quant_tgs/{run}/{cell}.log"
    shell:
        """
        nasctools QuantTGS -t {input.bed} -i {input.bam} -s 2 -c 2 -o {output.out} &> {log}
        """
