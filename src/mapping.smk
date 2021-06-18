
rule SampleTrim:
    input:
        r1_fastq = sample_dir + '{sample}'+postfix[0],
        r2_fastq = sample_dir + '{sample}'+postfix_PE[0],
    output:
        r1_paired = 'results/{sample}/1.processed/{sample}'+postfix[0].split("fastq")[0]+'trimmed.paired.fastq.gz',
        r1_unpaired = 'results/{sample}/1.processed/{sample}'+postfix[0].split("fastq")[0]+'trimmed.unpaired.fastq.gz',
        r2_paired = 'results/{sample}/1.processed/{sample}'+postfix_PE[0].split("fastq")[0]+'trimmed.paired.fastq.gz',
        r2_unpaired = 'results/{sample}/1.processed/{sample}'+postfix_PE[0].split("fastq")[0]+'trimmed.unpaired.fastq.gz',
    threads: config["SampleInfo"]["threads"]
    params:
        trimmingStd=config["trimmomatic"]["trimmer"],
        trimmingAdapter=config["Run"]["SequencingInfo"]["adapter"]
    log:
        "results/logs/trimmomatic/{sample}.log"
    message:
        "\n### DSC GATK WES somatic variant pipeline ###\n"\
        "### Trimmomatic start ###\n"\
        "### Error report: jgjeong-intern@insilicogen.com ###\n",
    shell:
        """
        trimmomatic PE -phred33 -threads {threads} \
        {input.r1_fastq} {input.r2_fastq} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        ILLUMINACLIP:{params.trimmingAdapter} \
        {params.trimmingStd} 2> {log}\
        """


rule map_reads:
    input:
        **get_sample_readingtype()
    output:
        temp("results/{sample}/2.aligned/{sample}-Analysis-aligned.bam")
    threads: config["SampleInfo"]["threads"]
    log:
        "results/logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -t {threads} {ref_dir} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule add_readgroups:
    input:
        "results/{sample}/2.aligned/{sample}-Analysis-aligned.bam"
    output:
        temp("results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.bam")
    params:
        LibraryName = config["Run"]["SequencingInfo"]["LibraryName"]+"_{sample}",
        Platform = config["Run"]["SequencingInfo"]["Platform"]+"_{sample}",
        PlatformBarcode = config["Run"]["SequencingInfo"]["Barcode"]+"_{sample}",
        SampleName = config["Run"]["SequencingInfo"]["SampleName"]+"_{sample}"
    log:
        "results/logs/bwa_mem/rg/{sample}.log"
    shell:
        """
        gatk AddOrReplaceReadGroups -I {input} -O {output} \
        -LB {params.LibraryName} -PL {params.Platform} -PU {params.PlatformBarcode} -SM {params.SampleName} \
        --TMP_DIR {tmp_folder} 2> {log}
        """


rule sort_bam:
    input:
        "results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.bam"
    output:
        "results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.bam"
    log:
        "results/logs/samtools/sort/{sample}.log"
    threads: config['SampleInfo']['threads']
    shell:
        """
        samtools sort -@ {threads} {input} \
         -o {output} 2> {log}
        """

rule mark_duplicates:
    input:
        "results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.bam"
    output:
        bam="results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.dedup.bam",
        metrics="results/{sample}/2.aligned/{sample}-Analysis-dedup.metrics"
    threads: config['SampleInfo']['threads']
    log:
        "results/logs/picard/dedup/{sample}.log"
    shell:
        "gatk MarkDuplicatesSpark -I {input} -O {output.bam} -M {output.metrics} --tmp-dir {tmp_folder} --spark-master local[{threads}] 2> {log}"


rule recalibrate_base_qualities:
    input:
        Sample_dedup_bam="results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.dedup.bam",
    output:
        recal_table="results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.dedup.recal.table"
    log:
        "results/logs/gatk/bqsr/{sample}.log"
    params:
        recalibrationDB = "DoTest" #config["params"]["gatk"]["BaseRecalibrator"]
    shell:
        #"gatk BaseRecalibrator -I {input.Sample_dedup_bam} --known-sites {dbsnp} --known-sites {Mills} --known-sites {G1000} --known-sites {gnomad} --known-sites {G1000phase3} --known-sites {Hapmap} -R {ref_dir} -O {output.recal_table}"
        "gatk BaseRecalibrator -I {input.Sample_dedup_bam} --known-sites {dbsnp} --known-sites {Mills} --known-sites {G1000} --known-sites {gnomad} --known-sites {G1000phase3} -R {ref_dir} -O {output.recal_table}"


rule apply_bqsr:
    input:
        bam="results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.dedup.bam",
        recal_table="results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.dedup.recal.table"
    output:
        final_bam=protected("results/{sample}/2.aligned/{sample}-Analysis-aligned.rg.sorted.dedup.recal.bam")
    log:
        "results/logs/gatk/bqsr/applybqsr/{sample}.log"
    shell:
        "gatk ApplyBQSR -I {input.bam} -R {ref_dir} --bqsr-recal-file {input.recal_table} "
        "-O {output.final_bam} 2> {log}"
