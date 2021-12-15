import os
import sys
import subprocess
import re
import argparse

Get_Options = argparse.ArgumentParser()
Get_Options.add_argument("-r1","--read1", required = True)
Get_Options.add_argument("-r2","--read2", required = True)
Get_Options.add_argument("-pr1","--p_read1", required = False)
Get_Options.add_argument("-pr2","--p_read2", required = False)

Get_Options.add_argument("-t","--thread", required = True)
Get_Options.add_argument("-o","--output", required = True)

Options = Get_Options.parse_args()
read1 = Options.read1; read2 = Options.read2; p_read1=Options.p_read1; p_read2=Options.p_read2; thread=str(Options.thread); output=Options.output

maxcpu = 40

outputlog_path = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/fastq/out/temp/log.log"
out = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/fastq/out/"
tmp = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/fastq/out/temp/"

print(read1)
print(read2)

reference_genome =  "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/gatkgoogle/Homo_sapiens_assembly38.fasta"
dbsnp = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/gatkgoogle/Homo_sapiens_assembly38.dbsnp138.vcf"
mills = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/gatkgoogle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
g1000phase1_snp = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/gatkgoogle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
g1000phase1_indel = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/gatkgoogle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
g1000phase3 = "/BiO2/users/jjg/tools/Somatic_Pipeline_DB/gatkgoogle/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf"

# if you want to define library info on bam file, edit it!
library_name = "SGI"
platform = "Illumina"
platform_barcode = "jjg_Dev"



def Process(cmd):
    print(cmd)
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        stdout, stderr = process.communicate()
    except:
        print('error')


def getsamplename(filepath):
    samplename = filepath.split("/")[-1].split(".")[0]
    return samplename

def trimming(read1, read2, samplename, outputfolder):
    trimmomatic = " ".join(["trimmomatic PE -phred33 -threads", thread,
        read1, read2,
        "".join([outputfolder,"/",samplename,"_R1.trimmed.paired.fastq.gz"]), "".join([outputfolder,"/",samplename,"_R1.trimmed.unpaired.fastq.gz"]),
        "".join([outputfolder,"/",samplename,"_R2.trimmed.paired.fastq.gz"]), "".join([outputfolder,"/",samplename,"_R2.trimmed.unpaired.fastq.gz"]),
        "ILLUMINACLIP:/BiO2/users/jjg/tools/Somatic_Pipeline/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80"])
    print(trimmomatic)
    Process(trimmomatic)

def preprocessing(samplename, outputfolder):
    bwamapping = " ".join(["bwa mem -t", thread, reference_genome, 
        "".join([outputfolder,"/",samplename,"_R1.trimmed.paired.fastq.gz"]), "".join([outputfolder,"/",samplename,"_R2.trimmed.paired.fastq.gz"]),
        "| samtools view -Sb - > ", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.bam"])])
    addbamrg = " ".join(["gatk AddOrReplaceReadGroups -I", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.bam"]),
    "-O", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.bam"]),
        "-LB", library_name, "-PL", platform, "-PU", platform_barcode, "-SM ", samplename,
        "--TMP_DIR", tmp])
    sortbam = " ".join(["samtools sort -@", thread, "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.bam"]), "-o", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.bam"])])
    dedup = " ".join(["gatk MarkDuplicatesSpark -I", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.bam"]),
        "-O", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.bam"]),
        "-M", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.metrics"]),
        "--tmp-dir", tmp, "--spark-master local[" + thread + "]"])
    recal = " ".join(["gatk BaseRecalibrator -I", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.bam"]),
        "--known-sites", dbsnp, "--known-sites", mills, "--known-sites", g1000phase1_snp, "--known-sites", g1000phase1_indel,
        "--known-sites", g1000phase3, "-R", reference_genome,
        "-O", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.table"])])
    applyrecal = " ".join(["gatk ApplyBQSR -I", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.bam"]),
        "-R", reference_genome, "--bqsr-recal-file", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.table"]),
        "-O", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.bam"])])
    return [bwamapping, addbamrg, sortbam, dedup, recal, applyrecal]

print(p_read1)
if p_read1!=None and p_read2!=None:
    sample = getsamplename(read1)
    #trimming(read1, read2, sample , output)
    for i_cmd in preprocessing(sample, output):
        Process(i_cmd)
else:
    paired_sample = getsamplename(p_read1)


























