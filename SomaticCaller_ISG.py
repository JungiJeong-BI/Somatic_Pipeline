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


print(read1)
print(read2)
db_path = os.path.dirname(os.path.realpath(__file__)) + "/DB/"

reference_genome =  db_path + "/Homo_sapiens_assembly38.fasta"
#dbsnp = db_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills = db_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
g1000phase1_snp = db_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
g1000phase1_indel = db_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

g1000_pon = db_path + "1000g_pon.hg38.vcf.gz"
gnomad = db_path + "af-only-gnomad.hg38.vcf.gz"
# if you want to define library info on bam file, edit it!
library_name = "SGI"
platform = "Illumina"
platform_barcode = "jjg_Dev"
tmp = db_path+"tmp/"


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
        "ILLUMINACLIP:"+db_path+"/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80"])
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
        "--known-sites", mills, "--known-sites", g1000phase1_snp, "--known-sites", g1000phase1_indel,
        "-R", reference_genome,
        "-O", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.table"])])
    applyrecal = " ".join(["gatk ApplyBQSR -I", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.bam"]),
        "-R", reference_genome, "--bqsr-recal-file", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.table"]),
        "-O", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.bam"])])
    return [bwamapping, addbamrg, sortbam, dedup, recal, applyrecal]



def paired_vaf_call_recalibration(samplename, paired_samplename, outputfolder):
    mutect2_paired_call = " ".join(["gatk Mutect2 -R", reference_genome, "-I",
        "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.bam"]), "-I",
        "".join([outputfolder,"/",paired_samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.bam"]), "-normal", 
        paired_samplename, "-O",
        "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.vcf.gz"])])
    variant_recal_pileup = " ".join(["gatk GetPileupSummaries -I",
        "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.bam"]),
        "-L", mills, "-O", "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.Pileup.table"]),
        "-V", mills])
    variant_recal_contam = " ".join(["gatk CalculateContamination -I",
        "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.Pileup.table"]),
        "-O", "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.Pileup.contam.table"])])
    F1R2_source = " ".join(["gatk CollectF1R2Counts -I", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.bam"]),
        "-O", "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.F1R2.tar.gz"]),
        "-R", reference_genome])
    F1R2_model = " ".join(["gatk LearnReadOrientationModel -I",
        "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.F1R2.tar.gz"]),
        "-O", "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.F1R2.model.tar.gz"])])
    filt_mutect2 = " ".join(["gatk FilterMutectCalls -R", reference_genome, "-V",
        "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.vcf.gz"]),
        "-O", "".join([outputfolder,"/",samplename,"_aligned.trimmomatic.bwa.rg.sorted.dedup.recal.filtered.vcf.gz"]),
        "--contamination-table", "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.Pileup.contam.table"]),
        "-ob-priors", "".join([outputfolder,"/",samplename,"-Analysis-aligned.rg.sorted.dedup.recal.F1R2.model.tar.gz"])])
    return [mutect2_paired_call, variant_recal_pileup, variant_recal_contam, F1R2_source, F1R2_model, filt_mutect2]

Process("mkdir" + tmp)
sample = getsamplename(read1)
trimming(read1, read2, sample , output)
for i_cmd in preprocessing(sample, output):
    Process(i_cmd)
print(read1)
paired_sample = getsamplename(p_read1)
trimming(p_read1, p_read2, paired_sample , output)
for i_cmd in preprocessing(paired_sample, output):
    Process(i_cmd)


for i_cmd in paired_vaf_call_recalibration(sample, paired_sample, output):
    Process(i_cmd)




