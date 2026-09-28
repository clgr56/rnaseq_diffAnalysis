import glob
import numpy as np
import subprocess


def fastQC():
    f_list = []
    out_dir = "results/fastqc"
    for fname in glob.glob('data/SLR24_15MioSeqDepth/*'):
        f_list.append(fname)
    print("fastQC")
    print(len(f_list))
    for f in f_list:
        print(f)
        subprocess.run(["fastqc", f, "--outdir", out_dir])


def multiQC():
    files_dir = "results/"
    out_dir = "results/multiqc"
    print("multiQC")
    subprocess.run(["multiqc", files_dir, "--outdir", out_dir])


def fastp():
    f_list = []
    out_dir = "results/fastp_trimmed/"
    for fname in glob.glob('data/SLR24_15MioSeqDepth/*1.fq.gz'):
        f_list.append(fname)
    print("trimming")
    for f in f_list:
        print(f)
        f2 = f[:-7]
        f2 = f2 + '2.fq.gz'
        out_1 = out_dir+f[:-7].split('/')[2]+"1_trimmed.fq.gz"
        out_2 = out_dir+f[:-7].split('/')[2]+"2_trimmed.fq.gz"
        out_report = out_dir+f[:-7].split('/')[2]+"report"
        out_html = out_report+'.html'
        out_json = out_report+'.json'
        subprocess.run(["fastp", "-i", f, "-I", f2, "-o", out_1, "-O", out_2, "--thread", "16", "--detect_adapter_for_pe", "--trim_poly_x", "-h", out_html, "-j", out_json])


def trim():
    f_list = []
    out_dir = "data/trimmed/"
    for fname in glob.glob('data/SLR24_15MioSeqDepth/*1.fq.gz'):
        f_list.append(fname)
    print("trimming")
    for f in f_list:
        print(f)
        f2 = f[:-7]
        f2 = f2 + '2.fq.gz'
        subprocess.run(["trim_galore", "--paired", f, f2, "-o", out_dir])


def trimmed_fastqc():
    f_list = []
    out_dir = "results/fastqc_trimmed/"
    for fname in glob.glob('results/fastp_trimmed/*.fq.gz'):
        f_list.append(fname)
    print("fastQC")
    print(len(f_list))
    for f in f_list:
        print(f)
        subprocess.run(["fastqc", f, "--outdir", out_dir])


def iso_quant():
    print("salmon indexing")
    subprocess.run(["salmon", "index", "-t", "data/ref/gencode.v50.transcripts.fa", "-i", "results/salmon/index"])
    f_list = []
    for fname in glob.glob('data/trimmed/*1.fq.gz'):
        f_list.append(fname)
    print("salmon quantification")
    for f in f_list:
        print(f)
        f2 = f[:-13]
        f2 = f2 + '2_val_2.fq.gz'
        f_base = 'results/salmon/' +f2[13:-12] + '_trimmmed'
        subprocess.run(["salmon", "quant", "-i", "results/salmon/index", "-l", "A", "-1", f, "-2", f2,"--validateMappings", "-o", f_base])
        #--validateMapping depricated


def star_align_idx():
    f_list = []
    for fname in glob.glob('data/trimmed/*1.fq.gz'):
        f_list.append(fname)
    print("STAR indexing")   #trying data\ref\gencode.v50.chr_patch_hapl_scaff.annotation.gtf instead of data\ref\gencode.v50.annotation.gtf and w/o "--genomeSAindexNbases", "11",
    subprocess.run(["STAR", "--runMode", "genomeGenerate", "--genomeDir", "results/star/index/", "--genomeFastaFiles", "data/ref/GRCh38.p14.genome.fa", "--sjdbGTFfile", "data/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf", "--sjdbOverhang", "149", "--runThreadN", '16'])
    #--sjdbOverhang max.length.read -1


def star_align():
    f_list = []
    for fname in glob.glob('data/trimmed/*1.fq'):#trying gz again
        f_list.append(fname)
    print("STAR align")
    for f in f_list:
        print(f)
        f2 = f[:-10] #12
        f2 = f2 + '2_val_2.fq'#.gz
        print(f2)
        f_base = 'results/star/' + f2[13:-11] + '_trimmmed' #f2[13:-13]
        RG = 'ID:' + f2[13:-11]
        SM = 'SM:' + f2[13:-11]
        print(f_base) #possible sam outout: --outSAMaatributes NH HI AS nM NM MD jM jI MC ch uT and possible unstranded option for cufflinks/cuffdiff: --outSAMstrandField intronMotif  if cufflinks you should remove non-canonical junctions with --outFilterIntronMotifs RemoveNoncanonical
        subprocess.run(["STAR", "--genomeDir", "results/star/index/", "--runThreadN", '16', "--readFilesIn", f, f2, "--outFileNamePrefix", f_base, "--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within", "--outSAMattributes", "All", "--outSAMattrRGline", RG, SM, "--quantMode", "GeneCounts"]) #'--readFIlesCommand', 'gunzip', '-c',   '--readFilesCommand', 'gunzip', '-c', outSamattributes Standard


def bam_bai():
    f_list = []
    for fname in glob.glob(f"results/star/*.bam"):
        f_list.append(fname)
    for f in f_list:
        cmd = ['samtools', 'index', '-M', '--bai', '--threads', '16'] #     Interpret all filename arguments as alignment files to be indexed individually
        cmd.append(f)
        subprocess.run(cmd, check=True)



def sam_depth():
    f_list = []
    out_dir = "results/sam/"
    for fname in glob.glob(f"results/star/*.bam"):
        f_list.append(fname)
    for f in f_list:
        out_f = f.split('/')[2]
        out = out_dir+out_f[:-4]+'.depth.txt'
        cmd = ['samtools', 'depth', f, '-o', out]
        subprocess.run(cmd, check=True)



def sam_QC():
    f_list = []
    out_dir = "results/sam/"
    for fname in glob.glob(f"results/star/*.bam"):
        f_list.append(fname)
    for f in f_list:
        out_f = f.split('/')[2]
        out_flagstat = out_dir+out_f[:-4]+'.flagstats.tsv'
        cmd = ['samtools', 'flagstats', '-@', '16', '-O', 'tsv', f]#, '>', out_flagstat
        with open(out_flagstat, "w",encoding="utf-8") as f:
            subprocess.call(cmd,stdout=f)
        #subprocess.run(cmd, check=True)
        out_stats = out_dir+out_f[:-4]+'.stats.txt'
        cmd_stats = ['samtools', 'stats', '--threads', '16', '--ref-seq', 'data/ref/GRCh38.p14.genome.fa', f]#, '>', out_stats
        with open(out_stats,"w",encoding="utf-8") as f:
            subprocess.call(cmd_stats, stdout=f)
        #subprocess.run(cmd_stats, check=True)


def picard_markdup():
    f_list = []
    out_dir = "results/sam/"
    for fname in glob.glob(f"results/star/*.bam"):
        f_list.append(fname)
    for f in f_list:
        out_f = f.split('/')[2]
        out = out_dir+out_f[:-4]+".markdup.bam"
        out_metrics = out_dir+out_f[:-4]+".metrics.txt"
        cmd = ['picard', 'MarkDuplicates', '--INPUT', f, '--OUTPUT', out, '--METRICS_FILE', out_metrics]
        subprocess.run(cmd, check=True)


def quant_mapper_gene():
    mode = ['all']#'UC', 'GM', 'IR'
    for mod in mode:
        f_list = []
        if mode != 'all':
            for fname in glob.glob(f"results/star/*{mod}*.bam"):
                f_list.append(fname)
        else:
            for fname in glob.glob('results/star/*.bam'):
                f_list.append(fname)
        string_list = f_list#" ".join(f_list)
        print(string_list)
        outfile = f'results/feature_counts/counts_{mod}_feature.txt'
        cmd = ['featureCounts', '-p', '--countReadPairs', '-M', '-t', 'exon', '-g', 'gene_id', '-a', 'data/ref/gencode.v50.annotation.gtf']
        cmd.extend(['-o', outfile])
        cmd.extend(string_list)
        subprocess.run(cmd, check=True)
        #--countReadPairs for counting features; -M for countMultiMappingReads: ully count every alignment reported for a multi-mapping read (each alignment carries 1 count)
        #-p paired-end;
        outfile2 = f'results/feature_counts/counts_{mod}_reads.txt'
        del cmd[2]
        del cmd[10]
        cmd.insert(10,outfile2)
        subprocess.run(cmd, check=True)


def deseq2():
    cmd = ['python3', 'deseq.py']
    subprocess.run(cmd,check=True)


def drimseq():
    pass


if __name__=="__main__":
    #fastQC()
    #fastp()
    #trimmed_fsatqc()
    #multiQC()
    #trim()
    #iso_quant()
    #star_align_idx()
    star_align()
    bam_bai()
    sam_depth()
    sam_QC()
    picard_markdup()
    quant_mapper_gene()
    #deseq2()
    #multiQC()
