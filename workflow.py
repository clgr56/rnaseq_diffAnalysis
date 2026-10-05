import glob
import numpy as np
import subprocess
import config


def fastQC(file):
    f_list = []
    out_dir = f"results/{file}fastqc"
    for fname in glob.glob(f'data/{file}*'):
        f_list.append(fname)
    print("fastQC")
    print(len(f_list))
    for f in f_list:
        print(f)
        subprocess.run(["fastqc", f, "--outdir", out_dir])


def multiQC(file):
    files_dir = f"results/{file}"
    out_dir = f"results/{file}multiqc"
    print("multiQC")
    subprocess.run(["multiqc", files_dir, "--outdir", out_dir])


def fastp(files):
    f_list = []
    out_dir = f"results/{files}fastp_trimmed/"
    for fname in glob.glob(f'data/{files}*1.fq.gz'):
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
        print(out_1)
        print(out_2)
        print(out_report)
        subprocess.run(["fastp", "-i", f, "-I", f2, "-o", out_1, "-O", out_2, "--thread", "16", "--detect_adapter_for_pe", "--trim_poly_x", "-h", out_html, "-j", out_json])


def trim(files):
    print("trim_galore")
    f_list = []
    out_dir = f"results/{files}trimmed/"
    for fname in glob.glob(f'data/{files}*1.fq.gz'):
        f_list.append(fname)
    print("trimming")
    for f in f_list:
        print(f)
        f2 = f[:-7]
        f2 = f2 + '2.fq.gz'
        subprocess.run(["trim_galore", "--paired", f, f2, "-o", out_dir])


def trimmed_fastqc(files):
    print("trimmed_fastqc")
    f_list = []
    out_dir = f"results/{files}fastqc_trimmed/"
    for fname in glob.glob(f'results/{files}fastp_trimmed/*.fq.gz'):
        f_list.append(fname)
    print("fastQC")
    print(len(f_list))
    for f in f_list:
        print(f)
        subprocess.run(["fastqc", f, "--outdir", out_dir])


def iso_quant(files):
    print("salmon indexing")
    subprocess.run(["salmon", "index", "-p", "16", "-t", config.ref_transcript_file, "-i", f"results/{files}salmon/index"], check=True)#ref anpassen
    f_list = []
    for fname in glob.glob(f'results/{files}trimmed/*1.fq.gz'):
        f_list.append(fname)
    print("salmon quantification")
    for f in f_list:
        print(f)
        f2 = f[:-13]
        f2 = f2 + '2_val_2.fq.gz'
        print(f2[36:-12])
        f_base = f'results/{files}salmon/' +f2[36:-12] + '_trimmed'
        subprocess.run(["salmon", "quant", "-p", "16", "-i", f"results/{files}salmon/index", "-l", "A", "-1", f, "-2", f2,"--validateMappings", "-o", f_base], check=True)
        #--validateMapping depricated


def star_align_idx(files):
    print("star_align_idx")
    f_list = []
    for fname in glob.glob(f'results/{files}trimmed/*1.fq.gz'):
        f_list.append(fname)
    print("STAR indexing")   #trying data\ref\gencode.v50.chr_patch_hapl_scaff.annotation.gtf instead of data\ref\gencode.v50.annotation.gtf and w/o "--genomeSAindexNbases", "11",
    subprocess.run(["STAR", "--runMode", "genomeGenerate", "--genomeDir", f"results/{files}star/index/", "--genomeFastaFiles", config.ref_gene_fasta, "--sjdbGTFfile", config.ref_anno_gtf, "--sjdbOverhang", "149", "--runThreadN", '16'], check=True)
    #--sjdbOverhang max.length.read -1


def star_align(files):
    print("star_align")
    cmd = ['unpigz', f'results/{files}trimmed/*.fq.gz']
    subprocess.run(cmd,check=True)
    f_list = []
    for fname in glob.glob(f'results/{files}trimmed/*1.fq'):#trying gz again
        f_list.append(fname)
    print("STAR align")
    for f in f_list:
        print(f)
        f2 = f[:-10] #12
        f2 = f2 + '2_val_2.fq'#.gz
        print(f2)
        f_base = f'results/{files}star/' + f2[36:-11] + '_trimmed' #f2[13:-13]
        RG = 'ID:' + f2[36:-11]
        SM = 'SM:' + f2[36:-11]
        print(RG)
        print(SM)
        print(f_base) #possible sam outout: --outSAMaatributes NH HI AS nM NM MD jM jI MC ch uT and possible unstranded option for cufflinks/cuffdiff: --outSAMstrandField intronMotif  if cufflinks you should remove non-canonical junctions with --outFilterIntronMotifs RemoveNoncanonical
        subprocess.run(["STAR", "--genomeDir", f"results/{files}star/index/", "--runThreadN", '16', "--readFilesIn", f, f2, "--outFileNamePrefix", f_base, "--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within", "--outSAMattributes", "All", "--outSAMattrRGline", RG, SM, "--quantMode", "GeneCounts"],check=True) #'--readFIlesCommand', 'gunzip', '-c',   '--readFilesCommand', 'gunzip', '-c', outSamattributes Standard


def bam_bai(files):
    print("bam_bai")
    f_list = []
    for fname in glob.glob(f"results/{files}star/*.bam"):
        f_list.append(fname)
    for f in f_list:
        cmd = ['samtools', 'index', '-@', '16', '-M', '--bai'] #     Interpret all filename arguments as alignment files to be indexed individually
        cmd.append(f)
        subprocess.run(cmd, check=True)



def sam_depth(files):
    print("sam_depth")
    f_list = []
    out_dir = f"results/{files}sam/"
    for fname in glob.glob(f"results/{files}star/*.bam"):
        f_list.append(fname)
    for f in f_list:
        out_f = f.split('/')[3]
        out = out_dir+out_f[:-4]+'.depth.txt'
        cmd = ['samtools', 'depth', '-@', '16', f, '-o', out]
        subprocess.run(cmd, check=True)



def sam_QC(files):
    print("sam_QC")
    f_list = []
    out_dir = f"results/{files}sam/"
    for fname in glob.glob(f"results/{files}star/*.bam"):
        f_list.append(fname)
    for fil in f_list:
        out_f = fil.split('/')[3]
        out_flagstat = out_dir+out_f[:-4]+'.flagstats.tsv'
        cmd = ['samtools', 'flagstats', '-@', '16', '-O', 'tsv', fil]#, '>', out_flagstat
        with open(out_flagstat, "w",encoding="utf-8") as fi:
            subprocess.call(cmd,stdout=fi)
        #subprocess.run(cmd, check=True)
        out_stats = out_dir+out_f[:-4]+'.stats.txt'
        cmd_stats = ['samtools', 'stats', '-@', '16', '--ref-seq', config.ref_gene_fasta, fil]#, '>', out_stats
        with open(out_stats,"w",encoding="utf-8") as f:
            subprocess.call(cmd_stats, stdout=f)
        #subprocess.run(cmd_stats, check=True)


def picard_markdup(files):
    print("picard_markdup")
    f_list = []
    out_dir = f"results/{files}sam/"
    for fname in glob.glob(f"results/{files}star/*.bam"):
        f_list.append(fname)
    for f in f_list:
        out_f = f.split('/')[3]
        out = out_dir+out_f[:-4]+".markdup.bam"
        out_metrics = out_dir+out_f[:-4]+".metrics.txt"
        cmd = ['picard', 'MarkDuplicates', '--INPUT', f, '--OUTPUT', out, '--METRICS_FILE', out_metrics]
        subprocess.run(cmd, check=True)


def quant_mapper_gene(files):
    print("quant_mapper_gene")
    mode = config.mode#['all']#'UC', 'GM', 'IR'
    for mod in mode:
        f_list = []
        if mod != 'all':
            for fname in glob.glob(f"results/{files}star/*{mod}*.bam"):
                f_list.append(fname)
        else:
            for fname in glob.glob(f'results/{files}star/*.bam'):
                f_list.append(fname)
        string_list = f_list#" ".join(f_list)
        print(string_list)
        outfile = f'results/{files}feature_counts/counts_{mod}_feature.txt'
        cmd = ['featureCounts', '-T', '16', '-p', '--countReadPairs', '-M', '-t', 'exon', '-g', 'gene_id', '-a', config.ref_anno_gtf]
        cmd.extend(['-o', outfile])
        cmd.extend(string_list)
        subprocess.run(cmd, check=True)
        #--countReadPairs for counting features; -M for countMultiMappingReads: ully count every alignment reported for a multi-mapping read (each alignment carries 1 count)
        #-p paired-end;
        outfile2 = f'results/{files}feature_counts/counts_{mod}_reads.txt'
        del cmd[4]
        del cmd[12]
        cmd.insert(12,outfile2)
        subprocess.run(cmd, check=True)


def quant_mapper_gene_all(all):
    print("quant_mapper_gene_all")
    f_list = []
    for l in all:
        for fname in glob.glob(f'results/{l}star/*.bam'):
            f_list.append(fname)
    string_list = f_list#" ".join(f_list)
    print(string_list)
    outfile = 'results/feature_counts/counts_all_feature.txt'
    cmd = ['featureCounts', '-T', '16', '-p', '--countReadPairs', '-M', '-t', 'exon', '-g', 'gene_id', '-a', config.ref_anno_gtf]
    cmd.extend(['-o', outfile])
    cmd.extend(string_list)
    subprocess.run(cmd, check=True)
    #--countReadPairs for counting features; -M for countMultiMappingReads: ully count every alignment reported for a multi-mapping read (each alignment carries 1 count)
    #-p paired-end;
    outfile2 = 'results/feature_counts/counts_all_reads.txt'
    del cmd[4]
    del cmd[12]
    cmd.insert(12,outfile2)
    subprocess.run(cmd, check=True)


def deseq2():
    cmd = ['conda', 'activate', 'pydeseq2']
    subprocess.run(cmd,check=True)
    cmd = ['python3', 'deseq.py']
    subprocess.run(cmd,check=True)
    cmd = ['conda', 'deactivate', 'pydeseq2']
    subprocess.run(cmd,check=True)


def drimseq():
    pass


if __name__=="__main__":
    #for file in config.data_folder:
    file = config.data_folder
    fastQC(file)
    fastp(file)
    trimmed_fastqc(file)
    #multiQC(file)
    trim(file)
    iso_quant(file)
    star_align_idx(file)
    star_align(file)
    bam_bai(file)
    sam_depth(file)
    sam_QC(file)
    picard_markdup(file)
    quant_mapper_gene(file)
    #quant_mapper_gene_all(config.data_folder)
    deseq2()
    multiQC(file)
