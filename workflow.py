import glob
import subprocess


def fastQC():
    f_list = []
    out_dir = "results/fastqc"
    for fname in glob.glob('data/Berlin_AKI_tubuloid_bulk_RNA_seq/*'):
        f_list.append(fname)
    print("fastQC")
    print(len(f_list))
    for f in f_list:
        print(f)
        subprocess.run(["fastqc", f, "--outdir", out_dir])


def multiQC():
    files_dir = "results/"#fastqc/
    out_dir = "results/multiqc"
    print("multiQC")
    subprocess.run(["multiqc", files_dir, "--outdir", out_dir])


def fastp():
    f_list = []
    out_dir = "results/fastp_trimmed/"
    for fname in glob.glob('data/Berlin_AKI_tubuloid_bulk_RNA_seq/*1.fq.gz'):
        f_list.append(fname)
    print("trimming")
    for f in f_list:
        print(f)
        f2 = f[:-7]
        f2 = f2 + '2.fq.gz'
        out_1 = out_dir+f[:-7].split('/')[2]+"1_trimmed.fq.gz"
        out_2 = out_dir+f[:-7].split('/')[2]+"2_trimmed.fq.gz"
        subprocess.run(["fastp", "-i", f, "-I", f2, "-o", out_1, "-O", out_2, "--thread", 8, "--detect_adaüter_for-pe", "--trim_poly_x"])


def trim():
    f_list = []
    out_dir = "data/trimmed/"
    for fname in glob.glob('data/Berlin_AKI_tubuloid_bulk_RNA_seq/*1.fq.gz'):
        f_list.append(fname)
    print("trimming")
    for f in f_list:
        print(f)
        f2 = f[:-7]
        f2 = f2 + '2.fq.gz'
        subprocess.run(["trim_galore", "--paired", f, f2, "-o", out_dir])


def trimmed_fsatqc():
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

def align():
    f_list = []
    for fname in glob.glob('data/trimmed/*1.fq'):
        f_list.append(fname)
    print("STAR indexing")
    subprocess.run(["STAR", "--runMode", "genomeGenerate", "--genomeDir", "results/star/index/", "--genomeFastaFiles", "data/ref/GRCh38.p14.genome.fa", "--sjdbGTFfile", "data/ref/gencode.v50.annotation.gtf", "--sjdbOverhang", "149", "--genomeSAindexNbases", "11"])
    #--sjdbOverhang max.length.read -1
    print("STAR align")
    for f in f_list:
        print(f)
        f2 = f[:-10]
        f2 = f2 + '2_val_2.fq'
        print(f2)
        f_base = 'results/star/' + f2[13:-11] + '_trimmmed'
        print(f_base)
        subprocess.run(["STAR", "--genomeDir", "results/star/index/", "--runThreadN", '8', "--readFilesIn", f, f2, "--outFileNamePrefix", f_base, "--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within", "--outSAMattributes", "Standard", "--quantMode", "GeneCounts"])


def quant_mapper_gene():
    mode = ['UC', 'GM', 'IR', 'all']
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
    #do smth with pydeseq2
    pass


def drimseq():
    pass


if __name__=="__main__":
    fastQC()
    fastp()
    trimmed_fsatqc()
    multiQC()
    #trim()
    #iso_quant()
    #align()
    #quant_mapper_gene()
