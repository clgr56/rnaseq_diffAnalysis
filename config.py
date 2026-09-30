data_folder=["SLR24_15MioSeqDepth/","SLR24_15MioSeqDepth/"]
ref_transcript_file="data/ref/gencode.v50.transcripts.fa"
ref_gene_fasta="data/ref/GRCh38.p14.genome.fa"
ref_anno_gtf="data/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf"
#star_align: max(length)-1
sjdbOverhang=149
#add_read_group=["--outSAMattrRGline", RG, SM,]
#overal config
threads=16
#feature_counts
#['UC', 'GM', 'IR', 'all']
mode=['A3','K3','all']
#deseq
##guy
#condition = ['B','C','C', 'A', 'A', 'B', 'C', 'B', 'A', 'A', 'A', 'A'] #A=IR, B=GM, C=UC
#time = ['Y', 'Y', 'Y', 'X', 'Y', 'Y', 'Y', 'Y', 'Y', 'X', 'Y', 'X'] #x=3h, y=24h
#batch = ['pA', 'pC', 'pB', 'pA', 'pA', 'pB', 'pA', 'pC', 'pC', 'pB', 'pB', 'pC']
#metadata = pd.DataFrame({'condition':condition,'time':time, 'batch':batch}, index=sample_names)
#met = metadata.loc[(metadata['time']=='Y')]  # without IR 3h #ignore 3h IR
#design= "~batch + condition"
##15
#condition=['carcinoma','healthy']
#metadata = pd.DataFrame({'condition':condition}, index=sample_names)
#design= "~condition"
##30