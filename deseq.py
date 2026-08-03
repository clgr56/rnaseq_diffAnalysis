import os
#import pickle as okl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats



def read_in():
    DATA_PATH = "results/feature_counts/"
    counts_df = pd.read_csv(os.path.join(DATA_PATH, "counts_all_reads.txt"), skiprows=1, sep='\t', index_col=0)
    in_deseq2 = counts_df.iloc[:,-12:].T
    print(in_deseq2)
    condition = ['B','C','C', 'A', 'A', 'B', 'C', 'B', 'A', 'A', 'A', 'A'] #A=IR, B=GM, C=UC
    time = ['Y', 'Y', 'Y', 'X', 'Y', 'Y', 'Y', 'Y', 'Y', 'X', 'Y', 'X'] #x=3h, y=24h
    batch = ['pA', 'pC', 'pB', 'pA', 'pA', 'pB', 'pA', 'pC', 'pC', 'pB', 'pB', 'pC']
    sample_names = in_deseq2.index
    metadata = pd.DataFrame({'condition':condition,'time':time, 'batch':batch}, index=sample_names)
    print(metadata)
    met = metadata.loc[(metadata['time']=='Y')] # without IR 3h
    del met['time']
    i_d = in_deseq2.loc[in_deseq2.index.isin(met.index)]
    print(met)
    print(i_d)
    return i_d, met # in_deseq2, metadata


def save(SAVE):
    if SAVE:
        OUTPUT_PATH = "results/deseq2/"
        os.makedirs(OUTPUT_PATH, exist_ok=True)


def deseq2(data_in, metadata):
    counts_df = data_in
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design= "~batch + condition",
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()
    return dds


def stats(dds):
    ds = DeseqStats(
        dds,
        contrast=['condition', 'A', 'B'], #nimmt irgwie trotzdem B und C????
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
    )
    ds.run_wald_test()
    print(ds.p_values)
    if ds.independent_filter:
        ds._independent_filtering()
    else:
        ds._p_value_adjustment()
    return ds


def show_res(dds, ds, metadata):
    print(dds)
    print(dds.var["dispersions"])
    print(dds.varm["LFC"])
    dds.vst()
    vst = dds.layers["vst_counts"]

    x = np.asarray(vst)
    pca = PCA(n_components=3).fit_transform(x)  # PCs across samples

    # Plot, color by condition
    cond = metadata["condition"].values
    plt.scatter(pca[:, 0], pca[:, 1], c=(cond == cond[0]))
    for i, sid in enumerate(metadata.index):
        plt.text(pca[i, 0], pca[i, 1], sid, fontsize=8)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("VST-PCA by condition")
#    plt.show()
    plt.savefig('results/deseq2/VST-PCA.png',format='png')

    print(ds.p_values)
    print(ds.padj)
    ds.summary()
    #ds.plot_MA(s=20)
    re = ds.results_df
    re_sig = re.query("padj<0.05").sort_values("padj")#.head()
    print(re_sig)


if __name__=='__main__':
    #save(False)
    data_in, meta = read_in()
    dds = deseq2(data_in, meta)
    ds = stats(dds)
    show_res(dds, ds, meta)
