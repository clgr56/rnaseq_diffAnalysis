import os
import re
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
    #write file for names, condition, time, batch
    print(in_deseq2)
    condition = ['B','C','C', 'A', 'A', 'B', 'C', 'B', 'A', 'A', 'A', 'A'] #A=IR, B=GM, C=UC
    #time = ['Y', 'Y', 'Y', 'X', 'Y', 'Y', 'Y', 'Y', 'Y', 'X', 'Y', 'X'] #x=3h, y=24h
    #batch = ['pA', 'pC', 'pB', 'pA', 'pA', 'pB', 'pA', 'pC', 'pC', 'pB', 'pB', 'pC']
    sample_names = [re.split(r"[/\s.]+",x)[2] for x in in_deseq2.index]
    in_deseq2.index = sample_names
    metadata = pd.DataFrame({'condition':condition,'time':time, 'batch':batch}, index=sample_names)
    print(metadata)
    #ignore 3h IR
    met = metadata.loc[(metadata['time']=='Y')] # without IR 3h
    del met['time']
    i_d = in_deseq2.loc[in_deseq2.index.isin(met.index)]
    print(met)
    print(i_d)
    return i_d, met


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
        contrast=['condition', 'A', 'B'], #confusing namespace in output table
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
    #make own method for plotting etc
    print(dds)
    print(dds.var["dispersions"])
    print(dds.varm["LFC"])
    dds.vst()
    vst = dds.layers["vst_counts"]

    x = np.asarray(vst)
    pca_model = PCA(n_components=3)
    pca = pca_model.fit_transform(x)  # PCAs across samples
    print(f"Explained variance ratio: {pca_model.explained_variance_ratio_}")
    print(f"Singular values: {pca_model.singular_values_}")
    # Plot, color by condition, w/o batch norm
    cond = metadata["condition"].values
    plt.scatter(pca[:, 0], pca[:, 1], c=(cond == cond[0]))
    for i, sid in enumerate(metadata.index):
        plt.text(pca[i, 0], pca[i, 1], sid, fontsize=8)
    plt.xlabel(f"PC1 ({pca_model.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca_model.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title("VST-PCA by condition")
    plt.savefig('results/deseq2/VST-PCA_2.png',format='png')

    print(ds.p_values)
    print(ds.padj)
    ds.summary()
    re = ds.results_df
    re_sig = re.query("padj<0.05").sort_values("padj")
    print(re_sig)


if __name__=='__main__':
    #save(False)
    data_in, meta = read_in()
    dds = deseq2(data_in, meta)
    ds = stats(dds)
    show_res(dds, ds, meta)
