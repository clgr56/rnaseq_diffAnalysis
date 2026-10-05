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
    counts_df = pd.read_csv(os.path.join(DATA_PATH, "counts_all_feature.txt"), skiprows=1, sep='\t', index_col=0)#config with counts_{mod}_{reads_feature}
    in_de = counts_df.iloc[:,-12:].T
    in_deseq2 = in_de.filter(regex='^results', axis=0)
    #write file for names, condition, time, batch
    print(in_deseq2)
    condition=['carcinoma','healthy']#config
    sample_names = [f"30_{re.split(r'[/\s.]+',x)[3]}" if "30" in re.split(r'[/\s.]+',x)[1] else f"15_{re.split(r'[/\s.]+',x)[3]}" for x in in_deseq2.index]
    in_deseq2.index = sample_names
    metadata = pd.DataFrame({'condition':condition}, index=sample_names)#config 'seq_batch':seq_batch
    print(metadata)
    met = metadata
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
    print(counts_df)
    #print(counts_df[pd.to_numeric(counts_df['SLR24_A3_trimmmedAligned'],errors='coerce').isna()])
    print(counts_df.shape)
    print(counts_df.index)
    print(counts_df.columns[:10])
    print(counts_df.head())
    print(counts_df.dtypes)
    print(counts_df.dtypes.unique())
    counts_df = counts_df.apply(pd.to_numeric, errors='raise')
    print(counts_df.dtypes.unique())
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design= "~condition",#config
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()
    return dds


def stats(dds):
    ds = DeseqStats(
        dds,
        contrast=['condition', 'carcinoma','healthy'], #confusing namespace in output table
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


def volcano_plot(ds, output_path="results/volcano_condtion.png",padj_cutoff=0.05, lfc_cutoff=1):
    res = ds.results_df.copy()
    # Remove rows without statistics
    res = res.replace([np.inf, -np.inf], np.nan)
    res = res.dropna(subset=["log2FoldChange", "padj"])
    # -log10 adjusted p-value
    res["neg_log10_padj"] = -np.log10(res["padj"].clip(lower=np.finfo(float).tiny))
    # Significance categories
    res["significance"] = "Not significant"
    res.loc[(res["padj"] < padj_cutoff) & (res["log2FoldChange"] > lfc_cutoff),"significance"] = "Higher in A3"
    res.loc[(res["padj"] < padj_cutoff) & (res["log2FoldChange"] < -lfc_cutoff),"significance"] = "Higher in K3"

    # Plot
    plt.figure(figsize=(9, 7))

    for category in ["Not significant", "Higher in A3", "Higher in K3"]:
        subset = res[res["significance"] == category]
        plt.scatter(
            subset["log2FoldChange"],
            subset["neg_log10_padj"],
            s=12,
            alpha=0.6,
            label=category
        )

    # Cutoff lines
    plt.axhline(
        -np.log10(padj_cutoff),
        linestyle="--",
        linewidth=1
    )

    plt.axvline(
        lfc_cutoff,
        linestyle="--",
        linewidth=1
    )

    plt.axvline(
        -lfc_cutoff,
        linestyle="--",
        linewidth=1
    )

    plt.xlabel("log2 fold change (A3 vs K3)")
    plt.ylabel("-log10 adjusted p-value")
    plt.title("DESeq2: A3 vs K3 sequencing model")

    plt.legend()
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Volcano plot saved to: {output_path}")


if __name__=='__main__':
    #save(False)
    data_in, meta = read_in()
    dds = deseq2(data_in, meta)
    ds = stats(dds)
    show_res(dds, ds, meta)
    volcano_plot(ds)
