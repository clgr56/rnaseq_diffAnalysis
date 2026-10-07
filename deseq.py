import os
import re
import config
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats



def read_in():
    DATA_PATH = f"results/{config.data_folder}feature_counts/"
    counts_df = pd.read_csv(os.path.join(DATA_PATH, "counts_all_feature.txt"), skiprows=1, sep='\t', index_col=0)#config with counts_{mod}_{reads_feature}
    in_de = counts_df.iloc[:,-12:].T
    in_deseq2 = in_de.filter(regex='^results', axis=0)
    #write file for names, condition, time, batch
    print(in_deseq2)
    sample_names = [f"{re.split('/',x)[3][:-37]}" for x in in_deseq2.index]
    in_deseq2.index = sample_names
    metadata = pd.DataFrame(index=sample_names)#config 'seq_batch':seq_batch
    metadata['condition'] = sum(metadata.index.str.extract(r"^(Aza|Gua|Control)").values.tolist(),[])
    print(metadata)
    met = metadata
    i_d = in_deseq2.loc[in_deseq2.index.isin(met.index)]
    print(met)
    print(i_d)
    return i_d, met


def save(SAVE):
    if SAVE:
        OUTPUT_PATH = f"results/{config.data_folder}deseq2/"
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


def stats(dds, group):
    ds = DeseqStats(
        dds,
        contrast=group, #confusing namespace in output table
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
    )
    ds.run_wald_test()
    if ds.independent_filter:
        ds._independent_filtering()
    else:
        ds._p_value_adjustment()
    print(ds.p_values)
    return ds


def show_res(dds, comp, ds, metadata):
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
    plt.savefig(f'results/{config.data_folder}deseq2/VST_PCA_{comp[2]}_{comp[3]}.png',format='png')

    print(ds.p_values)
    print(ds.padj)
    ds.summary()
    re = ds.results_df
    re_sig = re.query("padj<0.05").sort_values("padj")
    print(re_sig)


def volcano_plot(ds, comp, padj_cutoff=0.05, lfc_cutoff=1, label_top=10):
    output_path = f"results/{config.data_folder}deseq2/volcano_plot_{comp[2]}_{comp[3]}.png"
    res = ds.results_df.copy()
    # Remove invalid values
    res = res.replace([np.inf, -np.inf], np.nan)
    res = res.dropna(subset=["log2FoldChange", "padj"])
    # -log10 adjusted p-value
    min_padj = np.finfo(float).tiny
    res["neg_log10_padj"] = -np.log10(res["padj"].clip(lower=min_padj))
    # Significance categories
    res["significance"] = "Not significant"
    up = (res["padj"] < padj_cutoff) & (res["log2FoldChange"] >= lfc_cutoff)
    down = (res["padj"] < padj_cutoff) & (res["log2FoldChange"] <= -lfc_cutoff)
    res.loc[up,"significance"] = f"Higher in {comp[2]}"
    res.loc[down,"significance"] = f"Higher in {comp[3]}"

    # Plot
    plt.figure(figsize=(9, 7))

    categories = ["Not significant", f"Higher in {comp[2]}", f"Higher in {comp[3]}"]

    for category in categories:
        subset = res[res["significance"] == category]
        plt.scatter(
            subset["log2FoldChange"],
            subset["neg_log10_padj"],
            s=14,
            alpha=0.6,
            label=category
        )

    # Cutoff lines
    #FDR Threshold
    plt.axhline(
        -np.log10(padj_cutoff),
        linestyle="--",
        linewidth=1
    )
    #LFC Threshold
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

    # Label top genes
    top_genes = (res[res['padj']< padj_cutoff].sort_values('padj').head(label_top))
    for gene, row in top_genes.iterrowas():
        plt.annotate(str(gene),
                     (
                         row["log2FoldChange"],
                         row["neg_log10_padj"]
                     ),
                     xytext=(5,5),
                     textcoords="offset points",
                     fontsize=8
                     )

    plt.xlabel(f"log2 fold change ({comp[2]} vs {comp[3]})")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(f"DESeq2: {comp[2]} vs {comp[3]} sequencing model")
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    #stats
    print(f"\n{comp[2]} + {comp[3]}")
    print("significatn up:", up.sum())
    print("signifikant down:", down.sum())
    print("Total significants:", (up | down).sum())
    print(f"Volcano plot saved to: {output_path}")


def ma_plot(ds, comp, padj_cutoff = 0.05, lfc_cutoff=1):
    output_path = f"results/{config.data_folder}deseq2/MA_{comp[2]}_{comp[3]}.png"
    res = ds.results_df.copy()
    res = res.replace([np.inf, -np.inf], np.nan)
    res = res.dropna(subset=["log2FoldChange", "baseMean", "padj"])

    significant = ((res["padj"] > padj_cutoff) & (res["log2FoldChange"].abs() >= lfc_cutoff))

    plt.figure(figsize=(9, 7))

    non_sig = res[~significant]

    plt.scatter(
        np.log2(non_sig["baseMean"] + 1),
        non_sig["log2FoldChange"],
        s=12,
        alpha=0.5,
        label="Not significant"
    )

    sig = res[~significant]

    plt.scatter(
            np.log2(sig["baseMean"] + 1),
            sig["log2FoldChange"],
            s=12,
            alpha=0.5,
            label="Not significant"
        )

    
    #zero-line
    plt.axhline(0, linestyle="--", linewidth=1)
    #LFC threshold
    plt.axhline(lfc_cutoff, linestyle="--", linewidth=1)
    plt.axhline(-lfc_cutoff, linestyle="--", linewidth=1)

    plt.xlabel("log2(baseMean + 1)")
    plt.ylabel(f"log2 fold change ({comp[2]} vs {comp[3]})")
    plt.title(f"DESeq2 MA plot: {comp[2]} vs {comp[3]}")
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"MA plot: {significant.sum()} significant genes")


def heatmap_significant_genes(
    dds,
    ds,
    metadata,
    comp,
    padj_cutoff=0.05,
    lfc_cutoff=1.0,
    top_n=100):
    output_path = f"results/{config.data_folder}deseq2/heatmap_{comp[2]}_{comp[3].png}"
    # ------------------------------------------------
    # Get significant DESeq2 genes
    # ------------------------------------------------

    res = ds.results_df.copy()

    res = res.replace(
        [np.inf, -np.inf],
        np.nan
    )

    res = res.dropna(
        subset=[
            "log2FoldChange",
            "padj"
        ]
    )

    significant = res[
        (res["padj"] < padj_cutoff) &
        (res["log2FoldChange"].abs() >= lfc_cutoff)
    ].copy()

    if significant.empty:

        print(
            f"No significant genes for "
            f"{comp[2]} vs {comp[3]}"
        )

        return

    # Sort by adjusted p-value
    significant = significant.sort_values(
        "padj"
    )

    # Limit number of genes
    if top_n is not None:

        significant = significant.head(top_n)

    genes = significant.index.tolist()

    print(
        f"Heatmap {comp[2]} vs {comp[3]}: "
        f"{len(genes)} genes"
    )

    # ------------------------------------------------
    # VST transformation
    # ------------------------------------------------

    # Only calculate VST if not already available
    if "vst_counts" not in dds.layers:

        dds.vst()

    vst = np.asarray(
        dds.layers["vst_counts"]
    )

    # PyDESeq2 stores:
    # rows    = samples
    # columns = genes

    vst_df = pd.DataFrame(
        vst,
        index=dds.obs_names,
        columns=dds.var_names
    )

    # ------------------------------------------------
    # Select significant genes
    # ------------------------------------------------

    heatmap_data = vst_df.loc[
        metadata.index,
        genes
    ]

    # ------------------------------------------------
    # Z-score each gene
    # ------------------------------------------------

    heatmap_z = (
        heatmap_data
        .sub(heatmap_data.mean(axis=0), axis=1)
        .div(heatmap_data.std(axis=0), axis=1)
    )

    # Remove genes with zero variance
    heatmap_z = heatmap_z.dropna(
        axis=1
    )

    # ------------------------------------------------
    # Transpose:
    # genes x samples
    # ------------------------------------------------

    heatmap_z = heatmap_z.T

    # ------------------------------------------------
    # Plot
    # ------------------------------------------------

    fig_width = max(
        8,
        0.8 * len(heatmap_z.columns)
    )

    fig_height = max(
        8,
        0.22 * len(heatmap_z.index)
    )

    plt.figure(
        figsize=(fig_width, fig_height)
    )

    plt.imshow(
        heatmap_z.values,
        aspect="auto",
        interpolation="nearest"
    )

    plt.colorbar(
        label="VST expression\n(row Z-score)"
    )

    plt.xticks(
        range(len(heatmap_z.columns)),
        heatmap_z.columns,
        rotation=90
    )

    plt.yticks(
        range(len(heatmap_z.index)),
        heatmap_z.index,
        fontsize=7
    )

    plt.xlabel("Sample")

    plt.ylabel("Gene")

    plt.title(
        f"Significant genes: "
        f"{comp[2]} vs {comp[3]}"
    )

    plt.tight_layout()

    os.makedirs(
        os.path.dirname(output_path),
        exist_ok=True
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


if __name__=='__main__':
    #save(False)
    data_in, meta = read_in()
    dds = deseq2(data_in, meta)
    group= [["condition", "Aza", "Gua"], ["condition", "Aza", "Control"], ["condition", "Gua", "Control"]]
    for g in group:
        ds = stats(dds, g)
        show_res(dds, g, ds, meta)
        volcano_plot(ds, g)
        ma_plot(ds, g)
        heatmap_significant_genes(dds, ds, meta, g)
