import os
import re
import config
import pandas as pd
import numpy as np
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats



def read_in():
    data_path = f"results/{config.data_folder}feature_counts/"
    counts_df = pd.read_csv(os.path.join(data_path, "counts_all_feature.txt"), skiprows=1, sep='\t', index_col=0)#config with counts_{mod}_{reads_feature}
    print(counts_df.iloc[:,:12].T)
    in_de = counts_df.iloc[:,-12:].T
    in_deseq2 = in_de.filter(regex='^results', axis=0)
    #write file for names, condition, time, batch
    print(in_deseq2)
    sample_names = [f"{re.split('/',x)[3][:-37]}" for x in in_deseq2.index]
    in_deseq2.index = sample_names
    metadata = pd.DataFrame(index=sample_names)#config 'seq_batch':seq_batch
    metadata['condition'] = sum(metadata.index.str.extract(r"^(Aza|Gua|Control)").values.tolist(),[])
    print(metadata)
    met = metadata.sort_index()
    i_d = in_deseq2.loc[in_deseq2.index.isin(met.index)]
    i_d = i_d.loc[met.index]
    print(met)
    print(i_d)
    return i_d, met


def save(res, comp):
    output_dir = f"results/{config.data_folder}deseq2/"
    os.makedirs(output_dir, exist_ok=True)
    res_df = res.results_df.copy()
    contrast_name = f"{comp[1]}_vs_{comp[2]}"
    res_path  = os.path.join(output_dir, f"DeSeq2_{contrast_name}.tsv")
    res_df.to_csv(res_path, sep='\t', index=True, index_labels="gene_id")


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
    plt.scatter(pca[:, 0], pca[:, 1], c= cond == cond[0])
    for i, sid in enumerate(metadata.index):
        plt.text(pca[i, 0], pca[i, 1], sid, fontsize=8)
    plt.xlabel(f"PC1 ({pca_model.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca_model.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title("VST-PCA by condition")
    plt.savefig(f'results/{config.data_folder}deseq2/VST_PCA_{comp[1]}_{comp[2]}.png',format='png')

    print(ds.p_values)
    print(ds.padj)
    ds.summary()
    re = ds.results_df
    re_sig = re.query("padj<0.05").sort_values("padj")
    print(re_sig)


def volcano_plot(ds, comp, padj_cutoff=0.05, lfc_cutoff=1, label_top=10):
    output_path = f"results/{config.data_folder}deseq2/volcano_plot_{comp[1]}_{comp[2]}.png"
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
    res.loc[up,"significance"] = f"Higher in {comp[1]}"
    res.loc[down,"significance"] = f"Higher in {comp[2]}"

    # Plot
    plt.figure(figsize=(9, 7))

    categories = ["Not significant", f"Higher in {comp[1]}", f"Higher in {comp[2]}"]

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
    for gene, row in top_genes.iterrows():
        plt.annotate(str(gene),
                     (
                         row["log2FoldChange"],
                         row["neg_log10_padj"]
                     ),
                     xytext=(5,5),
                     textcoords="offset points",
                     fontsize=8
                     )

    plt.xlabel(f"log2 fold change ({comp[1]} vs {comp[2]})")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(f"DESeq2: {comp[1]} vs {comp[2]} sequencing model")
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    #stats
    print(f"\n{comp[1]} + {comp[2]}")
    print("significatn up:", up.sum())
    print("signifikant down:", down.sum())
    print("Total significants:", (up | down).sum())
    print(f"Volcano plot saved to: {output_path}")


def ma_plot(ds, comp, padj_cutoff = 0.05, lfc_cutoff=1):
    output_path = f"results/{config.data_folder}deseq2/MA_{comp[1]}_{comp[2]}.png"
    res = ds.results_df.copy()
    res = res.replace([np.inf, -np.inf], np.nan)
    res = res.dropna(subset=["log2FoldChange", "baseMean", "padj"])

    significant = ((res["padj"] < padj_cutoff) & (res["log2FoldChange"].abs() >= lfc_cutoff))
    non_sig = res[~significant]
    sig = res[significant]

    plt.figure(figsize=(9, 7))

    plt.scatter(
        np.log2(non_sig["baseMean"] + 1),
        non_sig["log2FoldChange"],
        s=12,
        alpha=0.5,
        label="Not significant"
    )
    plt.scatter(
            np.log2(sig["baseMean"] + 1),
            sig["log2FoldChange"],
            s=12,
            alpha=0.5,
            label="significant"
        )

    #zero-line
    plt.axhline(0, linestyle="--", linewidth=1)
    #LFC threshold
    plt.axhline(lfc_cutoff, linestyle="--", linewidth=1)
    plt.axhline(-lfc_cutoff, linestyle="--", linewidth=1)

    plt.xlabel("log2(baseMean + 1)")
    plt.ylabel(f"log2 fold change ({comp[1]} vs {comp[2]})")
    plt.title(f"DESeq2 MA plot: {comp[1]} vs {comp[2]}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"MA plot: {significant.sum()} significant genes")


def ma_plot_pds(ds,comp):
    output_path = f"results/{config.data_folder}deseq2/MA_pds_{comp[1]}_{comp[2]}.png"
    ds.plot_MA(log=True, save_path=output_path, s=15)
    print(f"MA_pds plot saved to: {output_path}")


def dist_heatmap(dds, metadata):
    output_path = f"results/{config.data_folder}deseq2/dist_heatmap.png"
# Only calculate VST if not already available
    if "vst_counts" not in dds.layers:
        print("calculating VST")
        dds.vst()
    vst = np.asarray(dds.layers["vst_counts"])

    # PyDESeq2 stores:
    # rows    = samples
    # columns = genes

    vst_df = pd.DataFrame(
        vst,
        index=dds.obs_names,
        columns=dds.var_names
    )

    sample_metadata = metadata.loc[vst_df.index]

    samp_dist = pdist(vst_df.to_numpy(), metric="euclidean")
    samp_dist_matrix = squareform(samp_dist)

    sample_labels = [
        f"{sample_metadata.loc[s, "condition"]}-{s}"
        for s in vst_df.index
    ]

    samp_dist_df = pd.DataFrame(
        samp_dist_matrix,
        index=sample_labels,
        columns=sample_labels
    )

    hc = linkage(samp_dist, method="complete")

    grid = sns.clustermap(
        samp_dist_df,
        row_linkage=hc,
        col_linkage=hc,
        cmap="Blues",
        figsize=(10, 10),
        linewidths=0.3,
        dendrogram_ratio=(0.15, 0.15),
        cbar_kws={"label": "Euclidean distance"}
    )

    grid.ax_heatmap.set_title(
        "Heatmap of sample-to-sample distances",
        pad=15,
    )

    grid.fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(grid.fig)

    print(f"Sample distance heatmap saved: {output_path}")


def heatmap_significant_genes(
    dds,
    ds,
    metadata,
    comp,
    padj_cutoff=0.05,
    lfc_cutoff=1.0,
    top_n=100,
    linkage_method="complete",
    distance_metric="euclidean"):
    output_path = f"results/{config.data_folder}deseq2/heatmap_{comp[1]}_{comp[2]}.png"
    # ------------------------------------------------
    # Get significant DESeq2 genes
    # ------------------------------------------------

    res = ds.results_df.copy()
    res = res.replace([np.inf, -np.inf], np.nan)
    res = res.dropna(subset=["log2FoldChange", "padj"])

    significant = res.loc[(res["padj"] < padj_cutoff) & (res["log2FoldChange"].abs() >= lfc_cutoff)].sort_values("padj")
    genes = significant.index.tolist()
    print(
        f"Heatmap {comp[1]} vs {comp[2]}: "
        f"{len(genes)} genes"
    )

    if top_n is not None:
        significant = significant.head(top_n)
    if significant.empty:
        print(f"No significant genes for {comp[1]} vs {comp[2]}")
        return

    # ------------------------------------------------
    # VST transformation
    # ------------------------------------------------

    # Only calculate VST if not already available
    if "vst_counts" not in dds.layers:
        print("calculating VST")
        dds.vst()

    # PyDESeq2 stores:
    # rows    = samples
    # columns = genes

    vst_df = pd.DataFrame(
        vst = np.asarray(dds.layers["vst_counts"]),
        index=dds.obs_names,
        columns=dds.var_names
    )

    # 3. Align metadata and select genes
    # Preserve the order of samples in the PyDESeq2 object.
    missing_samples = vst_df.index.difference(metadata.index)
    if len(missing_samples):
        raise ValueError(
            f"Samples missing from metadata: {missing_samples.tolist()}"
        )

    sample_metadata = metadata.loc[vst_df.index].copy()

    genes = [
        gene for gene in significant.index
        if gene in vst_df.columns
    ]

    if not genes:
        print(
            "No significant result gene IDs match the VST matrix. "
            "Check that both use the same gene identifiers."
        )
        return
    # ------------------------------------------------
    # Select significant genes
    # ------------------------------------------------

    heatmap_data = vst_df.loc[:, genes].copy()

    # ------------------------------------------------
    # Z-score each gene
    # ------------------------------------------------

    # 4. Z-score each gene across samples
    gene_means = heatmap_data.mean(axis=0)
    gene_stds = heatmap_data.std(axis=0, ddof=1)

    # Exclude constant or invalid genes before clustering
    valid_genes = gene_stds[
        np.isfinite(gene_stds) & (gene_stds > 0)
    ].index

    heatmap_data = heatmap_data.loc[:, valid_genes]

    if heatmap_data.empty:
        print("No variable genes remain for the heatmap.")
        return

    if heatmap_data.empty:
        print("No variable genes remain for the heatmap.")
        return

    heatmap_z = (
        heatmap_data - gene_means.loc[valid_genes]
    ) / gene_stds.loc[valid_genes]

    # Transpose: genes = rows, samples = columns
    heatmap_z = heatmap_z.T

    if heatmap_z.shape[0] < 2 or heatmap_z.shape[1] < 2:
        print("At least two genes and two samples are needed for clustering.")
        return

    # 5. Calculate hierarchical clustering from expression profiles
    # Row clustering: distances between genes
    gene_distances = pdist(
        heatmap_z.to_numpy(),
        metric=distance_metric,
    )
    gene_linkage = linkage(
        gene_distances,
        method=linkage_method,
    )

    # Column clustering: distances between samples
    sample_distances = pdist(
        heatmap_z.to_numpy().T,
        metric=distance_metric,
    )
    sample_linkage = linkage(
        sample_distances,
        method=linkage_method,
    )

    # 6. Add condition colors above the heatmap
    condition_colors = {
        "Aza": "#4C72B0",
        "Gua": "#DD8452",
        "Control": "#55A868",
    }

    col_colors = sample_metadata["condition"].map(condition_colors)
    if col_colors.isna().any():
        missing_conditions = (
            sample_metadata.loc[col_colors.isna(), "condition"]
            .unique()
            .tolist()
        )
        raise ValueError(
            f"No color assigned to conditions: {missing_conditions}"
        )

    # 7. Plot heatmap with both dendrograms
    grid = sns.clustermap(
        heatmap_z,
        row_linkage=gene_linkage,
        col_linkage=sample_linkage,
        row_cluster=True,
        col_cluster=True,
        col_colors=col_colors,
        cmap="vlag",
        center=0,
        vmin=-2,
        vmax=2,
        figsize=(
            max(10, 0.8 * heatmap_z.shape[1]),
            max(8, 0.22 * heatmap_z.shape[0]),
        ),
        linewidths=0,
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.15, 0.15),
        cbar_kws={"label": "VST expression (gene-wise Z-score)"},
        tree_kws={"linewidths": 1.0},
    )

    grid.ax_heatmap.set_xlabel("Sample")
    grid.ax_heatmap.set_ylabel("Gene")
    grid.ax_heatmap.set_title(
        f"Significant genes: {comp[1]} vs {comp[2]}",
        pad=15,
    )
    grid.ax_heatmap.tick_params(axis="x", labelrotation=90)
    grid.ax_heatmap.tick_params(axis="y", labelsize=7)

    legend_handles = [
        Patch(facecolor=color, label=condition)
        for condition, color in condition_colors.items()
        if condition in sample_metadata["condition"].values
    ]

    grid.ax_heatmap.legend(
        handles=legend_handles,
        title="Condition",
        loc="upper left",
        bbox_to_anchor=(1.35, 1.0),
        borderaxespad=0,
    )

    grid.fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(grid.fig)

    print(
        f"Clustered heatmap saved: {output_path} "
        f"({heatmap_z.shape[0]} genes, "
        f"{heatmap_z.shape[1]} samples)"
    )


if __name__=='__main__':
    #save(False)
    data_in, meta = read_in()
    dds = deseq2(data_in, meta)
    group= [["condition", "Aza", "Gua"], ["condition", "Aza", "Control"], ["condition", "Gua", "Control"]]
    dist_heatmap(dds,meta)
    for g in group:
        ds = stats(dds, g)
        save(ds, g)
        show_res(dds, g, ds, meta)
        volcano_plot(ds, g)
        ma_plot(ds, g)
        ma_plot_pds(ds,g)
        heatmap_significant_genes(dds, ds, meta, g)
        dist_heatmap(dds, meta)

