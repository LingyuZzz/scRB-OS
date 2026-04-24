import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from gseapy import enrichr

def load_data(cancer_type):
    prefix = f"data/{cancer_type.lower()}_mal"
    meta = pd.read_csv(f"{prefix}_meta.csv", index_col=0)
    mat = sc.read_mtx(f"{prefix}_data.mtx").T
    genes = pd.read_csv(f"{prefix}_var.csv")
    cells = pd.read_csv(f"{prefix}_obs.csv")
    mat.var_names = genes['gene']
    mat.obs_names = cells['cell']
    mat.obs = meta.loc[mat.obs_names]
    
    adata = sc.AnnData(X=mat.X, obs=mat.obs, var=mat.var)
    adata.obs = meta.loc[adata.obs_names]
    
    auc_mtx = pd.read_csv(f'data/scenic/{cancer_type.lower()}_tf_activity_matrix.csv', index_col=0)
    auc_mtx = auc_mtx.loc[adata.obs_names]
    adata.obsm['X_aucell'] = auc_mtx.values
    adata.uns['aucell_regulons'] = auc_mtx.columns.tolist()
    
    return adata, auc_mtx

def correlate_tfs_with_module(adata, auc_df, module_col="ME_module1", alpha=0.05):
    module_scores = adata.obs[module_col]
    results = []
    for tf in auc_df.columns:
        r, p_val = spearmanr(module_scores, auc_df[tf], nan_policy="omit")
        results.append({"TF": tf, "rho": r, "p_value": p_val})
    
    result_df = pd.DataFrame(results)
    result_df["fdr"] = fdrcorrection(result_df["p_value"])[1]
    result_df["significant"] = result_df["fdr"] < alpha
    return result_df.sort_values("rho", ascending=False)

def split_pos_neg(corr_df, corr_col="rho", tf_col="TF"):
    pos = set(corr_df.loc[(corr_df[corr_col] > 0.3) & (corr_df['fdr'] < 0.05), tf_col])
    neg = set(corr_df.loc[(corr_df[corr_col] < -0.3) & (corr_df['fdr'] < 0.05), tf_col])
    return pos, neg

adata_os, auc_mtx_os = load_data('os')
adata_rb, auc_mtx_rb = load_data('rb')

os_me_tf_corr = correlate_tfs_with_module(adata_os, auc_mtx_os, module_col="ME_module1")
os_rs_tf_corr = correlate_tfs_with_module(adata_os, auc_mtx_os, module_col="RS_module1")
rb_me_tf_corr = correlate_tfs_with_module(adata_rb, auc_mtx_rb, module_col="ME_module1")
rb_rs_tf_corr = correlate_tfs_with_module(adata_rb, auc_mtx_rb, module_col="RS_module1")

os_me_pos, os_me_neg = split_pos_neg(os_me_tf_corr)
rb_me_pos, rb_me_neg = split_pos_neg(rb_me_tf_corr)
os_rs_pos, os_rs_neg = split_pos_neg(os_rs_tf_corr)
rb_rs_pos, rb_rs_neg = split_pos_neg(rb_rs_tf_corr)

me_pos_shared = os_me_pos & rb_me_pos
me_neg_shared = os_me_neg & rb_me_neg
rs_pos_shared = os_rs_pos & rb_rs_pos
rs_neg_shared = os_rs_neg & rb_rs_neg

shared_dict = {
    "ME_pos_shared": list(me_pos_shared),
    "ME_neg_shared": list(me_neg_shared),
    "RS_pos_shared": list(rs_pos_shared),
    "RS_neg_shared": list(rs_neg_shared),
}

if not os.path.exists("results"): os.makedirs("results")
for key, genes in shared_dict.items():
    pd.DataFrame(genes, columns=["TF"]).to_csv(f"results/{key}.csv", index=False)

try:
    os_adjacencies = pd.read_csv('data/scenic/os_adj.csv')
    rb_adjacencies = pd.read_csv('data/scenic/rb_adj.csv')

    me_pos_tfs = {tf.replace("(+)", "") for tf in me_pos_shared}
    rs_pos_tfs = {tf.replace("(+)", "") for tf in rs_pos_shared}
    all_tfs = list(me_pos_tfs) + list(rs_pos_tfs)

    common_targets = {}
    for tf in all_tfs:
        rb_targets = set(rb_adjacencies[rb_adjacencies['TF'] == tf]['target'].tolist())
        os_targets = set(os_adjacencies[os_adjacencies['TF'] == tf]['target'].tolist())
        common_targets[tf] = list(rb_targets.intersection(os_targets))

    go_results = {}
    for tf, targets in common_targets.items():
        if len(targets) > 0:
            enr = enrichr(gene_list=targets,
                          gene_sets=['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'GO_Cellular_Component_2023'],
                          organism='Human', cutoff=0.05)
            go_results[tf] = enr.results
            go_results[tf].to_csv(f"results/module_tf_GO_{tf}.csv", index=False)

    def build_matrix(go_res, category, alpha=0.05):
        tf_terms = {}
        for tf, df in go_res.items():
            terms = df.loc[(df["Gene_set"] == category) & (df["Adjusted P-value"] < alpha), "Term"].tolist()
            tf_terms[tf] = terms
        all_terms = sorted(set(term for terms in tf_terms.values() for term in terms))
        mat = pd.DataFrame(0, index=tf_terms.keys(), columns=all_terms)
        for tf, terms in tf_terms.items():
            mat.loc[tf, terms] = 1
        return mat

    bp_matrix = build_matrix(go_results, "GO_Biological_Process_2023")
    sim = cosine_similarity(bp_matrix)
    sim_df = pd.DataFrame(sim, index=bp_matrix.index, columns=bp_matrix.index)
    
    if not os.path.exists("plots"): os.makedirs("plots")
    sns.clustermap(sim_df, cmap="viridis", metric="euclidean", method="ward", figsize=(6,6))
    plt.title("Overall similarity (BP)")
    plt.savefig("plots/module_TF_similarity_clustermap.pdf", dpi=300, bbox_inches="tight")
except Exception as e:
    print("Could not complete target enrichment analysis:", e)
