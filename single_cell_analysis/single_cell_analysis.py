import scanpy as sc
import seaborn as sns
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colormaps as cm
import numpy as np
from scipy.sparse import issparse
import pandas as pd
import os

# adata raw counts post qc
counts_postqc = adata.layers['counts_postqc'] # compressed sparse matrix
counts_postqc = counts_postqc.toarray()

# extract number of cells expressing gene X
ebf2 = adata.var_names.get_loc('Ebf2') # get ebf2 index
ebf2_count = adata.layers['counts_postqc'][:, ebf2].toarray().flatten()
(ebf2_count > 0).sum() # 6773 cells expressing ebf2
ebf2_count.sum() # number of transcripts of ebf2 across all cells, 20796

adata.X = adata.layers['counts_postqc']

# cell types annotations
markers = {
        'Fibroblast progenitors': ['Hic1', 'Pdgfra'],
        'Mitotic cells': ['Brd4', 'Ctcf'],
        'Dermal Fibroblast': ['Twist2', 'Dpt', 'Crabp1'],
        'Neural progenitors': ['Sox2', 'Sox3', 'Lhx3'],
        'Smooth Muscle': ['Cnn1', 'Acta2'],
        #'Muscle connective tissue': ['Ngfr', 'Osr1', 'Osr2'],
        'Mast cells': ['Kit', 'Srgn'],
        'Neural cells': ['Stmn2', 'Stmn3', 'Hoxb8'],
        'Macrophages': ['Fcer1g', 'C1qa'],
        'Epidermis': ['Krt5', 'Krt14'],
        'Cartilage': ['Cnmd', 'Col2a1'],
        'Meninges': ['Foxc1', 'Cldn11', 'Aldh1a2'],
        'Endothelial Cells': ['Cdh5', 'Kdr'],
        'Schwann Cells': ['Sox10', 'Mpz'],
        'Skeletal Muscle': ['Myog', 'Ttn', 'Actn2', 'Myh3'],
        'Red Blood Cells': ['Hba-a1', 'Hba-a2']
         }

marker_genes = [
        ['Hic1', 'Pdgfra'],
        ['Brd4', 'Ctcf'],
        ['Twist2', 'Dpt', 'Crabp1'],
        ['Sox2', 'Sox3', 'Lhx3'],
        ['Cnn1', 'Acta2'],
        ['Kit', 'Srgn'],
        ['Stmn2', 'Stmn3', 'Hoxb8'],
        ['Fcer1g', 'C1qa'],
        ['Krt5', 'Krt14'],
        ['Cnmd', 'Col2a1'],
        ['Foxc1', 'Cldn11', 'Aldh1a2'],
        ['Cdh5', 'Kdr'],
        ['Sox10', 'Mpz'],
        ['Myog', 'Ttn', 'Actn2', 'Myh3'],
        ['Hba-a1', 'Hba-a2']
]

var_group_labels = [
                    'Fibroblast progenitors',
                    'Mitotic cells',
                    'Dermal Fibroblast',
                    'Neural progenitors', 
                    'Smooth Muscle', 
                    'Mast cells',
                    'Neural cells',
                    'Macrophages',
                    'Epidermis', 
                    'Cartilage',
                    'Meninges',
                    'Endothelial cells',
                    'Schwann cells',
                    'Skeletal muscle',
                    'Red Blood Cells'
                    ]

cluster_mapping = {
    "0": "Fibroblast progenitors",
    "1": "Fibroblast progenitors",
    "2": "Mitotic cells",
    "3": "Dermal fibroblasts",
    "4": "Neural progenitors",
    "5": "Smooth muscle",
    "6": "Fibroblast progenitors",
    "7": "Mast cells",
    "8": "Neural cells",
    "9": "Macrophages",
    "10": "Neural cells",
    "11": "Neural cells",
    "12": "Epidermis",
    "13": "Cartilage",
    "14": "Meninges",
    "15": "Neural cells",
    "16": "Neural cells",
    "17": "Endothelial cells",
    "18": "Schwann cells",
    "19": "Neural cells",
    "20": "Skeletal muscle",
    "21": "Neural cells",
    "22": "Red blood cells",
}

adata.obs["cell_type"] = adata.obs["leiden_0.4"].map(cluster_mapping)

# plot UMAP with annotations
sc.pl.umap(adata, color=['cell_type'], save='umap_annotations.png')
    
# plot for marker genes for cell type justification
sc.pl.dotplot(adata, 
              var_names=markers,
              groupby='cell_type',
              return_fig=False,
              color_map='Blues', 
              standard_scale='var',
              use_raw=False, 
              swap_axes=False, 
              save='cell_type_annotations_dotplot.png', 
              )

# heatmap of cell types across sample
adata.obs['day'] = adata.obs['sample'].replace({'JJ003': 'E13.5', 'JJ004': 'E12.5', 'JJ005': 'E11.5'})

sns.set(font_scale=0.85)
cmtx = sc.metrics.confusion_matrix("cell_type", "day", adata.obs)
cmtx = cmtx.iloc[:, [2,1,0]] # reorder for E11.5 to come first
sns_plot = sns.heatmap(cmtx, cmap='Blues')
sfig = sns_plot.get_figure()
sfig.tight_layout()
sfig.savefig('heatmap_cell_types_across_samples.png',  orientation="landscape")

# bat genes expressions
bat_genes = ['Ebf2', 'Sox9', 'Col2a1', 'Cdh4', 'Gata6', 'Cebpa', 'Pparg', 'Cidea']
sc.pl.umap(adata, color=bat_genes, color_map='Blues', save='ebf2_sox9_col2a1_exp.png')







######################################################################################
file = '/Users/rikac/Documents/project_lab/spring_25/scrna_seq_analysis/bc3306_mt_filter_sc_updated.h5ad'
adata = sc.read(file)

# figure 1C: UMAP of ebf2+/sox9+/col2a1- cells and UMAP of ebf2+/sox9+/col2a1+ cells

# ebf2+/sox9+/col2a1+ == 980 cells
coex = (adata.raw[:, '{}'.format('Ebf2')].X.todense() > 0) & (adata.raw[:, '{}'.format('Sox9')].X.todense() > 0) & (adata.raw[:, '{}'.format('Col2a1')].X.todense() > 0)
coex_list = [item for sublist in coex.tolist() for item in sublist]
adata.obs['ebf2+_sox9+_col2a1+'] = pd.Categorical(coex_list, categories=[True, False])
sc.pl.umap(adata, color='ebf2+_sox9+_col2a1+', groups=[True])

# ebf2+/sox9+/col2a1- == 1613 cells
coex1 = (adata.raw[:, '{}'.format('Ebf2')].X.todense() > 0) & (adata.raw[:, '{}'.format('Sox9')].X.todense() > 0) & (adata.raw[:, '{}'.format('Col2a1')].X.todense() == 0)
coex_list1 = [item for sublist in coex1.tolist() for item in sublist]
adata.obs['ebf2+_sox9+_col2a1-'] = pd.Categorical(coex_list1, categories=[True, False])
sc.pl.umap(adata, color='ebf2+_sox9+_col2a1-', groups=[True])

# the two figures on one page
sc.pl.umap(adata, color=['ebf2+_sox9+_col2a1+', 'ebf2+_sox9+_col2a1-'], groups=[True], legend_loc=None, frameon=False, 
           save='_ebf2_sox9_col2a1.png')

# plot the two conditions onto one figure with different colors
adata.obs['ebf2_sox9_col2a1'] = np.select(
    [coex_list, coex_list1],
    ['ebf2+_sox9+_col2a1+', 'ebf2+_sox9+_col2a1-'],
    default='other'
)
adata.obs['ebf2_sox9_col2a1'] = pd.Categorical(adata.obs['ebf2_sox9_col2a1'])
sc.pl.umap(adata, color='ebf2_sox9_col2a1', palette=['#E41A1C', '#377EB8', '#cacccc'], 
           frameon=False, save='_ebf2_sox9_col2a1_redblue.png')
sc.pl.umap(adata, color='ebf2_sox9_col2a1', palette=['#4DAF4A', '#E7298A', '#cacccc'], 
           frameon=False, save='_ebf2_sox9_col2a1_greenpink.png')

# figure 1D: chart of ebf2+/sox9+/col2a1- cells found in each cluster
cluster = adata.obs['leiden_0.4'].unique()
count_dict = {}
for i in cluster:
    counts = adata.obs[(adata.obs['ebf2+_sox9+_col2a1-'] == True) & (adata.obs['leiden_0.4'] == i)]
    count_dict[i] = len(counts)

df = pd.DataFrame(list(count_dict.items()), columns=['Cluster', 'Count'])
df['Cluster'] = pd.to_numeric(df['Cluster'], errors='coerce')
df = df.sort_values(by='Cluster').reset_index(drop=True)
plt.figure(figsize=(8,5))
sns.barplot(data=df, x='Cluster', y='Count')
plt.xlabel('Cluster')
plt.ylabel('Ebf2+/Sox9+/Col2a1- cells')
plt.show()

# figure 2A: umap FP subset (0,1,6)
# subset cluster 0,1,6
FP_mask = adata.obs['leiden_0.4'].isin(['0', '1', '6'])
FP_adata_subset = adata[FP_mask, :]     # shape: 7819 × 18973
# plot just the 0,1,6 umap
sc.pl.umap(FP_adata_subset)
# plot ebf2+_sox9+_col2a1- cells on FP_adata_subet
sc.pl.umap(FP_adata_subset, color='ebf2+_sox9+_col2a1-', groups=[True], frameon=False, 
           legend_loc=None, save='_ebf2+_sox9+_col2a1-_FPclusters.png')



# figure 2B: expression plots of BAT markers in ebf2+_sox9+_col2a1- cells within FP clusters
# 2B.1 = bat markers in FP clusters
bat_genes = ['Ebf2', 'Sox9', 'Cdh4', 'Pparg', 'Cidea', 'Hoxa5', 'Gata6', 'Cebpa']

with rc_context({"figure.figsize": (3, 3)}):
    sc.pl.umap(FP_adata_subset, color=bat_genes, color_map = 'plasma', s=50, frameon=False, ncols=5, vmax="p99")
    plt.savefig('_bat_markers_FP.png')

# 2B.2 = bat markers in ebf2+_sox9+_col2a1- cells within FP clusters
# this is going to have a different shape than the one above, since it is subsetting within FP
# for just ebf2+_sox9+_col2a1- cells
bat_genes1 = ['Ebf2', 'Sox9', 'Cdh4', 'Pparg', 'Cidea', 'Hoxa5', 'Gata6', 'Cebpa']    # not plotting Col2a1 bc 0
FP_subset_cond = FP_adata_subset[FP_adata_subset.obs['ebf2+_sox9+_col2a1-'] == True, :] # shape: 1255x18973
with rc_context({"figure.figsize": (3, 3)}):
    sc.pl.umap(FP_subset_cond, color=bat_genes1, color_map = 'plasma', s=50, frameon=False, vmax="p99")
    plt.savefig('_bat_markers_FPcondition.png')



# figure 2D: using fig 2B.2, label by sample of origin
# plot FP clusters and FP clusters+conditons and label by samples
FP_subset_cond.obs['day'] = FP_subset_cond.obs['sample'].replace({'JJ003': 'E13.5', 'JJ004': 'E12.5', 'JJ005': 'E11.5'})
FP_adata_subset.obs['day'] = FP_adata_subset.obs['sample'].replace({'JJ003': 'E13.5', 'JJ004': 'E12.5', 'JJ005': 'E11.5'})

sc.pl.umap(FP_adata_subset, color='day', palette=['#e6564c', '#1b9110', '#145ea3'], frameon=False,
           save='_samples_FP.png')
sc.pl.umap(FP_subset_cond, color='day', palette=['#e6564c', '#1b9110', '#145ea3'], frameon=False, 
           save='_samples_FPcondition.png')

# if I want to plot bat markers in FP+condition, need to convert bat markers expression to categorical
# meaning, e.g., if pparg is expressed or not, labeled/colored if Yes+day

days = FP_subset_cond.obs['day']

for i in bat_genes1:
    bat_day = (FP_subset_cond.raw[:, '{}'.format(i)].X.todense() > 0)
    bat_day_list = [item for sublist in bat_day.tolist() for item in sublist]
    FP_subset_cond.obs[i] = pd.Categorical(bat_day_list, categories=[True, False])

    FP_subset_cond.obs[i+'_day'] = np.where(
        bat_day_list,
        FP_subset_cond.obs['day'],
        'None'
    )

custom_palette = {
    'E11.5': '#145ea3',
    'E12.5': '#1b9110', 
    'E13.5': '#e6564c',  
    'None': '#d3d3d3'  
}
bat_genes_day= ['Ebf2_day', 'Sox9_day', 'Cdh4_day', 'Pparg_day', 'Cidea_day', 'Hoxa5_day', 'Gata6_day', 'Cebpa_day'] 
sc.pl.umap(FP_subset_cond, color=bat_genes_day, palette=custom_palette, frameon=False, 
           save='_bat_markers_FPconditions_days.png')

##################################################################
# reclustering 0,1,6 and output DEGs

res_list = ['leiden_0.2', 'leiden_0.4', 'leiden_0.6', 'leiden_0.8', 'leiden_1.0', 'leiden_1.2', 'leiden_1.4']

def plot_umap_leiden_markers(adata, resolutions_list):
    for resolution in resolutions_list:
        sc.tl.leiden(adata, resolution=float(resolution.split("_")[1]), key_added=resolution)

    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(adata, color=resolutions_list, frameon=False, legend_loc='on data', show=False)
        # plt.title(name)
        plt.savefig('_umap_colorbyleiden.png')

    for resolution in resolutions_list:
        res_key = "DGE" + resolution.replace("_", "")
        res_key_filename = res_key.replace(".", "")
        resolution_csv_name = resolution.replace("_", "").replace(".", "")

        # Find marker genes
        with rc_context({'figure.figsize': (5, 5)}):
            sc.tl.rank_genes_groups(adata, resolution, method='wilcoxon', use_raw=True, key_added=res_key)
            sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key=res_key, show=False)
            # plt.tight_layout()
            # plt.title(name)
            plt.savefig('_rank_genes_groups_' + res_key_filename + '.png')

        # Saving markers with pvalues to csv
        result = adata.uns[res_key]
        groups = result['names'].dtype.names
        degs_names = dict()
        degs_names[resolution_csv_name] = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals_adj']}).head(500)
        degs_names[resolution_csv_name].to_csv('_degs_names_' + resolution_csv_name + '.csv', index=False)

        degs = pd.read_csv('_degs_names_' + resolution_csv_name + '.csv')

        # Saving markers only to excel
        sel_col = []
        for colname in degs.columns.tolist():
            if "_n" in colname:
                sel_col.append(colname)

        sel = degs[sel_col]
        sel.columns = range(sel.columns.size)

        sel.to_excel('_degs_names_' + resolution_csv_name + '.xlsx',
                     sheet_name="Sheet1", startrow=1, header=True, index=False)

plot_umap_leiden_markers(FP_adata_subset, res_list)

# making the same set of figures using reclustered data (leiden_0.4)

# figure 1C: UMAP of ebf2+/sox9+/col2a1- cells and UMAP of ebf2+/sox9+/col2a1+ cells
sc.pl.umap(FP_adata_subset, color=['ebf2+_sox9+_col2a1+', 'ebf2+_sox9+_col2a1-'], groups=[True], legend_loc=None, frameon=False)

# pparg is significant in cluster 0 of leiden_0.4, Pparg - 1.86694757242082E-09 
