import numpy as np
import matplotlib.pyplot as plt

# This code assumes that you have a pandas dataframe of your differentially expressed genes
# The df should contain at least GENEID (gene names), log2FoldChange, and padj
# Note: indexing GENEID is case-sensitive

# Calculate negative log10 p-value and log2 fold change
df['neg_log10_pval'] = -np.log10(df['padj'])
df['log2_fold_change'] = df['log2FoldChange']

# Define colors for depleted, enriched, and non-significant genes
depleted_color = 'blue'
enriched_color = 'red'
non_sig_color = 'gray'

# Define cutoff for significant genes
pval_cutoff = 0.05

# Create boolean arrays for depleted, enriched, and non-significant genes
is_depleted = (df['log2_fold_change'] < 0) & (df['padj'] < pval_cutoff)
is_enriched = (df['log2_fold_change'] > 0) & (df['padj'] < pval_cutoff)
is_non_sig = (df['padj'] >= pval_cutoff)

# Create volcano plot
plt.scatter(df['log2_fold_change'][is_depleted], df['neg_log10_pval'][is_depleted], color=depleted_color, label='Depleted')
plt.scatter(df['log2_fold_change'][is_enriched], df['neg_log10_pval'][is_enriched], color=enriched_color, label='Enriched')
plt.scatter(df['log2_fold_change'][is_non_sig], df['neg_log10_pval'][is_non_sig], color=non_sig_color, label='Non-significant')

# Add axis labels and legend
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 p-value')
plt.legend()

# Set y-axis limit
plt.ylim(ymin, ymax)
plt.xlim(xmin, xmax)

def get_gene_data(gene_id, df):
    """
    Returns x and y values for a given gene ID.
    Assumes that the input DataFrame has 'GENEID', 'log2FoldChange', and 'padj' columns.
    """
    gene_data = df.loc[df['GENEID'] == gene_id].iloc[0].to_dict()
    x = gene_data['log2FoldChange']
    y = -np.log10(gene_data['padj'])
    return x, y

def create_gene_labels(gene_ids, df):
    """
    Creates a dictionary of gene labels for the given list of gene IDs and DataFrame.
    Assumes that the input DataFrame has 'GENEID', 'log2FoldChange', and 'padj' columns.
    """
    gene_labels = {}
    for gene_id in gene_ids:
        x, y = get_gene_data(gene_id, df)
        gene_labels[gene_id.upper()] = (x, y)  # use uppercase gene ID as the key
    return gene_labels

# Example usage:
genes_of_interest = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6']
gene_labels = create_gene_labels(genes_of_interest, df)

# Add gene labels to plot
for gene, (log2fc, neglog10p) in gene_labels.items():
    plt.text(log2fc, neglog10p, gene, ha='center', va='center', fontsize=8, fontweight='bold')

