import pandas as pd

df = pd.read_csv(
    r"C:\Users\97250\Downloads\gene_tpm_2017-06-05_v8_artery_coronary.gct\gene_tpm_artery_coronary.gct",
    sep="\t",
    skiprows=2
)

print("✅ File loaded successfully!")
print(df.shape)
print(df.head())
# 1) Genes catalogue – two columns only
genes_catalogue = df[['Name', 'Description']].rename(
    columns={'Name': 'gene_id', 'Description': 'gene_name'}
)
print("\n--- Genes catalogue ---")
print(genes_catalogue.head())

# 2) Gene-expression matrix – rows = gene IDs, columns = samples
expression_data = df.drop(columns=['Description'])
expression_data.set_index('Name', inplace=True)
print("\n--- Expression data (first 5 genes × 5 samples) ---")
print(expression_data.iloc[:5, :5])

# Optional: save them for later labs
genes_catalogue.to_csv(r"C:\Users\97250\Downloads\genes_catalogue.csv", index=False)
expression_data.to_csv(r"C:\Users\97250\Downloads\gene_expression_data.csv")
print("\n✅ Saved both CSV files to Downloads.")
