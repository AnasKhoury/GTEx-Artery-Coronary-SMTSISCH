import pandas as pd
import matplotlib.pyplot as plt

# Load the filtered metadata file
meta_path = r"C:\Users\97250\Downloads\GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
meta = pd.read_csv(meta_path, sep="\t")

# Filter for good RNA integrity (RIN ≥ 6)
filtered_meta = meta[meta["SMRIN"] >= 6]

# Plot histogram
plt.figure(figsize=(7, 5))
plt.hist(filtered_meta["SMRIN"], bins=20, color="skyblue", edgecolor="black")
plt.title("Distribution of RNA Integrity Number (RIN) After Filtering", fontsize=13)
plt.xlabel("RIN Value")
plt.ylabel("Number of Samples")
plt.grid(axis='y', alpha=0.6)
plt.tight_layout()

# Save figure to file
plt.savefig(r"C:\Users\97250\Downloads\RIN_After_Filtering.png", dpi=150)
plt.show()
