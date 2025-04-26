import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load the CSV data
df = pd.read_csv('pdb_water_analysis.csv')

# Clean up the data - remove rows with 0 water molecules (likely errors or missing data)
df_clean = df[df['water_count'] > 0].copy()

# Basic statistics
total_structures = len(df)
structures_with_water = len(df_clean)
avg_water = df_clean['water_count'].mean()
median_water = df_clean['water_count'].median()
min_water = df_clean['water_count'].min()
max_water = df_clean['water_count'].max()
std_water = df_clean['water_count'].std()

# Print statistics
print(f"Total structures analyzed: {total_structures}")
print(f"Structures with water molecules: {structures_with_water}")
print(f"Structures without water molecules: {total_structures - structures_with_water}")
print(f"Average water molecules per structure: {avg_water:.2f}")
print(f"Median water molecules per structure: {median_water:.2f}")
print(f"Minimum water molecules: {min_water}")
print(f"Maximum water molecules: {max_water}")
print(f"Standard deviation: {std_water:.2f}")

# Create a histogram of water molecule counts
plt.figure(figsize=(12, 6))
plt.hist(df_clean['water_count'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
plt.title('Distribution of Water Molecules in PDB Structures (1.01Å < Resolution < 1.06Å)', fontsize=14)
plt.xlabel('Number of Water Molecules', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.grid(True, alpha=0.3)
plt.savefig('water_distribution_histogram.png', dpi=300, bbox_inches='tight')

# Create a scatter plot of resolution vs water count
plt.figure(figsize=(12, 6))
plt.scatter(df_clean['resolution'], df_clean['water_count'], alpha=0.5, color='blue')
plt.title('Resolution vs Water Molecule Count', fontsize=14)
plt.xlabel('Resolution (Å)', fontsize=12)
plt.ylabel('Number of Water Molecules', fontsize=12)
plt.grid(True, alpha=0.3)
plt.savefig('resolution_vs_water.png', dpi=300, bbox_inches='tight')

# If Rfree values are available, create a scatter plot of Rfree vs water count
if df_clean['rfree'].notna().any():
    df_with_rfree = df_clean[df_clean['rfree'].notna()]
    plt.figure(figsize=(12, 6))
    plt.scatter(df_with_rfree['rfree'], df_with_rfree['water_count'], alpha=0.5, color='green')
    plt.title('Rfree vs Water Molecule Count', fontsize=14)
    plt.xlabel('Rfree Value', fontsize=12)
    plt.ylabel('Number of Water Molecules', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.savefig('rfree_vs_water.png', dpi=300, bbox_inches='tight')

# Create bins for water molecule counts and calculate statistics for each bin
bins = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, float('inf')]
bin_labels = ['1-100', '101-200', '201-300', '301-400', '401-500', '501-600', 
              '601-700', '701-800', '801-900', '901-1000', '1000+']

df_clean['water_bin'] = pd.cut(df_clean['water_count'], bins=bins, labels=bin_labels)
bin_counts = df_clean['water_bin'].value_counts().sort_index()

print("\nDistribution of structures by water molecule count:")
for bin_name, count in bin_counts.items():
    percentage = (count / structures_with_water) * 100
    print(f"{bin_name}: {count} structures ({percentage:.2f}%)")

# Save the cleaned data with bins to a new CSV
df_clean.to_csv('water_analysis_with_bins.csv', index=False)

print("\nAnalysis complete. Visualizations saved as PNG files.")
