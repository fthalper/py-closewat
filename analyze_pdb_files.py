import os
import re
import csv
from pathlib import Path

# Directory containing PDB files
pdb_dir = Path("pdb_files")
output_file = "pdb_water_analysis.csv"

# Regular expressions to extract information
resolution_pattern = re.compile(r"REMARK\s+2\s+RESOLUTION\.\s+(\d+\.\d+)\s+ANGSTROMS")
rfree_pattern = re.compile(r"REMARK\s+3\s+FREE R VALUE\s+:\s+(\d+\.\d+)")
water_pattern = re.compile(r"FORMUL\s+\d+\s+HOH\s+\*(\d+)\(H2 O\)")

# Function to extract information from a PDB file
def analyze_pdb_file(file_path):
    pdb_id = os.path.basename(file_path).replace("pdb", "").replace(".ent", "")
    resolution = None
    rfree = None
    water_count = 0
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
            # Extract resolution
            res_match = resolution_pattern.search(content)
            if res_match:
                resolution = float(res_match.group(1))
            
            # Extract Rfree
            rfree_match = rfree_pattern.search(content)
            if rfree_match:
                rfree = float(rfree_match.group(1))
            
            # Extract water count
            water_match = water_pattern.search(content)
            if water_match:
                water_count = int(water_match.group(1))
            else:
                # If no water formula record, count HOH records
                water_count = content.count("HOH")
                
        return {
            "pdb_id": pdb_id,
            "resolution": resolution,
            "rfree": rfree,
            "water_count": water_count
        }
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return {
            "pdb_id": pdb_id,
            "resolution": None,
            "rfree": None,
            "water_count": 0
        }

# Main function
def main():
    results = []
    
    # Get all PDB files
    pdb_files = list(pdb_dir.glob("*.ent"))
    total_files = len(pdb_files)
    print(f"Found {total_files} PDB files")
    
    # Process each file
    for i, file_path in enumerate(pdb_files):
        if i % 50 == 0:
            print(f"Processing file {i+1}/{total_files}...")
        
        data = analyze_pdb_file(file_path)
        results.append(data)
        

    # Filter results for resolution between 1.01 and 1.06 Å
    results = [r for r in results if r['resolution'] is not None and 1.01 < r['resolution'] < 1.06]
    print(f"Filtered structures count: {len(results)} (1.01Å < resolution < 1.06Å)")

    # Sort results by PDB ID
    results.sort(key=lambda x: x["pdb_id"])
    
    # Write results to CSV
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['pdb_id', 'resolution', 'rfree', 'water_count']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for data in results:
            writer.writerow({
                'pdb_id': data['pdb_id'],
                'resolution': data['resolution'],
                'rfree': data['rfree'],
                'water_count': data['water_count']
            })
    
    # Calculate some statistics
    water_counts = [data['water_count'] for data in results if data['water_count'] > 0]
    if water_counts:
        avg_water = sum(water_counts) / len(water_counts)
        min_water = min(water_counts)
        max_water = max(water_counts)
        print(f"\nStatistics for water molecules in processed structures:")
        print(f"Average water molecules per structure: {avg_water:.2f}")
        print(f"Minimum water molecules: {min_water}")
        print(f"Maximum water molecules: {max_water}")
    
    print(f"\nResults written to {output_file}")

if __name__ == "__main__":
    main()
