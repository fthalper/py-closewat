#!/usr/bin/env python3
import requests
import csv
import re
import sys

# RCSB PDB search and entry APIs
def get_pdb_ids(min_res, max_res):
    url = 'https://search.rcsb.org/rcsbsearch/v2/query?json='
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "greater", "value": min_res}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less", "value": max_res}}
            ]
        },
        "return_type": "entry"
    }
    r = requests.post(url, json=query)
    r.raise_for_status()
    data = r.json()
    return [item['identifier'] for item in data.get('result_set', [])]


def get_entry_info(pdb_id):
    url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'
    r = requests.get(url)
    if r.status_code != 200:
        return None, None
    j = r.json()
    info = j.get('rcsb_entry_info', {})
    # resolution_combined is a list
    res = info.get('resolution_combined', [None])[0]
    # r_free or r_factor_R_free
    rfree = info.get('r_free', [None])[0] or info.get('r_free', [None])[0]
    return res, rfree


def count_water(pdb_id):
    pdb_url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    r = requests.get(pdb_url)
    if r.status_code != 200:
        return 0
    text = r.text
    m = re.search(r"^FORMUL\s+\d+\s+HOH\s+\*(\d+)\(H2 O\)", text, re.MULTILINE)
    if m:
        return int(m.group(1))
    # fallback: count HOH lines
    return len([line for line in text.splitlines() if line.startswith('HETATM') and ' HOH ' in line])


def main():
    if len(sys.argv) >= 3:
        try:
            min_res = float(sys.argv[1])
            max_res = float(sys.argv[2])
        except:
            print("Usage: fetch_and_count_pdb_info.py <min_resolution> <max_resolution>")
            sys.exit(1)
    else:
        min_res, max_res = 1.01, 1.06

    print(f"Querying PDB for X-ray structures with {min_res}Å < resolution < {max_res}Å...")
    pdb_ids = get_pdb_ids(min_res, max_res)
    print(f"Found {len(pdb_ids)} entries")

    out_file = 'pdb_summary.csv'
    with open(out_file, 'w', newline='') as csvf:
        writer = csv.writer(csvf)
        writer.writerow(['pdb_id', 'resolution', 'r_free', 'water_count'])
        # process a few (first 20) for initial test
        for pdb_id in pdb_ids[:20]:
            pdb = pdb_id.upper()
            res, rfree = get_entry_info(pdb)
            wcount = count_water(pdb)
            print(f"{pdb}: res={res}, Rfree={rfree}, waters={wcount}")
            writer.writerow([pdb, res, rfree, wcount])
    print(f"Summary written to {out_file}")

if __name__ == '__main__':
    main()
