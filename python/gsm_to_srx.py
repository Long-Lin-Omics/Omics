#!/usr/bin/env python3

from pysradb import SRAweb
import sys

if len(sys.argv) != 3:
    print("Usage: python gsm_to_srx.py gsm_ids.txt srx_ids.txt")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

db = SRAweb()

with open(input_file, 'r') as f:
    gsm_ids = [line.strip().split()[0] for line in f if line.strip()]

with open(output_file, 'w') as out:
    for gsm in gsm_ids:
        print(f"Fetching SRX for {gsm}...")
        try:
            df = db.gsm_to_srx(gsm)
            if not df.empty:
                for srx in df['experiment_accession']:
                    out.write(srx + '\n')
            else:
                print(f"  No SRX found for {gsm}.")
        except Exception as e:
            print(f"  Error processing {gsm}: {e}")

print("All SRX IDs written to", output_file)
