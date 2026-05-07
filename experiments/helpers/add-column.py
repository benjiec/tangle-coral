# add a column to a TSV

import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("tsv_file")
ap.add_argument("column")
ap.add_argument("value")
ap.add_argument("--output", "-o", help="Output file path (defaults to input)")
args = ap.parse_args()

df = pd.read_csv(args.tsv_file, sep='\t')
df[args.column] = args.value
        
output_path = args.output if args.output else args.tsv_file
df.to_csv(output_path, sep='\t', index=False)
