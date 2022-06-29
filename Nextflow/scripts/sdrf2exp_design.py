#!/usr/bin/env python3

import pandas as pd

sdrf = pd.read_csv("sdrf_local.tsv", sep='\t')
sdrf = sdrf.astype(str)
sdrf.columns = map(str.lower, sdrf.columns)  # convert column names to lower-case
changing_columns = list()
for col in sdrf.columns:
    if "factor" in col or "characteristics" in col:
        if sdrf[col].value_counts().size > 1:
            changing_columns.append(col)


changing_columns = [x for x in changing_columns if "factor" in x  ]

sdrf_out = pd.DataFrame()
sdrf_out["raw_file"] = sdrf["comment[file uri]"]
if (len(changing_columns) > 0):
	sdrf_out["exp_conditions"] = sdrf[changing_columns].agg('_'.join, axis=1)
else: 
	sdrf_out["exp_conditions"] = "A"

sdrf_out.to_csv("exp_design.txt", sep="\t", index=False)
    

