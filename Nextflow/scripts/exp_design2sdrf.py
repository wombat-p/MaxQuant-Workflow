#!/usr/bin/env python3

import pandas as pd


expd = pd.read_csv("exp_design_in.txt", sep='\t')
expd.columns = map(str.lower, expd.columns)  # convert column names to lower-case
sdrf_out = expd
sdrf_out.columns = ['comment[file uri]', 'characteristics[experimental samples]']
sdrf_out['comment[data file]'] = expd['comment[file uri]']
sdrf_out['comment[instrument]'] = "not available"
sdrf_out['comment[label]'] = 'AC=MS:1002038;NT=label free sample'
sdrf_out["comment[cleavage agent details]"] = 'NT=trypsin/P;AC=MS:1001313'
sdrf_out.insert(0, 'source name', "Sample")


sdrf_out.to_csv("sdrf_local.tsv", sep="\t", index=False)
    

