#!/usr/bin/env python

######
INFO = "check_zero_depth.py"
__version__ = "0.1"
######

"""
Title:          check_zero_depth.py
Author:         P.D. Koopman
Date:           11-06-2025 (dd-mm-yyyy)
Description:    Check for non-sequenced positions in the resistance database
"""

import pandas as pd
import os
import argparse
import simplejson as json
import datetime
import numpy as np


def parse_args():
    """
    Parse commandline arguments
    Input:
    Output:
    args - argparse object, contains arguments provided by user
    """
    parser = argparse.ArgumentParser(description=INFO,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-s", "--samplename", required=True,
                        help="Name of the sample being processed")
    parser.add_argument("-c", "--coverage", required=True,
                        help="Path to the coverage file containing sample depths per genome position")
    parser.add_argument("-r", "--resdb", required=True,
                        help="Path to resistance database including genome position for mutations")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Dir to store zero depth resistance mutations in", default=os.path.abspath("./"))
    parser.add_argument("-j", "--json", required=False, type=argparse.FileType('r'),
                        help="Path to meta json for tracking all functions in MyCodentifier pipeline")

    # Parse the arguments
    arguments = parser.parse_args()
    return arguments


def parse_input(cov, res_db):
    """
    This parse the input csv files.
    Input:
    cov - string, path to sample coverage file
    res_db - string, path to resistance db with genomic positions
    Output:
    cov_df - Pandas dataframe
    res_df - Pandas dataframe
    """

    cov_df = pd.read_csv(cov, sep="\t", names=["chromosome", "position", "depth"])
    res_df = pd.read_csv(res_db, sep="\t", header=0)

    return cov_df, res_df


def combine_depth_res(cov_df, res_df):
    """
    This will combine the sample depth and the genomic positions of resistance mutations.
    """
    mut_df = pd.merge(res_df, cov_df[["position", "depth"]], on='position', how='left')
    mut_df["depth"].fillna(value=0, inplace=True)

    return mut_df, mut_df.loc[mut_df["depth"] == 0].drop_duplicates(subset=['drug', 'position'])


def load_json(f):
    try:
        dict_ob = json.load(f)
        return dict_ob
    except json.JSONDecodeError:
        raise Exception('Provide a json file as run info file')


def write_results(meta_json, zero_json):
    """
    write the generated json with mutation info to pipeline json
    :param meta_json: meta json of mycodentifier pipeline
    :param zero_json: json with zero coverage positions
    """
    meta = load_json(meta_json)
    meta['zero_depth'] = zero_json
    print(meta)
    with open(meta_json.name, 'w') as dbout:
        json.dump(meta, dbout)


if __name__ == "__main__":
    args = parse_args()
    coverage_df, resdb_df = parse_input(args.coverage, args.resdb)

    # Combine coverage information and mutation genome position
    mutation_depth_df, zero_depth_df = combine_depth_res(coverage_df, resdb_df)
    zero_depth_df.to_csv(f"{args.outdir}/{args.samplename}_zerodepth.tsv", sep="\t", header=True, index=False)
    mutation_depth_df.to_csv(f"{args.outdir}/{args.samplename}.tsv", sep="\t", header=True, index=False)

    # Write to meta json
    if args.json:
        write_results(args.json, zero_depth_df.to_dict(orient="index"))
