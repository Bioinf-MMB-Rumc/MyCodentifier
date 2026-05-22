#!/usr/bin/env python

######
INFO = "convert_resdb_v1.py"
__version__ = "0.1"
######

"""
Title:          convert_resdb_v1.py
Author:         P.D. Koopman
Date:           25-09-2023 (dd-mm-yyyy)
Description:    Convert the WHO resistance database of TB to be suitable for mycodentifier
"""

import pandas as pd
import argparse
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
    parser.add_argument("-i", "--whodb", required=True,
                        help="Path to the excel file containing the WHO database")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Dir to store the new tb resistance database")

    # Parse the arguments
    arguments = parser.parse_args()
    return arguments


def parse_who_db(who_excel):
    """
    This will parse the input WHO TB resistance database
    Input:
    who_excel - pandas Excel file, dataframe cont
    Output:
    mut_df - pandas df, contains all the statistics of effects of the mutations
    ind_df - pandas df, contains all the information of the mutations themselves
    """
    mut_df = pd.read_excel(who_excel, "Mutation_catalogue")
    ind_df = pd.read_excel(who_excel, "Genome_indices")
    return mut_df, ind_df


def combine_res_information(mut_df, ind_df):
    """
    This will combine the relevant information from the WHO database for the Mycodentifier TB resistance database
    """
    res_mut = mut_df.loc[mut_df["FINAL CONFIDENCE GRADING"].isin(["1) Assoc w R", "2) Assoc w R - Interim"])]
    res_mut["variant"] = res_mut["variant (common_name)"].copy().str.split(" ").str[0]
    res_mut_hgvs = pd.merge(res_mut[["drug", "variant"]],
                            ind_df[["final_annotation.Gene", "variant", "final_annotation.TentativeHGVSNucleotidicAnnotation",
                                    "final_annotation.TentativeHGVSProteicAnnotation"]],
                            left_on="variant", right_on="variant")
    res_mut_hgvs.rename(columns={"final_annotation.TentativeHGVSNucleotidicAnnotation": "hgvs_nuc",
                                 "final_annotation.TentativeHGVSProteicAnnotation": "hgvs_prot"}, inplace=True)
    res_mut_hgvs["lit"] = ""
    return res_mut_hgvs[["final_annotation.Gene", "hgvs_nuc", "drug", "lit"]]


if __name__ == "__main__":
    args = parse_args()

    # Read and parse the who database
    who_file = pd.ExcelFile(args.whodb)
    mutation_df, index_df = parse_who_db(who_file)

    # Combine the relevant information for the Mycodentifier TB resistance database
    res_db = combine_res_information(mutation_df, index_df)
    date = datetime.datetime.now().strftime('%Y%m%d')
    res_db.to_csv(f"{args.outdir}/{date}_NC_000962.3_resistance_db.tsv", sep="\t", header=False, index=False)
