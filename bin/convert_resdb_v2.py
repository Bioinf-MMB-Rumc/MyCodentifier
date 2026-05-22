#!/usr/bin/env python

######
INFO = "convert_resdb_v2.py"
__version__ = "0.1"
######

"""
Title:          convert_resdb_v2.py
Author:         P.D. Koopman
Date:           13-12-2023 (dd-mm-yyyy)
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
    mut_df = pd.read_excel(who_excel, "Catalogue_master_file", skiprows=[0, 1])
    genomic_df = pd.read_excel(who_excel, "Genomic_coordinates")
    return mut_df, genomic_df


def combine_res_information(mut_df):
    """
    This will combine the relevant information from the WHO database for the Mycodentifier TB resistance database
    """
    res_mut = mut_df.loc[mut_df["FINAL CONFIDENCE GRADING"].isin(["1) Assoc w R", "2) Assoc w R - Interim"])]
    res_mut.drop(res_mut[res_mut["mutation"] == "LoF"].index, inplace=True)
    res_mut.drop(res_mut[res_mut["mutation"] == "deletion"].index, inplace=True)
    res_mut.drop(res_mut[res_mut["mutation"]
                 .str.contains(r"^[pn]\.-?[A-Za-z]{3}(\d){1,6}\_[A-Za-z]{3}(\d){1,6}(del|ins).*")].index, inplace=True)
    res_mut.drop(res_mut[res_mut["mutation"].str.contains("\?")].index, inplace=True)

    # Specifically implemented for FabG1 and inhA, check if alias is present
    res_mut[["Alias", "Commentcore", "Commentrest"]] = res_mut["Comment"].str.split(" ", n=2, expand=True)
    res_mut[["Alias_gene", "Alias_mut"]] = res_mut["Commentcore"].str.split("_", n=1, expand=True)
    res_mut["Alias_mut"] = res_mut["Alias_mut"].str.rstrip(".")
    res_mut.loc[res_mut["Alias"] == "Alias", "gene"] = res_mut.loc[
        res_mut["Alias"] == "Alias", "Alias_gene"]
    res_mut.loc[res_mut["Alias"] == "Alias", "mutation"] = res_mut.loc[
        res_mut["Alias"] == "Alias", "Alias_mut"]

    res_mut["lit"] = ""
    return res_mut[["gene", "mutation", "drug", "lit"]], res_mut


def attach_genome_coordinates(resistance_db, genome_df):
    """
    This will attach genome coordinates of mutations to the resistance database.
    """
    genome_res_db = pd.merge(resistance_db, genome_df, on='variant', how='left')

    return genome_res_db[["gene", "mutation", "drug", "lit", "position"]]


if __name__ == "__main__":
    args = parse_args()

    # Read and parse the who database
    who_file = pd.ExcelFile(args.whodb)
    mutation_df, genomic_df = parse_who_db(who_file)

    # Combine the relevant information for the Mycodentifier TB resistance database
    res_db, res_db_full = combine_res_information(mutation_df)

    # Build a separate file with the genomic coordinates of the resistance mutations
    genomic_res_db = attach_genome_coordinates(res_db_full, genomic_df)

    date = datetime.datetime.now().strftime('%Y%m%d')
    res_db.to_csv(f"{args.outdir}/{date}_NC_000962.3_resistance_db.tsv", sep="\t", header=False, index=False)
    genomic_res_db.to_csv(f"{args.outdir}/{date}_NC_000962.3_resistance_db_genomic.tsv", sep="\t", header=True, index=False)
