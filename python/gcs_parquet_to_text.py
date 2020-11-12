#!/usr/bin/env python
#
# Jeremy Schwartzentruber
#
# Gets the first few lines of a parquet file from Google cloud, and saves it
# in text format.
#
import sys
import os
import pandas as pd
import subprocess as sp
import argparse
import re
import gzip

#args.file = "gs://genetics-portal-staging/v2d/200207/toploci.parquet"

def main():
    args = parse_args()

    fname = args.file
    if re.match('^gs://', args.file):
        # First download the file with gsutil
        fname = "tmp." + os.path.basename(args.file)
        if not os.path.isfile(fname) or args.overwrite:
            cmd = "gsutil cp {} {}".format(args.file, fname)
            sp.call(cmd, shell=True)
    
    df = pd.read_parquet(fname, engine='pyarrow')
    if args.nlines is not None:
        nlines = min(args.nlines, df.shape[0])
        df = df[:nlines]
    
    if args.out_file is not None:
        if args.pretty:
            if re.search('.gz$', args.out_file):
                with gzip.open(args.out_file, "wt") as f:
                    f.write(df.to_string(index=None, na_rep='NA') + "\n")
            else:
                with open(args.out_file, "wt") as f:
                    f.write(df.to_string(index=None, na_rep='NA') + "\n")
        else:
            df.to_csv(args.out_file, sep='\t', index=None, na_rep='NA')
    else:
        if args.pretty:
            print(df.to_string(index=None, na_rep='NA'))
        else:
            df.to_csv(sys.stdout, sep='\t', index=None, na_rep='NA', compression='gzip')


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', metavar="<file>",
                        help=('A parquet file (may be on Google Cloud Storage)'), type=str, required=True)
    parser.add_argument('--nlines', metavar="<int>",
                        help=('Number of lines to get'), type=int, required=False)
    parser.add_argument('--out_file', metavar="<file>",
                        help=("Output file path (created if not present)"), type=str, required=False)
    parser.add_argument('--pretty',
                        help=("If True, then for a Google storage file the file will be re-downloaded even if the local copy exists"), action='store_true')
    parser.add_argument('--overwrite',
                        help=("If True, then for a Google storage file the file will be re-downloaded even if the local copy exists"), action='store_true')
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    main()
