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
# Example:
# gcs_parquet_to_text.py -f gs://genetics-portal-staging/v2d/200207/toploci.parquet
def main():
    args = parse_args()
    if args.mirror_path and args.out_file is not None:
        print("Arguments --out_file and --mirror_path are mutually exclusive")
        exit(1)

    fname = args.file
    downloaded_file = None
    if re.match('^gs://', args.file):
        # First download the file with gsutil
        fname = "tmp." + os.path.basename(args.file)
        if not os.path.isfile(fname) or args.overwrite:
            cmd = "gsutil cp {} {}".format(args.file, fname)
            result = sp.call(cmd, shell=True)
            if result != 0:
                print("gsutil call failed", file=sys.stderr)
                exit(1)
            downloaded_file = fname
    
    try:
        df = pd.read_parquet(fname, engine='pyarrow')
        shape_msg = ''
        if args.nlines is not None:
            nlines = min(args.nlines, df.shape[0])
            if nlines < df.shape[0]:
                shape_msg = '...{} total rows\n'.format(df.shape[0])
            df = df[:nlines]
        
        out_file = args.out_file
        if args.mirror_path:
            out_file = args.file.replace("gs://", "")
            out_file = os.path.splitext(out_file)[0] + '.tsv'
        if out_file is not None:
            # If out_file is empty, we create a relative path that is equivalent to
            # that in the input file
            dir = os.path.dirname(out_file)
            if dir and not os.path.exists(dir):
                os.makedirs(dir)
            if args.pretty:
                if re.search('.gz$', out_file):
                    with gzip.open(out_file, "wt") as f:
                        f.write(df.to_string(index=None, na_rep='NA') + "\n" + shape_msg)
                else:
                    with open(out_file, "wt") as f:
                        f.write(df.to_string(index=None, na_rep='NA') + "\n" + shape_msg)
            else:
                df.to_csv(out_file, sep='\t', index=None, na_rep='NA')
                if not re.search('.gz$', out_file):
                    with open(out_file, "at") as f:
                        f.write(shape_msg)
        else:
            if args.pretty:
                print(df.to_string(index=None, na_rep='NA') + shape_msg)
            else:
                df.to_csv(sys.stdout, sep='\t', index=None, na_rep='NA')
                print(shape_msg)
    finally:
        if downloaded_file is not None and not args.keep:
            os.remove(downloaded_file)
        


def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', metavar="<file>", type=str, required=True,
                        help=('A parquet file (may be on Google Cloud Storage)'))
    parser.add_argument('-n', '--nlines', metavar="<int>", type=int, required=False,
                        help=('Number of lines to get'))
    parser.add_argument('-o', '--out_file', metavar="<file>", type=str, required=False,
                        help=("Output file path (created if not present)"))
    parser.add_argument('-m', '--mirror_path', action='store_true',
                        help=("Save the output to a file with a relative filepath equivalent to the downloaded file."))
    parser.add_argument('--pretty', action='store_true',
                        help=("If True, then for a Google storage file the file will be re-downloaded even if the local copy exists"))
    parser.add_argument('--keep', action='store_true',
                        help=("If True, then keep the downloaded parquet file"))
    parser.add_argument('--overwrite', action='store_true',
                        help=("If True, then for a Google storage file the file will be re-downloaded even if the local copy exists"))
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    main()
