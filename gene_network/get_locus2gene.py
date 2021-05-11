# This script reads data from OT Genetics datasets with locus2gene scores,
# and saves a subset with the relevant studies.
from pyspark.sql import SparkSession
import pandas as pd

spark = SparkSession.builder.getOrCreate()
df = spark.read.load("/Users/jeremys/Downloads/l2g")
#df.head(2)

# The below is inefficient
#df = pd.read_parquet('/Users/jeremys/Downloads/l2g/')

# Katie de Lange's IBD study
dfs = df.filter(df["study_id"] == "GCST004131")

#print((df.count(), len(df.columns)))
#dfs.write.csv('/Users/jeremys/Downloads/l2g.IBD.deLange.tsv', sep='\t')
dfs.toPandas().to_csv('/Users/jeremys/Downloads/l2g.IBD.deLange.tsv', sep='\t', index=False)

# Jimmy Liu's IBD study
dfs = df.filter(df["study_id"] == "GCST003043")
#print((dfs.count(), len(dfs.columns)))
dfs.toPandas().to_csv('/Users/jeremys/Downloads/l2g.IBD.Liu.tsv', sep='\t', index=False)

