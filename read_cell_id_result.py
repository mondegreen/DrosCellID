import sys

import pandas as pd

from pyranges.pyranges import PyRanges

from io import StringIO

import pyranges as pr
from pyranges.version import __version__

def read_cell_id(f, w=500, as_df=False, nrows=None):

  """Return cell id ressult file as PyRanges.

  This is a reader for files that follow the bed format. They can have from
  5 or 14 columns which will be named like so:

  Chromosome1 Pos1 Num_Reads1 UniJunctionsLeft1 UniJunctionsRight1 sep skip skip sep
    Chromosome2 Pos2 Num_Reads2 UniJunctionsLeft2 UniJunctionsRight2

  Parameters
  ----------
  f : str

      Path to result file

  w : int, default 500

      Window out from chromosomal position to use

  as_df : bool, default False

      Whether to return as pandas DataFrame instead of PyRanges.

  nrows : int, default None

      Number of rows to return.
  """

  columns = "Chromosome1 Pos1 Num_Reads1 UniJunctionsLeft1 UniJunctionsRight1 sep1 skip1 skip2 sep2 Chromosome2 Pos2 Num_Reads2 UniJunctionsLeft2 UniJunctionsRight2".split(
  )

  if f.endswith(".gz"):
    import gzip
    fh = gzip.open(f)
  else:
    fh = open(f)

  lines = fh.read().splitlines()
  elements = []
  for line in lines:
    values = line.split("\t")
    if len(values) > 5 and values[5] == "x":
      entry1 = []
      if int(values[3]) > int(values[4]):
        entry1.append(values[0])
        b = int(values[1]) - w
        if b < 0:
          b = 0
        entry1.append(b)
        entry1.append(int(values[1])+1)
      else:
        entry1.append(values[0])
        entry1.append(int(values[1]))
        entry1.append(int(values[1])+w)
      elements.append(entry1)

      entry2 = []
      if int(values[12]) > int(values[13]):
        entry2.append(values[9])
        b = int(values[10]) - w
        if b < 0:
          b = 0
        entry2.append(b)
        entry2.append(int(values[10])+1)
      else:
        entry2.append(values[9])
        entry2.append(int(values[10]))
        entry2.append(int(values[10])+w)
      elements.append(entry2)
    else:
      entry = []
      if int(values[3]) > int(values[4]):
        entry.append(values[0])
        b = int(values[1]) - w
        if b < 0:
          b = 0
        entry.append(b)
        entry.append(int(values[1])+1)
      else:
        entry.append(values[0])
        entry.append(int(values[1]))
        entry.append(int(values[1])+w)
      elements.append(entry)

  df = pd.DataFrame(elements,columns=['Chromosome','Start','End'])
  if not as_df:
      return PyRanges(df)
  else:
      return df
