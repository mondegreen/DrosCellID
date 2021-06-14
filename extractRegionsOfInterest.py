import os
import sys
import argparse
import read_cell_id_result as cid
from Bio import SeqIO

# read a game directory, modify the files therein, update crc16 checks
def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

# read parameters
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--list', action="store", dest="list", help="file with list of inputs")
parser.add_argument('-n', '--name', action="store", dest="name", help="name of inputs")
parser.add_argument('-f', '--fasta', action="store", dest="fasta", help="genome fasta file")
parser.add_argument('-w', '--window', action="store", type=int, dest="window", help="size of window around target sites")
params = parser.parse_args()

list_files = params.list.split(",")
name_list = params.name.split(",")
if params.fasta.endswith(".gz"):
  import gzip
  handle = gzip.open(params.fasta,"rt")
  genome = SeqIO.to_dict(SeqIO.parse(handle,'fasta'))
else:
  genome = SeqIO.to_dict(SeqIO.parse(open(params.fasta),'fasta'))

#for i in genome:
#	print("xxx",genome[i].seq)

retro = {}
if len(list_files) == len(name_list):
  for i in range(0,len(list_files)):
    if list_files[i].endswith(".gz"):
      import gzip
      fh = gzip.open(list_files[i])
    elif list_files[i] == "-":
      fh = sys.stdin
    else:
      fh = open(list_files[i])
    lines = fh.read().splitlines()
    r = None
    for j in lines:
      if r is None:
        r = cid.read_cell_id(j,w=params.window)
      else:
        n = cid.read_cell_id(j,w=params.window)
        r.join(n)
    retro[name_list[i]] = r.merge()
else:
  sys.exit("Number of files and names does not match")

for i in retro:
  for j in retro[i].chromosomes:
    for k in retro[i][j].as_df().to_numpy():
      sub = genome[k[0]].seq[k[1]:k[2]]
      print(">",i,".",j,".",k[1],".",k[2],sep="")
      print(sub)
