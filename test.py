import numpy as np
from bwampy import RefIndex

x = RefIndex("/home/skovaka/Dropbox/data/refs/zymo_hmw/all.fasta", True)
x# = RefIndex("/home/skovaka/Dropbox/data/refs/hg38/chr19.fa", True)

st = 5000000

seq = x.get_base(np.arange(st, st+50), False)

bitvec = np.array(x.pan_kmer_bitvec(
    seq, 11, [
    ["Escherichia_coli_chromosome", "Escherichia_coli_plasmid"], 
    ["BS.pilon.polished.v3.ST170922"], 
    ["Enterococcus_faecalis_complete_genome"]]))


print(np.array(bitvec))
print(len(seq))
print(bitvec.shape)
