# importing the required packages for the functions script
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline

# function 1: to fetch the data from NCBI database
def DataFetch(id, format, file_name):
    Entrez.email = "g.elsaid2194@gmail.com"
    data = Entrez.efetch(db='nucleotide', id=id, rettype=format, retmode='text')
    seq = SeqIO.read(data, format=format)
    file = open(file_name + '.' + format, 'a')
    SeqIO.write(sequences=seq, handle=file, format=format)
    file.close()

# function 2: to calculate the GC percent of a given sequence
def gc_percent(seq):
    gc = seq.count('C') + seq.count('c') + seq.count('G') + seq.count('g')
    # formatting the GC content to 2 decimals only
    gc_per = float('{:.2f}'.format(gc / len(seq) * 100))
    return gc_per

# function 3: to calculate each nucleotide percent of a given sequence
def nuc_percent(seq):
    a_content = seq.count('A') + seq.count('a')
    # formatting the A content to 2 decimals only
    a_per = float('{:.2f}'.format(a_content / len(seq) * 100))
    c_content = seq.count('C') + seq.count('c')
    # formatting the C content to 2 decimals only
    c_per = float('{:.2f}'.format(c_content / len(seq) * 100))
    g_content = seq.count('G') + seq.count('g')
    # formatting the G content to 2 decimals only
    g_per = float('{:.2f}'.format(g_content / len(seq) * 100))
    t_content = seq.count('T') + seq.count('t')
    # formatting the T content to 2 decimals only
    t_per = float('{:.2f}'.format(t_content / len(seq) * 100))
    return a_per, c_per, g_per, t_per

# function 4: to perform ClustalW alignment of a given fasta file
def align(file):
    clustal_exe = 'C:\Program Files (x86)\ClustalW2\clustalw2'
    clustal_cline = ClustalwCommandline(clustal_exe, infile=file)
    clustal_cline()

# function 5: to get the consensus sequence of a given alignment file
def get_align_consensus(file, threshold):
    alignment = AlignIO.read(file, format='clustal')
    summary_align = AlignInfo.SummaryInfo(alignment)
    dissimilar = summary_align.dumb_consensus(threshold=threshold, require_multiple=10)
    return dissimilar
