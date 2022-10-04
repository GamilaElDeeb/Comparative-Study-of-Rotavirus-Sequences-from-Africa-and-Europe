# importing the required packages for the script
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np
import ProjectFunctions as f


# 1- Downloading the 10 RotaVirus sequences from Africa
# Country: Mali
# Trial link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4594120/
# 10 Africa Samples IDs
africa_ids = ['KP883218.1', 'KP883207.1', 'KP883196.1', 'KP883185.1', 'KP883174.1', 'KP883163.1',
              'KP883152.1 ', 'KP883141.1', 'KP883130.1', 'KP883119.1']
# fetching samples from NCBI virus database using Entrez
for h in africa_ids:
    # DataFetch is a function that takes 3 arguments 1-NCBI id 2-file format 3-file name
    # and return the sequence in the file format and save it the file name
    # Usage: f.DataFetch(id, format, file_name)
    f.DataFetch(h, 'fasta', 'africa_samples')


# 2- Making the consensus sequence from Africa samples
# reading the africa samples file and creating a list of sequences to get the consensus sequence
file1 = SeqIO.parse('africa_samples.fasta', format='fasta')
africa_seqs = []
for l in file1:
    africa_seqs.append(str(l.seq))
# creating profile matrix for nucleotides in each sequence
n = len(africa_seqs[0])
profile = {'T': [0]*n, 'G': [0]*n, 'C': [0]*n, 'A': [0]*n}
for seq in africa_seqs:
    for i, char in enumerate(seq):
        profile[char][i] += 1
# getting the consensus sequence
afr_consensus = ""
for i in range(n):
    max_count = 0
    max_nt = 'x'
    for nt in "ACGT":
        if profile[nt][i] > max_count:
            max_count = profile[nt][i]
            max_nt = nt
    afr_consensus += max_nt

# saving the africa consensus sequence in a fasta file
a = Seq(afr_consensus)
consensus_rec = SeqRecord(a, id='1', name='Consensus',  description='Africa_consensus',
                          annotations={'molecule_type': 'DNA'})
SeqIO.write(sequences=consensus_rec, handle='africa_consensus.fasta', format='fasta')

# 3- Downloading the 10 RotaVirus sequences from Europe
# Country : Hungary
# Trial link: https://www.sciencedirect.com/science/article/abs/pii/S156713481400344X?via%3Dihub
# 10 Europe Samples Ids
europe_ids = ['KJ919185.1', 'KJ919184.1', 'KJ919183.1', 'KJ919182.1', 'KJ919181.1', 'KJ919180.1',
              'KJ919179.1', 'KJ919178.1', 'KJ919177.1', 'KJ919176.1']
# fetching europe samples from NCBI virus database using Entrez
for h in europe_ids:
    # DataFetch is a function that takes 3 arguments 1-NCBI id 2-file format 3-file name
    # and return the sequence in the file format and save it the file name
    # Usage: f.DataFetch(id, format, file_name)
    f.DataFetch(h, 'fasta', 'europe_samples')

# 4- Preparing the alignment sequences
# combining Africa consensus sequence and Europe sequences in one file for alignment
# first reading files and appending sequences in one list
file2 = SeqIO.read('africa_consensus.fasta', format='fasta')
file3 = SeqIO.parse('europe_samples.fasta', format='fasta')
align_seqs = []
eur_align_ids = []
align_seqs.append(file2.seq)
for l in file3:
    align_seqs.append(l.seq)
    eur_align_ids.append(l.id)

# creating SeqRecord list from alignment sequences list
counter = 0
align_rec = []
for seq in align_seqs:
    if counter == 0:
        seq = SeqRecord(seq, id='AfricaReference', name='alignment', description='alignment seq' + str(counter),
                        annotations={'molecule_type': 'DNA'})
    else:
        seq = SeqRecord(seq, id=eur_align_ids[counter-1], name='alignment', description='alignment seq' + str(counter),
                        annotations={'molecule_type': 'DNA'})
    align_rec.append(seq)
    counter += 1

# creating the alignment sequences file
SeqIO.write(sequences=align_rec, handle='alignment_seqs.fasta', format='fasta')

# 5- performing the alignment between Africa consensus sequence and Europe samples
# align function takes a file of sequences in fasta format and performs ClustalW alignment between the sequences
# Usage: f.align(file)
f.align('alignment_seqs.fasta')

# 6- Making the phylogenetic tree between Africa consensus sequence and Europe samples
# building the phylogenetic tree
tree = Phylo.read('alignment_seqs.dnd', 'newick')
tree.ladderize()
Phylo.draw(tree)
# save the figure to continue running the code

# 7- Calculating the chemical constituents of Africa consensus and Europe samples
# calculating the nucleotides percents for Africa consensus and Europe samples
file4 = SeqIO.parse('alignment_seqs.fasta', format='fasta')
# opening a file to save the nucleotides percents
nuc_data = open('Nucleotides_content.csv', 'a')
# writing file header
nuc_data.write('Sample,GC%,A%,C%,G%,T%')
counter = 0
for l in file4:
    if counter == 0:
        # gc_percent is a function that takes a sequence string and calculates the GC percent of the sequence
        # Usage: f.gc_percent(seq)
        # nuc_percent is a function that take a sequence string, calculate each nucleotide percent
        # and returns the percentage of the four nucleotides in the order of (A, C, G, T)
        # Usage: f.nuc_percent(seq)
        # writing the Africa consensus Nucleotides percents to the file
        nuc_data.write('\nAfrica consensus,' + str(f.gc_percent(str(l.seq))) + ',' +
                       str(f.nuc_percent(str(l.seq))[0]) + ',' + str(f.nuc_percent(str(l.seq))[1]) +
                       ',' + str(f.nuc_percent(str(l.seq))[2]) + ',' + str(f.nuc_percent(str(l.seq))[3]))
    else:
        # writing the Europe samples nucleotides percents to the file
        nuc_data.write('\n' + eur_align_ids[counter - 1] + ',' + str(f.gc_percent(str(l.seq))) + ',' +
                       str(f.nuc_percent(str(l.seq))[0]) + ',' + str(f.nuc_percent(str(l.seq))[1]) +
                       ',' + str(f.nuc_percent(str(l.seq))[2]) + ',' + str(f.nuc_percent(str(l.seq))[3]))
    counter += 1

# 8- getting the consensus from Europe samples
# reading the Europe samples file and creating a list of sequences to get the consensus sequence
file6 = SeqIO.parse('europe_samples.fasta', format='fasta')
europe_seqs = []
for l in file6:
    europe_seqs.append(str(l.seq))
# creating profile matrix for nucleotides in each sequence
n1 = len(europe_seqs[0])
profile1 = {'T': [0]*n, 'G': [0]*n, 'C': [0]*n, 'A': [0]*n}
for seq in europe_seqs:
    for i, char in enumerate(seq):
        profile1[char][i] += 1
# getting the consensus sequence
eur_consensus = ""
for i in range(n1):
    max_count1 = 0
    max_nt1 = 'x'
    for nt in "ACGT":
        if profile1[nt][i] > max_count1:
            max_count1 = profile1[nt][i]
            max_nt1 = nt
    eur_consensus += max_nt1

# saving the europe consensus sequence in a fasta file
a1 = Seq(eur_consensus)
eur_cons_rec = SeqRecord(a1, id='1', name='Consensus',  description='Europe_consensus',
                         annotations={'molecule_type': 'DNA'})
SeqIO.write(sequences=eur_cons_rec, handle='europe_consensus.fasta', format='fasta')


# 9- Finding the dissimilar regions between Africa consensus and Europe consensus
# combining Africa consensus and Europe consensus in one file
# line 70 file2 = SeqIO.read('africa_consensus.fasta', format='fasta')
file7 = SeqIO.read('europe_consensus.fasta', format='fasta')
# making a list of both sequences as a SecRecord object to write them to a file
both_consens = []
afr_cons = SeqRecord(file2.seq, id='Africa Consensus', name='alignment', description='alignment seq 1',
                     annotations={'molecule_type': 'DNA'})
eur_cons = SeqRecord(file7.seq, id='Europe Consensus', name='alignment', description='alignment seq 2',
                     annotations={'molecule_type': 'DNA'})
both_consens.extend([afr_cons, eur_cons])
# writing the file of both consensuses
SeqIO.write(sequences=both_consens, handle='afr_eur_consensus.fasta', format='fasta')

# making the alignment between Africa consensus and Europe consensus
# align function takes a file of sequences in fasta format and performs ClustalW alignment between the sequences
# Usage: f.align(file)
f.align('afr_eur_consensus.fasta')

# getting the dissimilarity between the 2 consensus sequences
# get_align_consensus is a function that takes 2 arguments 1- alignment file name 2- threshold of the consensus
# and returns the consensus of the aligned sequences
# Usage: f.get_align_consensus(file, threshold)
print('Africa consensus and Europe consensus have ' +
      str(f.get_align_consensus('afr_eur_consensus.aln', 1).count('X')) + ' dissimilar nucleotides')

# getting the positions of the dissimilar nucleotides into a list
pos = 0
dis_pos = []
conserved = []
for nucl in f.get_align_consensus('afr_eur_consensus.aln', 1):
    pos += 1
    if nucl == 'X':
        dis_pos.append(pos)
    else:
        conserved.append(pos)
# saving the dissimilar nucleotides and conserved nucleotides in files
file8 = open('dissimilar_nucleotides.txt', 'a')
file9 = open('conserved_nucleotides.txt', 'a')
for m in dis_pos:
    file8.write('Dissimilar nucleotide no: ' + str(m) + '\n')
for n in conserved:
    file9.write('Conserved nucleotide no :' + str(n) + '\n')
print(' ')

# 10 - finding the functional interpretation of the dissimilar regions
# the coding region in the samples is from nucleotide 1 to nucleotide 528
# appending the coding region of Africa consensus sequence into triplets in a list
afr_amino_nuc = []
for aa in range(0, 528, 3):
    afr_amino_nuc.append(file2.seq[aa:aa+3])
# translating the coding region into a sequence of amino acids
afr_amino = []
for am in afr_amino_nuc:
    afr_amino.append(str(Seq.translate(am)))

# appending the coding region of Europe consensus sequence into triplets in a list
eur_amino_nuc = []
for bb in range(0, 528, 3):
    eur_amino_nuc.append(file7.seq[bb:bb + 3])
# translating the coding region into a sequence of amino acids
eur_amino = []
for ab in eur_amino_nuc:
    eur_amino.append(str(Seq.translate(ab)))
# turning the amino acids list into a Seq object
africa_ptn = Seq(''.join(afr_amino))
europe_ptn = Seq(''.join(eur_amino))
# appending the Africa amino acids and Europe amino acids in a list of SeqRecords to save them in a file
both_ptns = []
afr_ptn = SeqRecord(africa_ptn, id='Africa Consensus Protein', name='protein', description='Amino Acids',
                    annotations={'molecule_type': 'Protein'})
eur_ptn = SeqRecord(europe_ptn, id='Europe Consensus Protein', name='protein', description='Amino Acids',
                    annotations={'molecule_type': 'Protein'})
both_ptns.extend([afr_ptn, eur_ptn])
# saving the 2 consensuses proteins in one file
SeqIO.write(sequences=both_ptns, handle='afr_eur_proteins.fasta', format='fasta')
# making alignment between the 2 consensuses proteins to extract the dissimilarity between them
# align function takes a file of sequences in fasta format and performs ClustalW alignment between the sequences
# Usage: f.align(file)
f.align('afr_eur_proteins.fasta')

# getting the dissimilarity between the 2 protein consensus sequences
# get_align_consensus is a function that takes 2 arguments 1- alignment file name 2- threshold of the consensus
# and returns the consensus of the aligned sequences
# Usage: f.get_align_consensus(file, threshold)
print('Africa consensus protein and Europe consensus proteins have ' +
      str(f.get_align_consensus('afr_eur_proteins.aln', 1).count('X')) + ' dissimilar amino acids')

# getting the positions of the dissimilar amino acids into a list
position = 0
dis_amino = []
conserved_amino = []
dis_amino_afr = []
dis_amino_eur = []
for amino in f.get_align_consensus('afr_eur_proteins.aln', 1):
    position += 1
    if amino == 'X':
        dis_amino.append(position)
        dis_amino_afr.append(africa_ptn[position - 1])
        dis_amino_eur.append(europe_ptn[position - 1])
    else:
        conserved_amino.append(position)


# saving the dissimilar amino acids and conserved amino acids in files
file10 = open('dissimilar_amino_acids.txt', 'a')
file11 = open('conserved_amino_acids.txt', 'a')
for u in dis_amino:
    file10.write('Dissimilar amino acid no: ' + str(u) + '\n')
for k in conserved_amino:
    file11.write('Conserved amino acid no: ' + str(k) + '\n')

print(' ')
# printing to the screen the dissimilar amino acids between Africa and Europe
print('Africa dissimilar amino acids are', dis_amino_afr)
print('Europe dissimilar amino acids are', dis_amino_eur)


# 11- Drawing nucleotides percentage curves for africa reference sequence and europe samples
# determining the x-axis
x = ['A', 'C', 'G', 'T']
# creating a list of percentages to be used as y-axis
africa_nucs = []
europe_a = []
europe_c = []
europe_g = []
europe_t = []
file10 = SeqIO.parse('alignment_seqs.fasta', format='fasta')
count = 0
for k in file10:
    if count == 0:
        # nuc_percent is a function that take a sequence string, calculate each nucleotide percent
        # and returns the percentage of the four nucleotides in the order of (A, C, G, T)
        # appending the nucleotides percent to africa list
        africa_nucs.extend([f.nuc_percent(str(k.seq))[0], f.nuc_percent(str(k.seq))[1],
                            f.nuc_percent(str(k.seq))[2], f.nuc_percent(str(k.seq))[3]])
    else:
        # appending the nucleotides percents to each europe nucleotide list
        europe_a.append(f.nuc_percent(str(k.seq))[0])
        europe_c.append(f.nuc_percent(str(k.seq))[1])
        europe_g.append(f.nuc_percent(str(k.seq))[2])
        europe_t.append(f.nuc_percent(str(k.seq))[3])
    count += 1

# calculating the average nucleotide percent in europe samples and appending the result in a list
eur_nucs = []
eur_a = np.array(europe_a)
eur_c = np.array(europe_c)
eur_g = np.array(europe_g)
eur_t = np.array(europe_t)
eur_nucs.extend([eur_a.mean(), eur_c.mean(), eur_g.mean(), eur_t.mean()])

# plotting the average nucleotide percentage in europe samples
fig = plt.figure()
fig.set_size_inches(7, 6, forward=True)
# plotting the average nucleotide percentage in europe samples
eur_ax = fig.add_subplot(2, 1, 1)
eur_ax.plot(x, eur_nucs, 'ro-')
eur_ax.set_xlabel('Nucleotides')
eur_ax.set_ylabel('Percentage')
eur_ax.set_ylim(10, 45)
eur_ax.set_title('Average Europe Samples Nucleotides Percentages')
# plotting the africa reference nucleotides percentage
afr_ax = fig.add_subplot(2, 1, 2)
afr_ax.plot(x, africa_nucs, 'bo-')
afr_ax.set_xlabel('Nucleotides')
afr_ax.set_ylabel('Percentage')
afr_ax.set_ylim(10, 45)
afr_ax.set_title('Africa Reference Nucleotides Percentages')
fig.subplots_adjust(hspace=0.5)
plt.show()
# save the figure to continue running the code

# plotting Africa Reference and Europe Samples average in the same figure
plt.xlabel('Nucleotides')
plt.ylabel('Percentage')
plt.title('Africa Reference And Average Europe Samples Nucleotides Percentage')
plt.ylim(10, 45)
plt.plot(x, eur_nucs, 'ro-', label='Europe')
plt.plot(x, africa_nucs, 'bo-', label='Africa')
plt.legend()
plt.show()
# save the figure to continue running the code
