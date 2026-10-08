"""
Filter sequences from FASTA alignment files based on:
- Sequences <50% of longest ungapped length
- Species duplicates with lower % identity to REFERENCE
Outputs filtered sequences and logging information.
"""

import ntpath
import os
import sys
from difflib import SequenceMatcher

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fafile = sys.argv[1]

def similar(seq_a: str, seq_b: str) -> float:
    """Calculate sequence similarity ratio between two strings."""
    return SequenceMatcher(None, seq_a, seq_b).ratio()

def isclose(a: float, b: float, rel_tol: float = 1e-09, abs_tol: float = 0.0) -> bool:
    """Check if two floats are approximately equal."""
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

# Parse filename for gene and transcript info
fabase = ntpath.basename(fafile)
falist = fabase.split(".")
transcript = falist[0]
gene = falist[1]

# Create output filenames
base_dir = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/Summary_CDS_Phase3Alns"
protdir = "/share/ceph/wym219group/shared/data/filtered_alignments"

separator = '.'
outaln = os.path.join(protdir, separator.join([transcript, gene, 'filt', 'fa']))
outlog = os.path.join(base_dir, separator.join([transcript, gene, 'filtlog', 'txt']))
erlog = os.path.join(base_dir, separator.join([transcript, gene, 'erlog', 'txt']))
redospec = os.path.join(base_dir, 'FilesWithNoHumanRef.txt')
noseq = os.path.join(base_dir, 'EmptyAlignments.txt')

# Dictionaries to store info on each species
speciesseqdict = {}  # keep sequence for each record
specieslendict = {}  # track longest ungapped for dups
specieshumdistdict = {}  # track % identity with human for dups
speciesdescdict = {}  # track descriptions to keep

filteredlen = []  # save details of records filtered for length
filtereddup = []  # save details of duplicate species records filtered
longestungapped = 0

# Inform where the filtered alignment will be written
print(f"Writing filtered alignment to: {outaln}")
filt_fa = open(outaln, "w")

# Read alignment file
try:
    align = AlignIO.read(fafile, "fasta")
except Exception:
    with open(noseq, 'a') as nhandle:
        nhandle.write(fafile + '\n')
    sys.exit("Empty fafile: " + fafile)

# Parse first for human/reference sequence
humseq = None
for record in align:
    if record.id == "REFERENCE":
        humseq = str(record.seq)
        break

if humseq is None:
    with open(redospec, 'a') as rhandle:
        rhandle.write(fafile + '\n')
    sys.exit("No human reference: " + fafile)

# Process sequences
for record in align:
    seq_str = str(record.seq)
    sequence_length = len(seq_str)
    ungapped_length = len(seq_str.replace("-", ""))
    
    if ungapped_length > longestungapped:
        longestungapped = ungapped_length

    if float(ungapped_length) < 0.5 * sequence_length:
        # Filter short sequences
        filteredlen.append(record)
        continue

    if record.id in specieshumdistdict:
        # Compare seq as string to humseq (both plain strings)
        newsim = similar(seq_str, humseq)
        if newsim > specieshumdistdict[record.id]:
            specieshumdistdict[record.id] = newsim
            specieslendict[record.id] = ungapped_length
            speciesdescdict[record.id] = record.description
            speciesseqdict[record.id] = seq_str
        elif isclose(newsim, specieshumdistdict[record.id]):
            # Choose longest sequence if similarity is essentially equal
            if ungapped_length > specieslendict[record.id]:
                specieshumdistdict[record.id] = newsim
                specieslendict[record.id] = ungapped_length
                speciesdescdict[record.id] = record.description
                speciesseqdict[record.id] = seq_str
            else:
                # Log equivalent records
                with open(erlog, "a") as er_handle:
                    er_handle.write(
                        f'Multiple equivalent records for {record.id} in {fafile}: '
                        f'{speciesdescdict[record.id]} and {record.description}; '
                        f'including just the first.\n'
                    )
                filtereddup.append(record)
    else:
        specieshumdistdict[record.id] = similar(seq_str, humseq)
        specieslendict[record.id] = ungapped_length
        speciesdescdict[record.id] = record.description
        speciesseqdict[record.id] = seq_str

# Write output sequences after subsetting alignment
for key in specieshumdistdict.keys():
    nt_seq_str = speciesseqdict[key]
    # Replace uncommon characters with 'N'
    tmpdna = nt_seq_str.replace('!', 'N')
    nt_seq = Seq(tmpdna)
    
    nt_record = SeqRecord(nt_seq, id=key, description=speciesdescdict[key])
    SeqIO.write(nt_record, filt_fa, 'fasta')

filt_fa.close()

# Confirm the file was written
print(f"Wrote filtered alignment file: {outaln}")

# Write to output log
if len(filteredlen) > 0 or len(filtereddup) > 0:
    with open(outlog, 'a') as ohandle:
        if len(filteredlen) > 0:
            ohandle.write(f'Sequences filtered from {fafile} for <50% length:\n')
            for fl in filteredlen:
                ohandle.write(f'{fl.id} {fl.description}\n')
        if len(filtereddup) > 0:
            ohandle.write(f'Duplicate sequences filtered from {fafile}:\n')
            for fd in filtereddup:
                ohandle.write(f'{fd.id} {fd.description}\n')
