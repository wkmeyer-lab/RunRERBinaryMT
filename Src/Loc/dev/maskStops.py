import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from tqdm import tqdm

# --- CONFIGURATION ---
input_dir = "/share/ceph/wym219group/shared/data/filtered_alignments"
output_dir = "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/masked_genes"
log_file = "masked_genes_log.txt"

valid_extensions = (".fasta", ".fas", ".fna", ".fa", ".aln")

os.makedirs(output_dir, exist_ok=True)
modified_files = []
files_processed_count = 0 

def get_sanitized_seq_str(seq_obj):
    """
    Converts sequence to string and replaces illegal characters 
    (! and ?) with N so Biopython doesn't crash during translation.
    """
    s = str(seq_obj)
    #replace ! and ? with N (Unknown)
    if "!" in s or "?" in s:
        s = s.replace("!", "N").replace("?", "N")
    return s

def has_internal_stops(seq_obj):
    #sanitize the sequence string (replace ! with N)
    clean_seq_str = get_sanitized_seq_str(seq_obj)
    clean_seq = Seq(clean_seq_str)
    
    try:
        #translate
        protein = clean_seq.translate(to_stop=False)
    except Exception:
        #if it still fails (e.g. length not multiple of 3), skip or flag
        return False 
    
    #check for stop (*) anywhere before the final character
    if "*" in protein[:-1]:
        return True
    return False

def mask_stops(seq_obj):
    #sanitize first to get the protein coordinates
    clean_seq_str = get_sanitized_seq_str(seq_obj)
    clean_seq = Seq(clean_seq_str)
    
    #translate the CLEAN sequence to find stop positions
    protein = clean_seq.translate(to_stop=False)
    
    #create a mutable list of the original sequence to apply masks to
    mutable_seq = list(str(seq_obj))
    
    for i, aa in enumerate(protein):
        #if stop found before the end
        if aa == "*" and i < len(protein) - 1:
            dna_start = i * 3
            #mask with dashes in the output sequence
            mutable_seq[dna_start:dna_start+3] = ['-', '-', '-']
            
    return "".join(mutable_seq)

print(f"Scanning directory: {input_dir}")

all_files = os.listdir(input_dir)
print(f"Total files in folder: {len(all_files)}")

for filename in tqdm(all_files):
    #check extension (case insensitive)
    if not filename.lower().endswith(valid_extensions): 
        continue
        
    files_processed_count += 1
    filepath = os.path.join(input_dir, filename)
    outpath = os.path.join(output_dir, filename)
    
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        continue

    file_needs_masking = False
    cleaned_records = []

    for record in records:
        #check if masking is necessary
        if has_internal_stops(record.seq):
            file_needs_masking = True
            new_seq_str = mask_stops(record.seq)
            new_record = SeqRecord(
                Seq(new_seq_str),
                id=record.id,
                description=record.description
            )
            cleaned_records.append(new_record)
        else:
            cleaned_records.append(record)

    #write output if changes were made
    if cleaned_records:
        SeqIO.write(cleaned_records, outpath, "fasta")
    
    if file_needs_masking:
        modified_files.append(filename)

# --- REPORT ---
print("-" * 30)
print(f"Files actually processed: {files_processed_count}")
print(f"Files masked: {len(modified_files)}")

with open(log_file, "w") as f:
    f.write(f"Total files scanned: {len(all_files)}\n")
    f.write(f"Files processed: {files_processed_count}\n")
    f.write(f"Files masked: {len(modified_files)}\n")
    f.write("-" * 30 + "\n")
    for fname in modified_files:
        f.write(f"{fname}\n")

print(f"Log saved to: {log_file}")
