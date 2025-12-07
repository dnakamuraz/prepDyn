import sys

# CONFIGURATION
input_file = "21_domD1.fasta"   # Your file
chunk_size = 1000            # Size of window
output_prefix = "window"     # Prefix for output files

def read_fasta(filename):
    """Reads fasta into a dictionary {header: sequence}"""
    seqs = {}
    name = None
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    name = line[1:]
                    seqs[name] = []
                else:
                    seqs[name].append(line)
    except IOError:
        print("Error: Could not find or read file: " + filename)
        sys.exit(1)

    # Join lists into strings
    for name in seqs:
        seqs[name] = "".join(seqs[name])
    return seqs

# 1. Load Data
sequences = read_fasta(input_file)
if not sequences:
    print("Error: No sequences found. Check your input file format.")
    sys.exit(1)

taxa = list(sequences.keys())

# Get alignment length from the first sequence
aln_len = len(sequences[taxa[0]])
print("Alignment Length: {} bp | Taxa: {}".format(aln_len, len(taxa)))

# 2. Slice and Write
num_files = 0
for start in range(0, aln_len, chunk_size):
    end = min(start + chunk_size, aln_len)
    num_files += 1
    
    # Old style formatting
    filename = "{}_{}.fasta".format(output_prefix, num_files)
    
    with open(filename, 'w') as out:
        for taxon in taxa:
            # Extract the slice for this taxon
            seq_slice = sequences[taxon][start:end]
            out.write(">{}\n{}\n".format(taxon, seq_slice))
            
    print("Wrote {} (Sites {}-{})".format(filename, start+1, end))

print("Done.")
