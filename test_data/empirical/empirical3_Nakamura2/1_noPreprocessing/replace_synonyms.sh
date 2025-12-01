#!/bin/bash

# Define the synonym list file and the output file for sed commands
SYNONYM_LIST="Dendropsophus_synonym_list.txt"
SED_COMMANDS="replacements.sed"

# Check if the synonym list file exists
if [ ! -f "$SYNONYM_LIST" ]; then
    echo "Error: Synonym list file '$SYNONYM_LIST' not found."
    exit 1
fi

# 1. Generate the SED commands
echo "Generating replacement commands..."
# Use awk to read the synonym list
# $1 is the first string (the replacement name)
# The loop iterates from $2 to the last string (the names to be replaced)
awk '
{
    # Store the first string (the canonical name)
    canonical_name = $1
    
    # Iterate over the rest of the strings (the synonyms)
    for (i = 2; i <= NF; i++) {
        # Create a sed command: s/synonym/canonical_name/g
        # We use a non-standard delimiter (e.g., '#') for safety
        # in case names contain forward slashes, although this is rare in FASTA headers.
        print "s/" $i "/" canonical_name "/g"
    }
}' "$SYNONYM_LIST" > "$SED_COMMANDS"

# Check if any commands were generated
if [ ! -s "$SED_COMMANDS" ]; then
    echo "Warning: No replacement commands generated. Check if '$SYNONYM_LIST' is correctly formatted."
    exit 0
fi

# 2. Apply the replacements to all FASTA files
echo "Applying replacements to FASTA files..."
# Loop through all files ending in .fasta
for fasta_file in *.fasta; do
    # Check if the file actually exists (prevents running on literal string if no files found)
    if [ -f "$fasta_file" ]; then
        echo "Processing $fasta_file..."
        
        # Use sed with the generated commands. 
        # The -i flag performs the replacement "in-place" (i.e., modifies the original file).
        # We only apply the change to header lines (lines starting with >).
        sed -i.bak -f "$SED_COMMANDS" "$fasta_file"
        
        echo "Finished $fasta_file. A backup file (.bak) was created."
        # Optional: remove backup files after successful operation
        # rm "${fasta_file}.bak"
    fi
done

# 3. Cleanup
rm "$SED_COMMANDS"

echo "Replacement process complete! 🎉"
