#This script will find which contigs have telomeric repeats for the reference genome assembly.
# telomere_search.py

from Bio import SeqIO
print("Biopython is available.")

# Define the telomeric repeats and search parameters
repeat_patterns = ["CCCTAA", "TTAGGG"]
min_repeat_count = 5  # Minimum of 5 consecutive repeats for detection
repeat_length = 6     # Length of each repeat
search_region_length = 100  # Search within 100 bp of each contig end

# Paths to input and output files
input_file = "/nfs4/my-mgrid-s8/mykopat-proj1/mykopat-c-flacc/sort_ref_genome/renamed_sort_nucleargenome.fasta"
output_txt = "/nfs4/my-mgrid-s8/mykopat-proj1/mykopat-c-flacc/sort_ref_genome/telomere_repeats_report.txt"
output_csv = "/nfs4/my-mgrid-s8/mykopat-proj1/mykopat-c-flacc/sort_ref_genome/telomere_repeats_summary.csv"


# Function to count repeat occurrences in a given region
def count_telomeric_repeats(region, pattern):
    max_count = 0
    current_count = 0
    i = 0
    while i <= len(region) - repeat_length:
        if region[i:i + repeat_length] == pattern:
            current_count += 1
            i += repeat_length
        else:
            max_count = max(max_count, current_count)
            current_count = 0
            i += 1
    return max(max_count, current_count)


# Open the output files for writing results
with open(output_txt, "w") as report, open(output_csv, "w") as csv_report:
    # Write CSV header
    csv_report.write("Scaffold,Start_Telomere_Repeats,End_Telomere_Repeats\n")
    report.write("Telomeric Repeat Detection Report\n")
    report.write("================================\n")

    # Parse each scaffold in the FASTA file
    with open(input_file) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            scaffold_id = record.id
            sequence = str(record.seq)
            seq_len = len(sequence)

            # Define start and end regions to search
            start_region = sequence[:search_region_length]
            end_region = sequence[-search_region_length:]

            # Count repeats in start and end regions
            start_repeat_count = max(count_telomeric_repeats(start_region, pattern) for pattern in repeat_patterns)
            end_repeat_count = max(count_telomeric_repeats(end_region, pattern) for pattern in repeat_patterns)

            # Log to TXT report
            if start_repeat_count >= min_repeat_count or end_repeat_count >= min_repeat_count:
                report.write(f"Telomeric sequence detected in scaffold: {scaffold_id}\n")
                if start_repeat_count >= min_repeat_count:
                    report.write(
                        f" - Telomeric sequence found at start with {start_repeat_count} repeats within first {search_region_length} bp\n")
                if end_repeat_count >= min_repeat_count:
                    report.write(
                        f" - Telomeric sequence found at end with {end_repeat_count} repeats within last {search_region_length} bp\n")
            else:
                report.write(f"No telomeric sequence detected in scaffold: {scaffold_id}\n")

            # Write to CSV summary
            csv_report.write(f"{scaffold_id},{start_repeat_count},{end_repeat_count}\n")
