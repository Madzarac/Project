def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        while True:
            header = infile.readline().strip()  # Read the header line
            if not header:
                break  # End of file

            sequence = infile.readline().strip()  # Read the sequence line

            # Split the header to extract the required parts
            header_parts = header.split()
            if len(header_parts) < 2:
                continue  # Skip lines with unexpected format

            identifier = header_parts[0][1:]  # Remove the leading '@' from the identifier
            taxid_info_parts = header_parts[1].split(",")
            if len(taxid_info_parts) < 1:
                continue  # Skip lines with unexpected format

            taxid_info = taxid_info_parts[0]  # Extract the taxid info

            # Write the processed line to the output file
            outfile.write(f"{identifier}\t{taxid_info}\n")
            
input_file = "path/2bacteria4strains_c1.fastq"
output_file = "true_mappings.txt"

process_file(input_file, output_file)
