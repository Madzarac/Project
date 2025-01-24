def process_file(input_file, output_file):
    """
    Initialize the estimator by parsing the SAM file
    :param input_file: Path to the fastq file
    :param output_file: Path to the .txt file with output of this function which shows correct mapping for each read
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        while True:
            header = infile.readline().strip()  # Read the header line
            if not header:
                break

            sequence = infile.readline().strip()

            header_parts = header.split()
            if len(header_parts) < 2:
                continue

            identifier = header_parts[0][1:]
            taxid_info_parts = header_parts[1].split(",")
            if len(taxid_info_parts) < 1:
                continue

            taxid_info = taxid_info_parts[0]

            outfile.write(f"{identifier}\t{taxid_info}\n")
            
input_file = "path/2bacteria4strains_c1.fastq"
output_file = "true_mappings.txt"

process_file(input_file, output_file)
