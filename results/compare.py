def compare_mappings(true_file, output_file, comparison_file):
    """
    Initialize the estimator by parsing the SAM file
    :param true_file: Path to the .txt file with true, expected mappings
    :param output_file: Path to the .txt file with the results of assign.py chosen mappings
    :param comparison_file: Path to the .txt file where the function will write comparison of results
    """
    true_mappings = {}

    with open(true_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue
            identifier, taxid = parts
            true_mappings[identifier] = taxid

    correct_counts = {}
    incorrect_counts = {}
    no_alignment_counts = {}

    all_taxids = set(true_mappings.values())
    for taxid in all_taxids:
        correct_counts[taxid] = 0

    comparison_lines = []

    with open(output_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                identifier, taxid = parts

                if identifier in true_mappings:
                    true_taxid = true_mappings[identifier]
                    if true_taxid == taxid:
                        is_correct = "Correct"
                        correct_counts[true_taxid] += 1
                    else:
                        is_correct = f"Incorrect (Correct: {true_taxid})"
                        key = f"{taxid} -> {true_taxid}"
                        incorrect_counts[key] = incorrect_counts.get(key, 0) + 1
                else:
                    is_correct = "Incorrect (No matching identifier in true mappings)"

                comparison_lines.append(f"{identifier}\t{taxid}\t{is_correct}\n")

            elif len(parts) == 1 and parts[0].endswith("NO ALIGNMENT"):
                identifier = parts[0].split()[0]
                if identifier in true_mappings:
                    true_taxid = true_mappings[identifier]
                    is_correct = f"No Alignment (Should align to: {true_taxid})"
                    no_alignment_counts[true_taxid] = no_alignment_counts.get(true_taxid, 0) + 1
                else:
                    is_correct = "No Alignment (No matching identifier in true mappings)"

                comparison_lines.append(f"{identifier}\tNO ALIGNMENT\t{is_correct}\n")

    with open(comparison_file, 'w') as outfile:
        outfile.write("Taxid\tCorrect Mappings\n")
        for taxid, count in correct_counts.items():
            outfile.write(f"{taxid}\t{count}\n")

        outfile.write("\nIncorrect Situations and Counts:\n")
        for key, count in incorrect_counts.items():
            outfile.write(f"{key}: {count} occurrences\n")

        outfile.write("\nNo Alignment Situations and Counts:\n")
        for taxid, count in no_alignment_counts.items():
            outfile.write(f"NO ALIGNMENT -> {taxid}: {count} occurrences\n")

        outfile.write("\nComparison Results:\n")
        for line in comparison_lines:
            outfile.write(line)

true_file = "true_mappings.txt"
output_file = "../2bac4strain.txt"
comparison_file = "comparison.txt"

compare_mappings(true_file, output_file, comparison_file)
