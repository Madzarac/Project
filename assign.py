from collections import defaultdict
import pysam
from genericpath import samefile
import sys

class AbundanceEstimator:
    def __init__(self, sam_file):
        """
        Initialize the estimator by parsing the SAM file.
        :param sam_file: Path to the SAM file containing alignment information.
        """
        self.references = {}                       # {reference_name: reference_length (LN)}
        self.read_assignments = defaultdict(list)  # {read_id: [(reference_name, mapping_score (AS))]}
        self.strain_abundance = {}                 # {reference_name: initial abundance estimate}
        self.strain_coverage = {}                  # {reference_name: coverage fraction}

        self._parse_sam_file(sam_file)
        self._initialize_abundances()

    def _parse_sam_file(self, sam_file):
        """
        Parse the SAM file to extract reference lengths (from @SQ) and read mappings (AS field).
        """
        samfile = pysam.AlignmentFile(sam_file, "r")

        # checks header for references length
        for ref_info in samfile.header['SQ']:
            ref_name = ref_info['SN']
            ref_len = ref_info['LN']
            self.references[ref_name] = ref_len

        # extracts AS field
        for read in samfile:
            read_id = read.query_name
            ref_name = read.reference_name
            if read.flag != 4:
              if ref_name != None:  # Checks if read is mapped
                if read.has_tag('AS'):
                  as_field = read.get_tag('AS')
                  if as_field > 0:
                    if read_id not in self.read_assignments:
                      self.read_assignments[read_id] = []
                    self.read_assignments[read_id].append((ref_name, as_field))
                  else:
                    print(f"{read_id} AS <= 0: {as_field}")
                else:
                  print(f"{read_id} has NO AS tag")

        samfile.close()

    def _initialize_abundances(self):
        """
        Compute initial abundance estimates based on mapping scores and reference lengths.
        """
        strain_scores = defaultdict(float)

        for read_id, mappings in self.read_assignments.items():
            total_score = sum(score for _, score in mappings) # sums AS scores for all mapings of the read
            for ref_name, score in mappings:
                ref_len = self.references[ref_name]
                likelihood = score / ref_len  # Normalize by reference length, as longer reference could have better score because its longer
                strain_scores[ref_name] += likelihood / total_score # splits likelihood among all possible mappings

        # Normalize abundances to sum to 1
        total_abundance = sum(strain_scores.values())
        for ref_name, score in strain_scores.items():
            self.strain_abundance[ref_name] = score / total_abundance
            self.strain_coverage[ref_name] = strain_scores[ref_name] / total_abundance

    def em_algorithm(self, max_iter=300, eps=1e-3, min_abundance=0.1):
        """
        Run the EM algorithm to refine strain abundances.
        """
        strain_ids = list(self.strain_abundance.keys())            # list of all references (strains)
        new_abundance = {strain: 0.0 for strain in strain_ids}     # dict that stores updated abundance per strain in each iteration
        valid_strains = {strain: True for strain in strain_ids}    # a boolean flag for each strain

        for iteration in range(max_iter):
            # M-step: Update strain counts using mappings
            for strain in new_abundance:
                new_abundance[strain] = 0.0

            for read_id, mappings in self.read_assignments.items():
                total_prob = 0.0
                probs = []

                # Compute probability of this read coming from each reference
                for ref_name, score in mappings:
                    if valid_strains[ref_name]:
                        prob = score * self.strain_abundance[ref_name] * self.strain_coverage[ref_name]
                        probs.append((ref_name, prob))
                        total_prob += prob

                # Each read is fractionally assigned to all valid references it maps to, based on how probable each i
                if total_prob > 0:
                    for ref_name, prob in probs:
                        new_abundance[ref_name] += prob / total_prob

            # E-step: Normalize and check for convergence
            total_abundance = sum(new_abundance.values())
            max_diff = 0.0
            converged = True

            # Normalize new abundance values to sum to 1
            for strain in self.strain_abundance:
                if valid_strains[strain]:
                    diff = abs(new_abundance[strain] / total_abundance - self.strain_abundance[strain])
                    max_diff = max(max_diff, diff)
                    if diff > eps:  # If all strains have changed less than eps, we stop, othervise continue
                        converged = False
                    self.strain_abundance[strain] = new_abundance[strain] / total_abundance

            # Every 10 iterations, call apply_set_cover()
            if iteration % 10 == 0:
                self.apply_set_cover(valid_strains, min_abundance)

            if converged:
                break

    def apply_set_cover(self, valid_strains, min_abundance):
        """
        Eliminates strains with abundance smaller than min_abundance, unless they explain a read uniquely
        """
        removable_strains = {strain for strain, abundance in self.strain_abundance.items()
                             if valid_strains[strain] and abundance < min_abundance}

        # Prevent removing unique strains
        for read_id, mappings in self.read_assignments.items():
            unique_strains = {ref_name for ref_name, _ in mappings if valid_strains[ref_name]}
            if len(unique_strains) == 1:
                removable_strains.discard(next(iter(unique_strains)))

        for strain in removable_strains:
            valid_strains[strain] = False



class MoraAssignment:
    def __init__(self, references, reads):
        self.references = references                                                                                 # dictionary, {ref_name: abundance}
        self.reads = dict(sorted(reads.items(), key=lambda item: max(score for _, score in item[1]), reverse=True))  # sorted by highest scores
        self.assignments = {ref: [] for ref in references}                                                           # chosen mappings
        self.reference_capacity = {ref: round(references[ref] * len(self.reads)) for ref in references}              # reference capacities
        self.unassigned_reads = []                                                                                   # to store unassigned reads for reprocessing
        print("Reference Capacities Calculation:")
        for ref in references:
            calculated_capacity = round(references[ref] * len(self.reads))
            print(f"Reference {ref}: {references[ref]} * {len(self.reads)} = {calculated_capacity}")

    def calculate_priority(self, read_id):
        """
        Calculates priority for each read.
        """
        # Sort mappings of the read by score, best to worst
        read_mappings = sorted(self.reads[read_id], key=lambda x: x[1], reverse=True)
        best_score = read_mappings[0][1]
        if len(read_mappings) == 1:  # Priority 1 if only one mapping exists
            return 1
        second_best_score = read_mappings[1][1]
        if second_best_score / best_score < 0.5:  # Priority 2 if ratio is below 0.5
            return 2
        return 3  # Priority 3 for others

    def assign_read(self, read_id, priority):
        """
        Assigns reads with priority 1. 
        For 2 and 3 assigning to best reference if there is space, if not add to unassigned_reads for later
        """
        read_mappings = self.reads[read_id]
        if priority == 1:
            # Priority 1: Assign to the unique reference
            ref, score = read_mappings[0]
            if self.reference_capacity[ref] > 0:
                self.assignments[ref].append((read_id, score))
                self.reference_capacity[ref] -= 1
                return True
        else:
            # For Priority 2 or 3: Add mappings to the unassigned_reads list for reprocessing
            ref, score = read_mappings[0]
            if self.reference_capacity[ref] > 0: # and self.reference_capacity == 2:
                self.assignments[ref].append((read_id, score))
                self.reference_capacity[ref] -= 1
                return True
            else:
                # For Priority 2 or 3: Add mappings to the unassigned_mappings list for global reprocessing
                self.unassigned_reads.append(read_id)
                return False

    def reprocess_unassigned_mappings(self):
        """
        Using global best scores, assign read greediy, favoring high-score mappings as long as theres space for reference
        """
        # Collect all mappings for unassigned reads
        all_mappings = []
        for read_id in self.unassigned_reads:
            for ref, score in self.reads[read_id]:
                all_mappings.append((read_id, ref, score))

        # Sort all mappings globally by score (highest to lowest)
        all_mappings.sort(key=lambda x: x[2], reverse=True)
        # Try to assign the reads based on global sorted scores
        remaining_unassigned_reads = set(self.unassigned_reads)
        for read_id, ref, score in all_mappings:
            if read_id in remaining_unassigned_reads and self.reference_capacity[ref] > 0:
                self.assignments[ref].append((read_id, score))
                self.reference_capacity[ref] -= 1
                remaining_unassigned_reads.remove(read_id)

        # Update the list of unassigned reads
        self.unassigned_reads = list(remaining_unassigned_reads)


    def reassign_read(self, read_id):
        """
        Rearanges reads by removing lower quality mappings to another compatible reference.
        """
        # Try to "open up space" for the given read
        read_mappings = self.reads[read_id]
        for ref, score in read_mappings:
            # Check if the reference has any assigned reads
            if len(self.assignments[ref]) > 0:
                for assigned_read, assigned_score in self.assignments[ref]:
                    # If a lower-scoring read is found, reassign it
                    if assigned_score < score:
                        for new_ref, new_score in self.reads[assigned_read]:
                            if self.reference_capacity[new_ref] > 0:
                                # Reassign the lower-scoring read to a new reference
                                self.assignments[new_ref].append((assigned_read, assigned_score))
                                self.assignments[ref].remove((assigned_read, assigned_score))
                                self.reference_capacity[ref] += 1
                                self.reference_capacity[new_ref] -= 1
                                # Assign the current read to the original reference
                                self.assignments[ref].append((read_id, score))
                                return True
        return False

    def assign_reads(self):
        """
        Orchestrates assignment process.
        """
        # First pass: Attempt to assign reads based on priority
        for read_id in self.reads:
            priority = self.calculate_priority(read_id)
            self.assign_read(read_id, priority)

        # Second pass: Reprocess unassigned mappings
        self.reprocess_unassigned_mappings()

        # Third pass: Try to reassign unassigned reads by opening up space
        remaining_unassigned_reads = list(self.unassigned_reads)
        for read_id in remaining_unassigned_reads:
            if self.reassign_read(read_id):
                self.unassigned_reads.remove(read_id)



if __name__ == "__main__":
    samfile = input("Input path to SAM file: ").strip()
    sam_file = pysam.AlignmentFile(samfile, "r")
    
    filtered_sam = ""
    no_alignment = ""

    for read in sam_file:
        if read.flag != 4:
            filtered_sam += str(read) + "\n"
        else:
            no_alignment += str(read.query_name) + " NO ALIGNMENT" + "\n"

    estimator = AbundanceEstimator(samfile)
    estimator.em_algorithm()

    mora = MoraAssignment(estimator.strain_abundance, estimator.read_assignments)
    mora.assign_reads()

    result = ""
    for reference, reads in mora.assignments.items():
        for read, score in reads:
            result += f"{read}\t{reference}\n"

    result += no_alignment

    output_filename = input("Input path to txt file for storing results: ").strip()
    
    with open(output_filename, 'w') as file:
        file.write(result)

    print(f"Results saved in file '{output_filename}'.")
