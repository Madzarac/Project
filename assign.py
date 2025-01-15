import numpy as np
from collections import defaultdict
import pysam
from genericpath import samefile
import sys

class AbundanceEstimator:
    def __init__(self, sam_file):
        """
        Initialize the estimator by parsing the SAM file
        :param sam_file: Path to the SAM file containing alignment information
        """
        self.references = {}  # {reference_name: reference_length (LN)}
        self.read_assignments = defaultdict(list)  # {read_id: [(reference_name, mapping_score (AS))]}
        self.strain_abundance = {}  # {reference_name: initial abundance estimate}
        self.strain_coverage = {}  # {reference_name: coverage fraction}

        self._parse_sam_file(sam_file)
        self._initialize_abundances()

    def _parse_sam_file(self, sam_file):
        """
        Parse the SAM file to extract reference lengths (from @SQ) and read mappings (AS field)
        """
        samfile = pysam.AlignmentFile(sam_file, "r")

        #header
        for ref_info in samfile.header['SQ']:
            ref_name = ref_info['SN']
            ref_len = ref_info['LN']
            self.references[ref_name] = ref_len

        for read in samfile:
            read_id = read.query_name
            ref_name = read.reference_name
            if read.flag != 4:
              if ref_name != None:  # Checks if read is mapped
                  as_field = read.get_tag('AS') if read.has_tag('AS') else 0
                  if as_field > 0:  # Use only mapped reads with positive AS score
                      if read_id not in self.read_assignments:
                          self.read_assignments[read_id] = []
                      self.read_assignments[read_id].append((ref_name, as_field))

        samfile.close()

    def _initialize_abundances(self):
        """
        Compute initial abundance estimates based on mapping scores and reference lengths
        """
        strain_scores = defaultdict(float)

        for read_id, mappings in self.read_assignments.items():
            total_score = sum(score for _, score in mappings)
            for ref_name, score in mappings:
                ref_len = self.references[ref_name]
                likelihood = score / ref_len  # Normalize by reference length
                strain_scores[ref_name] += likelihood / total_score

        # Normalize abundances to sum to 1
        total_abundance = sum(strain_scores.values())
        for ref_name, score in strain_scores.items():
            self.strain_abundance[ref_name] = score / total_abundance
            self.strain_coverage[ref_name] = strain_scores[ref_name] / total_abundance

    def em_algorithm(self, max_iter=300, eps=1e-3, min_abundance=0.1):
        """
        Run the EM algorithm to refine strain abundances
        """
        strain_ids = list(self.strain_abundance.keys())
        new_abundance = {strain: 0.0 for strain in strain_ids}
        valid_strains = {strain: True for strain in strain_ids}

        for iteration in range(max_iter):
            # M-step: Update strain counts using mappings
            for strain in new_abundance:
                new_abundance[strain] = 0.0

            for read_id, mappings in self.read_assignments.items():
                total_prob = 0.0
                probs = []

                for ref_name, score in mappings:
                    if valid_strains[ref_name]:
                        prob = score * self.strain_abundance[ref_name] * self.strain_coverage[ref_name]
                        probs.append((ref_name, prob))
                        total_prob += prob

                if total_prob > 0:
                    for ref_name, prob in probs:
                        new_abundance[ref_name] += prob / total_prob

            # E-step: Normalize and check for convergence
            total_abundance = sum(new_abundance.values())
            max_diff = 0.0
            converged = True

            for strain in self.strain_abundance:
                if valid_strains[strain]:
                    diff = abs(new_abundance[strain] / total_abundance - self.strain_abundance[strain])
                    max_diff = max(max_diff, diff)
                    if diff > eps:
                        converged = False
                    self.strain_abundance[strain] = new_abundance[strain] / total_abundance

            if iteration % 10 == 0:
                self.apply_set_cover(valid_strains, min_abundance)

            if converged:
                break

    def apply_set_cover(self, valid_strains, min_abundance):
        """
        Reduce valid strains using a set cover heuristic
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
        """
        Initialize mora for read assignment
        :param references: reference and abundance pairs in dictionary format
        :param reads: reads and theri possible mappings
        """
        self.references = references
        self.reads = dict(sorted(reads.items(), key=lambda item: max(score for _, score in item[1]), reverse=True))  # sorted by highest scores
        self.assignments = {ref: [] for ref in references}  # chosen mappings
        self.reference_capacity = {ref: round(references[ref] * len(self.reads)) for ref in references}  # reference capacities
        self.unassigned_reads = []  # to store unassigned reads for reprocessing

    def calculate_priority(self, read_id):
        """
        Calculates priority (1, 2 or 3) for given read
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
        Assigns mappings to reads priority 1 and 2 if possible
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
            # Priority 2: reads are sorted and then assigned to the reference with the best mapping score if that reference has space. If the reference is at full capacity, the read is relabeled as a priority 3 read
            ref, score = read_mappings[0]
            if self.reference_capacity[ref] > 0 and self.reference_capacity == 2:
                self.assignments[ref].append((read_id, score))
                self.reference_capacity[ref] -= 1
                return True
            else:
                #Priority 3: Add mappings to the unassigned_mappings list
                self.unassigned_reads.append(read_id)
                return False

    def reprocess_unassigned_mappings(self):
        """
        Assigns mappings to reads priority 3 if possible
        """
        # Collect all mappings for unassigned reads
        all_mappings = []
        for read_id in self.unassigned_reads:
            for ref, score in self.reads[read_id]:
                all_mappings.append((read_id, ref, score))

        # Sort all mappings by score (highest to lowest)
        all_mappings.sort(key=lambda x: x[2], reverse=True)
        # Try to assign the reads based on sorted scores
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
        Opening up space to assign unassigned reads
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
        # First pass: Attempt to assign reads based on priority
        for read_id in self.reads:
            priority = self.calculate_priority(read_id)
            self.assign_read(read_id, priority)

        # Second pass: Reprocess unassigned mappings (priority 3)
        self.reprocess_unassigned_mappings()

        # Third pass: Try to reassign unassigned reads by opening up space
        remaining_unassigned_reads = list(self.unassigned_reads)
        for read_id in remaining_unassigned_reads:
            if self.reassign_read(read_id):
                self.unassigned_reads.remove(read_id)



if __name__ == "__main__":
    samfile = sys.argv[1]
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
    with open('2bac4strain.txt', 'w') as file:
        file.write(result)
