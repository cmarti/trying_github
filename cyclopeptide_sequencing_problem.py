#!/usr/bin/python
import sys

aa_masses = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L':113, 'N':114, 'D': 115,
	     'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
masses = sorted(list(set(aa_masses.values())))


def spectras_are_compatible(spectra1, spectra2):
    '''It checks if every fragment size in the spectra1 is present in spectra2'''
    for fragment in spectra1:
	if fragment not in spectra2:
	    return False
    return True


def get_theoretical_cyclic_spectra(peptide_mass_sequence):
    peptide_len = len(peptide_mass_sequence)
    cyclic_peptide = 2 * peptide_mass_sequence
    spectra = [0, sum(peptide_mass_sequence)]
    for fragment_size in range(1, peptide_len):
	for position in range(peptide_len):
	    sub_peptide = cyclic_peptide[position: position + fragment_size]
	    fragment_mass = sum(sub_peptide)
	    spectra.append(fragment_mass)
    spectra.sort()
    return spectra


def get_theoretical_linear_spectra(peptide_mass_sequence):
    peptide_len = len(peptide_mass_sequence)
    spectra = [0, sum(peptide_mass_sequence)]
    for fragment_size in range(1, peptide_len):
	for position in range(peptide_len):
	    sub_peptide = peptide_mass_sequence[position: position + fragment_size]
	    fragment_mass = sum(sub_peptide)
	    spectra.append(fragment_mass)
    spectra.sort()
    return spectra


def sequence_cyclopeptide_2(mass_spectra, candidates, present_aas):
    new_candidate_sequences = []
    for aa in present_aas:
	for candidate in candidates:
	    new_candidate = candidate + aa
	    candidate_spectra = get_theoretical_cyclic_spectra(new_candidate)
	    if spectras_are_compatible(candidate_spectra, mass_spectra):
		new_candidate_sequences.append(new_candidate)
    if new_candidate_sequences:
	new_candidate_sequences = sequence_cyclopeptide(mass_spectra)
	


def sequence_cyclopeptide(mass_spectra):
    present_aas = [aa for aa in masses if aa in mass_spectra]
    candidate_sequences = [[aa] for aa in present_aas]
    while candidate_sequences:
	new_candidate_sequences = []
	for candidate in candidate_sequences:
	    for aa in present_aas:
		new_candidate = candidate + [aa]
		candidate_spectra = get_theoretical_linear_spectra(new_candidate)
		#print new_candidate, candidate_spectra
		#raw_input()
		if spectras_are_compatible(candidate_spectra, mass_spectra):
		    new_candidate_sequences.append(new_candidate)
		    cyclic_spectra = get_theoretical_cyclic_spectra(new_candidate)
		    if cyclic_spectra == mass_spectra:
			yield '-'.join([str(x) for x in new_candidate])
	#print new_candidate_sequences
	candidate_sequences = new_candidate_sequences
    

if __name__ == '__main__':
    mass_spectra = [int(x) for x in open(sys.argv[1]).next().strip().split()]
    print ' '.join(list(sequence_cyclopeptide(mass_spectra)))
	    
    