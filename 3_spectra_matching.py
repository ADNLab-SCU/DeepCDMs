# -*- coding: utf-8 -*-

from rdkit import Chem
import numpy as np
import tensorflow as tf
import abc

class SimilarityProvider(object):
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, hparams=None):
        self.hparams = hparams
    
    @abc.abstractmethod
    def compute_similarity(self, library, queries):
      """Compute similarities."""

class CosineSimilarityProvider(SimilarityProvider):
    
    def _normalize_rows(self, tensor):
        return tf.nn.l2_normalize(tensor, axis=1)
    
    def compute_similarity(self, library, queries):
        similarities = tf.matmul(library, queries, transpose_b=True)
        return tf.transpose(similarities)

class GeneralizedCosineSimilarityProvider(CosineSimilarityProvider):
    """Custom cosine similarity that is popular for massspec matching."""
    
    def _make_weights(self, tensor):
      num_bins = tensor.shape[1].value
      weights = np.power(np.arange(1, num_bins + 1),
                         self.hparams.mass_power)[np.newaxis, :]
      return weights / np.sum(weights)
    
    def _normalize_rows(self, tensor):
      if self.hparams.mass_power != 0:
          tensor *= self._make_weights(tensor)
    
      return super(GeneralizedCosineSimilarityProvider,
                   self)._normalize_rows(tensor)
    
    def compute_similarity(self, library, queries):
        similarities = tf.matmul(library, queries, transpose_b=True)
        return tf.transpose(similarities)

def parse_peaks(pk_str):
    """constraction of peak_loc, peak_intnesity pair"""
    all_peaks = pk_str.split('\n')
    
    peak_locs = []
    peak_intensities = []
    
    for peak in all_peaks:
        loc, intensity = peak.split()
        peak_locs.append(int(float(loc)))
        peak_intensities.append(float(intensity))

    return peak_locs, peak_intensities

def get_parent_mass(info_str):
    """parent mass extraction"""
    lines = info_str.splitlines()
    pep_mass = float(lines[0].strip("PEPMASS="))
    
    return pep_mass

def make_dense_mass_spectra(peak_locs, peak_intensities, max_peak_loc):
    """matrix construction"""
    dense_spectrum = np.zeros(max_peak_loc)
    dense_spectrum[peak_locs] = peak_intensities
    
    return dense_spectrum

def make_query_spectra_array(ms_list):
    """convert query spectra to np.array"""
    mass_spec_spectra = np.zeros((len(ms_list), 1000))
    pep_mass_list = []
    for idx, spectra in enumerate(ms_list):
        spectra_information, spectra_str = spectra.split("min\n")
        parent_mass = get_parent_mass(spectra_information)
        spectral_locs, spectral_intensities = parse_peaks(spectra_str.rstrip("\n"))
        dense_mass_spec = make_dense_mass_spectra(
            spectral_locs, spectral_intensities, 1000)
        mass_spec_spectra[idx, :] = dense_mass_spec
        pep_mass_list.append(parent_mass)
        
    return pep_mass_list, mass_spec_spectra

def make_spectra_array(mol_list):
    """convert spectra in spectra library into np.array"""
    name_list = []
    smiles_list = []
    spectra_list = []
    exact_mass_list = []
    mass_spec_spectra = np.zeros((len(mol_list), 1000))
    for idx, mol in enumerate(mol_list):
        name_list.append(mol.GetProp("PRECURSOR_ION_NAME"))
        smiles_list.append(mol.GetProp("PRECURSOR_SMILES"))
        spectra_list.append(mol.GetProp("PREDICTED SPECTRUM"))
        exact_mass_list.append(mol.GetProp("EXACT_MASS"))
        spectral_locs, spectral_intensities = parse_peaks(mol.GetProp("PREDICTED SPECTRUM"))
        dense_mass_spec = make_dense_mass_spectra(
            spectral_locs, spectral_intensities, 1000)
        mass_spec_spectra[idx, :] = dense_mass_spec
    
    return name_list, smiles_list, mass_spec_spectra, exact_mass_list

def get_similarities(raw_spectra_library_array, raw_spectra_qurey_array, name_list):

    spec_library_array_var = tf.constant(raw_spectra_library_array)
    spec_qurey_array_var = tf.constant(raw_spectra_qurey_array)
    
    intensity_adjusted_spectra_library = tf.pow(spec_library_array_var, 0.5)
    intensity_adjusted_spectra_qurey = tf.pow(spec_qurey_array_var, 0.5)
    
    hparams = tf.contrib.training.HParams(mass_power=1.,)
    
    cos_similarity = GeneralizedCosineSimilarityProvider(hparams)
    norm_spectra_library = cos_similarity._normalize_rows(intensity_adjusted_spectra_library)
    norm_spectra_qurey = cos_similarity._normalize_rows(intensity_adjusted_spectra_qurey)
    similarity = cos_similarity.compute_similarity(norm_spectra_library, norm_spectra_qurey)
    
    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        dist = sess.run(similarity)
    
    return dist

def file_input_and_preprogressing(file_path):
    """load file, convert spectra and output parent ion list and spectra list"""
    content = []
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()
        spectras = content.split('\n\n')
    parent_mass_list, molecule_spectra_qurey_list = make_query_spectra_array(spectras)

    return parent_mass_list, molecule_spectra_qurey_list


def extract_spectra_library(spectra_library_mass):
    """extract parent mass and other infomation"""
    try:
        suppl = Chem.SDMolSupplier(spectra_library_mass)
        name_list, smiles_list, mass_spec_spectra_library, exact_mass_list = make_spectra_array(suppl)
        exact_result = True
    except:
        exact_result = False
        name_list, smiles_list, mass_spec_spectra_library, exact_mass_list = [], [], [], []
    return name_list, smiles_list, mass_spec_spectra_library, exact_result, exact_mass_list

def write_out(match_path, write_path, parent_mass, name_list, smiles_list, cos_similarity, exact_mass_list):
    """write name, smiles, cosine similarity to txt file"""
    f = open(write_path,"a+", encoding='utf-8')
    f_2 = open(match_path,"a+", encoding='utf-8')
    f.write("PEPMASS = " + parent_mass + "\n")
    match_result = []
    match_cosine = []
    for idx in range(len(name_list)):
        molecule_line = name_list[idx] + " " + smiles_list[idx] + " " + str(cos_similarity[0][idx]) + " " + str(exact_mass_list[idx]) +"\n"
        match_result.append(molecule_line)
        if cos_similarity[0][idx] >= 0.7:
            match_cosine.append(molecule_line)

    for line in match_result:
        f.write(line)
    if match_cosine != []:
        f_2.write("PEPMASS = " + parent_mass + "\n")
        for line in match_cosine:
            f_2.write(line)
        f_2.write("\n")
        f_2.close()
    f.write("\n")
    f.close()

def main():
    
    file_path = "_" # input file of spectra
    parent_mass_list, molecule_spectra_qurey_list = file_input_and_preprogressing(file_path)
    not_matched_list = []
    write_path = "_" # all match result of the input file
    match_path = "_" # match result of threshord (default 0.7)
    
    for idx, spectra in enumerate(molecule_spectra_qurey_list):
        spectra_library_mass = (
            "/spectra_library_dnscl/mass_"
            + str(int(float(parent_mass_list[idx])) - 1) + ".sdf") # self-constructed spectra library diveded based on the mass
        spectra = spectra.reshape(1, 1000)
        name_list, smiles_list, mass_spec_spectra_library, exact_result, exact_mass_list = extract_spectra_library(spectra_library_mass)
        if exact_result:
            cos_similarity = get_similarities(mass_spec_spectra_library, spectra, name_list)
            write_out(match_path, write_path, str(parent_mass_list[idx]), name_list, smiles_list, cos_similarity, exact_mass_list)
        else:
            not_matched_list.append(parent_mass_list[idx])

if __name__ == '__main__':
    main()
