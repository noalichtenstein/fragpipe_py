import numpy as np
import msfragger_runner
import os
import time
import shutil
import pandas as pd
import requests
from Bacteria import Bacteria
from bs4 import BeautifulSoup
import sys
import matplotlib.pyplot as plt
from comparing_distributions import compare_distributions, binary_classification, read_dna_result
import main_run_existing_fasta
import scipy.stats
import regression_pre

INPUT_NUMBER = 8
THRESHOLD = 0
BACTERIA_NAMES_FILE = "bacteria_list"


#  per raw
def get_prot_median_dict(unique_peptides_dictionary, results_folder, raw, cutoff, create_cake):
    df = pd.read_csv(f"{results_folder}\\{raw}\\peptide.tsv", sep="\t")
    bacteria = list(unique_peptides_dictionary.keys())
    proteins_intensities_dict = {}
    for bac in bacteria:
        if len(unique_peptides_dictionary[bac]) > THRESHOLD:  # number of unique peptides is above the chosen threshold
            for peptide in unique_peptides_dictionary[bac]:
                index = df[df['Peptide'] == peptide].index
                if not index.empty:
                    protein = df.loc[index, 'Protein'].iloc[0]
                    # fill proteins_intensities_dict with values
                    if protein in proteins_intensities_dict.keys():
                        proteins_intensities_dict[protein] = np.append(proteins_intensities_dict[protein],
                                                                    df.loc[index, 'Intensity'].sum()) # .sum() for extracting value from tuple
                    else:
                        proteins_intensities_dict[protein] = np.array(df.loc[index, 'Intensity'].sum()) # .sum() for extracting value from tuple

            for prot in proteins_intensities_dict.keys():
                # list_of_medians is a list of medians of specific proteins according to the order of
                # dict_protein_intensity.keys()
                proteins_intensities_dict[prot] = np.median(proteins_intensities_dict[prot])

    return proteins_intensities_dict


def normalise_intensity_values(prot_median_dict):
    values = np.array(list(prot_median_dict.values()))
    sum = np.sum(values)
    for key in prot_median_dict.keys():
        prot_median_dict[key] = prot_median_dict[key]/sum
    return prot_median_dict


def calculate_pearson(all_proteomes_folder, organisms, raw_list, cutoff_level, bacteria_list):
    """
    calculates pearson's correlations between each protein intensity and the dna relative amount
                # of dna of the bacteria contains it
    """
    # bacteria objects list according to taxonomy from the web
    bacteria_list = main_run_existing_fasta.create_bacteria_class_list(organisms)
    with open("bacteria_list", 'w') as f:
        for bac in bacteria_list:
            name = bac.genus + "_" + bac.species
            f.write(name+"\n")

    # creates dictionary {bacteria: proteins list}
    bacteria_protein_dict = main_run_existing_fasta.convert_file_to_bacteria_protein_dictionary(all_proteomes_folder)

    if cutoff_level == "species":
        proteins_dictionary_cutoff = bacteria_protein_dict
    else:
        # create dictionary according to cutoff level
        proteins_dictionary_cutoff = main_run_existing_fasta.create_cutoff_dict(bacteria_protein_dict, bacteria_list, cutoff_level)

    # creates unique bac:prot dict list for each raw
    unique_peptides_dictionary_list = main_run_existing_fasta.convert_file_to_unique_protein_dictionary(all_proteomes_folder, raw_list,
                                                                                proteins_dictionary_cutoff,
                                                                                bacteria_list, cutoff_level)
    # pearson calculation
    bacteria_names_list = list(unique_peptides_dictionary_list[0].keys())
    proteins_set = set()
    raws_dict = []

    for i, raw_name in enumerate(raw_list):
        prot_median_dict = get_prot_median_dict(unique_peptides_dictionary_list[i], all_proteomes_folder, raw_name,
                                                cutoff_level, create_cake=False)
        # change for normalising the median intensities
        prot_median_dict = normalise_intensity_values(prot_median_dict)

        raws_dict.append(prot_median_dict)


        prot_names = set(prot_median_dict.keys())
        proteins_set = proteins_set | prot_names

    proteins_list = list(proteins_set)
    pearson_correlations_per_prot = {}
    for prot in proteins_list:
        # results in order of raw
        x = []
        y = []
        for i, raw in enumerate(raw_list):
            dna_array = np.zeros(len(bacteria_list))
            if prot in raws_dict[i]:
                x.append(raws_dict[i][prot])
            else:
                x.append(0)
            sample_number = raw[raw.find('Sample') + 6]
            prenatal_num = int(raw[raw.find('Natal') + 5])
            dna_array = read_dna_result(bacteria_names_list, dna_array, sample_number, prenatal_num, "species")
            splitted = prot.split("_")
            bac_name = "_".join(splitted[1:-1])
            bacteria_index = bacteria_names_list.index(bac_name)
            y.append(dna_array[bacteria_index])

        if np.std(x) == 0 or np.std(y) == 0:
            print(f"Skipping protein {prot} due to constant input array.")
        else:
            pearson_coef = scipy.stats.pearsonr(x, y)
            pearson_correlations_per_prot[prot] = pearson_coef

        # print(f"x for prot {prot}:")
        # print(x)
        # print(f"y for prot {prot}:")
        # print(y)
        # print(f"pearson correlation for protein {prot}: {pearson_correlations_per_prot[prot]}")

    #sort the dict
    df = create_df(pearson_correlations_per_prot)


def create_df(dict):
    sorted_correlations = sorted(dict.items(), key=lambda item: item[1][0])
    protein_id, pearson_correlation, p_val = [], [], []
    for i in range(len(sorted_correlations)-1, -1, -1):
        protein_id.append(sorted_correlations[i][0])
        pearson_correlation.append(sorted_correlations[i][1][0])
        p_val.append(sorted_correlations[i][1][1])

    data = [protein_id, pearson_correlation, p_val]
    col_names = ["protein_id", "pearson_correlation", "p_val"]
    transposed_data = list(zip(*data))
    # Create a DataFrame from the sorted dictionary
    df = pd.DataFrame(transposed_data, columns=col_names)
    df.to_excel(f'pearson_correlations_ms_threshold_{THRESHOLD}_normalised.xlsx', index=False)
    return df


if __name__ == '__main__':
    """
    requires the files:
    "combined_proteome.fasta" (including human's top 100 proteins)
    "proteins_per_bacteria.txt"
    raw files
    cutoff_level:
    species/genus/family/order/class/phylum
    run msfragger value (true/false)

    input format:
    fragpipe_directory directory raw_file_list(format: 'raw1','raw2','raw3'...) organisms_list(format example:
                 'Actinomyces+naeslundii','Actinomyces+oris','Corynebacterium+amycolatum'... true/false
    """
    if len(sys.argv) != INPUT_NUMBER:
        sys.exit("Invalid input. the input should be in the format: "
                 "fragpipe_directory directory raw_file_list(format: 'raw1','raw2','raw3'...) organisms_list(format example:"
                 "'Actinomyces+naeslundii','Actinomyces+oris','Corynebacterium+amycolatum'... "
                 "species/genus/family/order/class/phylum true/false")

    fragpipe_directory = sys.argv[1]
    all_proteomes_folder = sys.argv[2]
    raw_list = sys.argv[3].split(',')
    organisms = sys.argv[4].replace("'", "").split(',')
    cutoff_level = sys.argv[5]
    fasta_name = sys.argv[6]
    run_msfragger = sys.argv[7]

    print(f"raw_list == {raw_list}")
    print(f"organisms == {organisms}")

    start_time = time.time()

    print(f"cutoff level is {cutoff_level}:")
    bacteria_names = main_run_existing_fasta.convert_file_to_list(all_proteomes_folder, BACTERIA_NAMES_FILE)
    calculate_pearson(all_proteomes_folder, organisms, raw_list, cutoff_level, bacteria_names)
    # regression_pre.create_df_of_correlative_proteins(bacteria_names,
    #                             rf"{all_proteomes_folder}\pearson_correlations_ms_threshold_0_normalised.xlsx")

    end_time = time.time()

