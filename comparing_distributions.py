import numpy as np
import pandas as pd
import os
import plotly.graph_objects as go
import sys
import matplotlib.pyplot as plt
import math
from scipy.stats import entropy
import seaborn as sns

INPUT_NUMBER = 4
DNA_ROW = 1
DNA_PATH_PRENATAL_1 = r"C:\Users\noa\Desktop\University\Year_3\semester_B\project\DNA_seq_results\Prenatal_1_DNA.xlsx"
DNA_PATH_PRENATAL_2 = r"C:\Users\noa\Desktop\University\Year_3\semester_B\project\DNA_seq_results\Prenatal_2_DNA.xlsx"
DNA_PATH_PRENATAL_3 = r"C:\Users\noa\Desktop\University\Year_3\semester_B\project\DNA_seq_results\Prenatal_3_DNA.xlsx"
DNA_THRESHOLD = 0


def read_dna_result(bacteria, dna_array, sample_number, prenatal_num, cutoff):
    dna_path = ''
    if prenatal_num == 1:
        dna_path = DNA_PATH_PRENATAL_1
    if prenatal_num == 2:
        dna_path = DNA_PATH_PRENATAL_2
    if prenatal_num == 3:
        dna_path = DNA_PATH_PRENATAL_3
    df_dna = pd.read_excel(dna_path)
    if cutoff == "genus":
        df_dna.iloc[:, 0] = df_dna.iloc[:, 0].str.split("|").str[0].str.split("_").str[-1]
    elif cutoff == "species":
        df_dna.iloc[:, 0] = df_dna.iloc[:, 0].str.split("_").str[-2] + "_" + df_dna.iloc[:, 0].str.split("_").str[-1]
    # g__Actinomyces|s__Actinomyces_naeslundii

    for index, row in df_dna.iterrows():
        bac = row.iloc[0]  # Extract the value in the first column
        # Convert sample_number to integer before using it as column index
        percentage = df_dna.iloc[index, int(sample_number)]
        if bac in list(bacteria):
            bec_index = list(bacteria).index(bac)
            dna_array[bec_index] += percentage

    return dna_array


def compare_distributions(bacteria, bacteria_intensities, results_folder, raw, cutoff):
    """
    gives a value for the similarity between the dna distribution and the ms distribution
    """
    dna_array = np.zeros(len(bacteria))
    sample_number = raw[raw.find('Sample')+6]
    prenatal_num = int(raw[raw.find('Natal') + 5])
    dna_array = read_dna_result(bacteria, dna_array, sample_number, prenatal_num, cutoff)/100  # /100 for normalizing to 1

    mass_vector = np.array(bacteria_intensities)
    sum_mass = np.sum(mass_vector)
    mass_vector = mass_vector/sum_mass

    epsilon = 1e-9  # Small constant to avoid division by zero or log(0)

    # Compute KL divergence with added epsilon
    kl_divergence = entropy(dna_array + epsilon, mass_vector + epsilon)

    print(f"kl_divergence for raw {raw}:", kl_divergence)
    return kl_divergence

    # print("dna_array:")
    # print(dna_array)
    # print("mass_vector:")
    # print(mass_vector)


def binary_classification(bacteria, bacteria_intensities, results_folder, raw, cutoff, dna_threshold):
    dna_array = np.zeros(len(bacteria))
    sample_number = raw[raw.find('Sample') + 6]
    prenatal_num = int(raw[raw.find('Natal') + 5])
    dna_array = read_dna_result(bacteria, dna_array, sample_number, prenatal_num, cutoff)
    mass_array = np.array(bacteria_intensities)

    return create_binary_excel_table(bacteria, mass_array, dna_array, raw, sample_number, results_folder, cutoff,
                                     dna_threshold)


def create_binary_excel_table(bacteria, mass_array, dna_array, raw_name, sample_number, results_folder, cutoff_level,
                              dna_threshold=DNA_THRESHOLD):
    # Ensure mass_array and dna_array are converted to DataFrame/Series if necessary
    mass_series = pd.Series(mass_array, index=bacteria)
    dna_series = pd.Series(dna_array, index=bacteria)

    # Create the DataFrame
    df = pd.DataFrame(index=bacteria)
    df['bacteria'] = bacteria
    df['mass spectrometry identifications'] = mass_series
    df['dna identifications'] = dna_series

    # Convert identifications to binary
    df['mass spectrometry identifications'] = (df['mass spectrometry identifications'] > 0).astype(int)
    df['dna identifications'] = (df['dna identifications'] > dna_threshold).astype(int)

    # Save to Excel
    output_filename = f"{raw_name}_{sample_number}_binary.xlsx"
    output_path = os.path.join(results_folder, output_filename)
    df.to_excel(output_path, index=False)

    # Calculate True Positives, False Negatives, False Positives, True Negatives
    TP = ((df['mass spectrometry identifications'] == 1) & (df['dna identifications'] == 1)).sum()
    FN = ((df['mass spectrometry identifications'] == 0) & (df['dna identifications'] == 1)).sum()
    FP = ((df['mass spectrometry identifications'] == 1) & (df['dna identifications'] == 0)).sum()
    TN = ((df['mass spectrometry identifications'] == 0) & (df['dna identifications'] == 0)).sum()

    return TP, FN, FP, TN

    # accuracy = (TP + TN) / (TP + TN + FP + TN)
    # precision = TP / (TP + FP)
    # sensitivity = TP / (TP + FN)

    # # Print the counts
    # print(f"Raw: {raw_name}")
    # print(f"True Positives: {TP}")
    # print(f"False Negatives: {FN}")
    # print(f"False Positives: {FP}")
    # print(f"True Negatives: {TN}")
    # print(f"accuracy: {accuracy}")
    # print(f"precision: {precision}")
    # print(f"sensitivity: {sensitivity}")


def create_heat_map(sensitivity_array, all_proteomes_folder):
    """
    creates a heat map of sensitivity to compare between different thresholds of DNA and MS
    """
    plt.figure(figsize=(12, 8))
    x_values = [round(i/3, 2) for i in range(10)]
    sns.heatmap(sensitivity_array, cmap='hot', xticklabels=x_values, yticklabels=range(10))
    plt.title(f'Sensitivity Heatmap')
    plt.xlabel('DNA Threshold')
    plt.ylabel('MS Threshold')

    # Save the figure
    output_path = os.path.join(all_proteomes_folder, f'sensitivity_heatmap.png')
    plt.savefig(output_path)
    plt.close()

