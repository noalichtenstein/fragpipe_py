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
from comparing_distributions import compare_distributions, binary_classification

INPUT_NUMBER = 8
THRESHOLD = 5


def download_html(url, output_file):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad responses (4xx and 5xx status codes)

        with open(output_file, 'w', encoding='utf-8') as file:
            file.write(response.text)

        print(f"HTML content successfully downloaded and saved to {output_file}")

    except requests.RequestException as e:
        print(f"Error downloading HTML content: {e}")


def multiple_run_msfragger_runner(fragpipe_directory, all_proteomes_folder, raw_name, fasta_name):
    #  raw_folder, fasta_folder, work_folder, fasta_file, raw_file = inputs
    counter = 0
    with open(fr"{all_proteomes_folder}\completed_runs.txt", 'w') as file:
        work_folder = f'{all_proteomes_folder}\\{raw_name}'  # folder with the name of the fasta
        if not os.listdir(work_folder):
            msfragger_runner.run_process(
                [fragpipe_directory, all_proteomes_folder, all_proteomes_folder, work_folder, f'{fasta_name}.fasta',
                 raw_name])
            counter += 1
        file.write(f"{raw_name}\n")
        print(f"-----------------------------------finished {counter} runs----------------------------------------")

    file.close()


def create_raw_folder(directory_path, raw_name):
    try:
        os.makedirs(f"{directory_path}\\{raw_name}")
    except OSError as e:
        print(f"Error creating folder: {str(e)}")


def convert_list_to_file(folder, organism_name, fasta_list):
    with open(fr"{folder}\{organism_name}.txt", 'w') as file:
        for prot in fasta_list:
            file.write(f"{prot}\n")
        file.close()


def convert_file_to_list(folder, organism_name):
    fasta_list = []
    with open(rf"{folder}\{organism_name}.txt", 'r') as file:
        line = file.readline()

        while line:
            fasta_list.append(line.strip())
            line = file.readline()

        file.close()
    return fasta_list


def create_fasta_folder_from_propeties(all_proteomes_folder, organism_name):
    # to convert file_path to a fasta file names file
    path = fr'{all_proteomes_folder}\{organism_name}_properties.txt'
    fasta_list = []
    with open(f"{path}", 'r') as file:
        line = file.readline()

        while line:
            fasta_list.append(line.strip().split("\\")[6][:-1])
            line = file.readline()

        file.close()

    with open(fr"{all_proteomes_folder}\{organism_name}.txt", 'w') as file:
        for prot in fasta_list:
            file.write(f"{prot}\n")
        file.close()


def combine_human_and_bacteria(new_proteomes_folder, proteomes_fasta_list):
    for prot in proteomes_fasta_list:
        new_folder_name = fr"{new_proteomes_folder}\{prot}"
        with open(fr"{new_folder_name}\{prot}.fasta", 'a') as new_file, open(
                fr"{new_proteomes_folder}\Human_Top_100.fasta", 'r') as human:
            to_copy = human.read()
            new_file.write(to_copy)
        new_file.close()
        human.close()


def combine_human_and_bacteria_and_move_fastas(old_proteomes_folder, new_proteomes_folder, proteomes_fasta_list):
    for prot in proteomes_fasta_list:
        new_folder_name = f"{new_proteomes_folder}\\{prot}"
        os.makedirs(new_folder_name)
        shutil.copy(f"{old_proteomes_folder}\\{prot}\\{prot}.fasta", f"{new_folder_name}\\{prot}.fasta")
        with open(f"{new_folder_name}\\{prot}.fasta", 'a') as new_file, open(
                f"{new_proteomes_folder}\\Human_Top_100.fasta", 'r') as human:
            to_copy = human.read()
            new_file.write(to_copy)
        new_file.close()
        human.close()


def proteins_from_bacteria(raw_list, results_folder, bacteria_protein_dict):
    """
    for each raw file, creates a dictionary of {bacteria name: list of proteins found}
    """
    dict_list = []

    for raw in raw_list:
        df = pd.read_csv(f"{results_folder}\\{raw}\\peptide.tsv", sep="\t")
        protein_col = 'Protein'
        prots_array = df[protein_col].to_list()
        mapped_proteins_array = make_list_of_mapped_proteins(df)
        found_prots = mapped_proteins_array + prots_array
        peptide_col = df['Peptide'].values.tolist()

        bacteria_proteins_raw_dict = {key: [] for key in bacteria_protein_dict}
        # check on Protein column
        for prot in found_prots:
            for key in bacteria_protein_dict.keys():
                if prot in bacteria_protein_dict[key]:
                    bacteria_proteins_raw_dict[key].append(prot)

        dict_list.append(bacteria_proteins_raw_dict)
        create_excel_table(bacteria_proteins_raw_dict, raw, results_folder, "", cutoff_level)
        # print(bacteria_proteins_raw_dict)

    return dict_list


def proteins_from_bacteria_as_peptides(raw_list, results_folder, bacteria_protein_dict):
    """
    for each raw file, creates a dictionary of {bacteria name: list of proteins found}
    """
    dict_list = []

    for raw in raw_list:
        df = pd.read_csv(f"{results_folder}\\{raw}\\peptide.tsv", sep="\t")
        protein_col = 'Protein'
        prots_array = df[protein_col].to_list()
        mapped_proteins_array = make_list_of_mapped_proteins(df)
        peptide_col = df['Peptide'].values.tolist()

        bacteria_proteins_raw_dict = {key: [] for key in bacteria_protein_dict}
        for i, peptide in enumerate(peptide_col):
            for key in bacteria_protein_dict.keys():
                array_of_mapped_proteins = get_mapped_proteins_of_peptide(i, df)
                prots_containing_peptide = array_of_mapped_proteins + [prots_array[i]]
                for prot in prots_containing_peptide:
                    if prot in bacteria_protein_dict[key]:
                        bacteria_proteins_raw_dict[key].append(peptide)

        dict_list.append(bacteria_proteins_raw_dict)
        # create_excel_table(bacteria_proteins_raw_dict, raw, results_folder, "", cutoff_level)
        # print(bacteria_proteins_raw_dict)

    return dict_list


def create_excel_table(bacteria_proteins_raw_dict, raw_name, results_folder, is_unique, cutoff_level):
    for bac in bacteria_proteins_raw_dict.keys():
        bacteria_proteins_raw_dict[bac].insert(0, len(bacteria_proteins_raw_dict[bac]))

    df = pd.DataFrame.from_dict(bacteria_proteins_raw_dict, orient='index')
    # Transpose the DataFrame to have names as columns
    df = df.transpose()
    df.loc['Count'] = df.count(axis=1)
    # Create an Excel writer
    excel_writer = pd.ExcelWriter(f'{results_folder}\\output_raw_{raw_name}_{is_unique}_by_{cutoff_level}.xlsx',
                                  engine='xlsxwriter')

    # Write the DataFrame to a sheet in the Excel file
    df.to_excel(excel_writer, sheet_name='Sheet1', index=False)
    # here: add blank row between rows 1 and 2
    # in each col write here the len of the col's values

    # Save the Excel file
    excel_writer.close()
    for bac in bacteria_proteins_raw_dict.keys():
        bacteria_proteins_raw_dict[bac].pop(0)


def get_mapped_proteins_of_peptide(i, df):
    mapped_proteins_df = df['Mapped Proteins'].values
    mapped_proteins_array = []
    cell = mapped_proteins_df[i]
    if type(cell) == str:
        tab = cell.split(',')
        for prot in tab:
            mapped_proteins_array.append(prot.strip())

    return mapped_proteins_array


def make_list_of_mapped_proteins(df):
    mapped_proteins_df = df['Mapped Proteins'].values
    mapped_proteins_array = []
    for tab in mapped_proteins_df:
        if type(tab) == str:
            tab = tab.split(',')
            for prot in tab:
                mapped_proteins_array.append(prot.strip())

    return mapped_proteins_array


def convert_file_to_bacteria_protein_dictionary(results_folder):
    """
    creates dictionary {bacteria name: proteins found}
    from the file: rf"{results_folder}\proteins_per_bacteria.txt"
    """
    bacteria_proteins_dictionary = {}
    with open(rf"{results_folder}\proteins_per_bacteria.txt") as file:
        cur_bacteria = ""
        for line in file:
            if line[0] == '*':
                cur_bacteria = line[1:-1]
                bacteria_proteins_dictionary[cur_bacteria] = []
            else:
                bacteria_proteins_dictionary[cur_bacteria].append(line[:-1])

        file.close()
    return bacteria_proteins_dictionary


def get_taxonomy_level_cutoff(bacteria, cutoff):
    taxonomy_level_name = ""
    if cutoff == "phylum":
        taxonomy_level_name = bacteria.phylum

    elif cutoff == "class":
        taxonomy_level_name = bacteria.bacteria_class

    elif cutoff == "order":
        taxonomy_level_name = bacteria.order

    elif cutoff == "family":
        taxonomy_level_name = bacteria.family

    elif cutoff == "genus":
        taxonomy_level_name = bacteria.genus

    elif cutoff == "species":
        taxonomy_level_name = bacteria.genus + "_" + bacteria.species

    return taxonomy_level_name


def convert_file_to_unique_protein_dictionary(results_folder, raw_list, cutoff_proteins_dictionary, bacteria_list,
                                              cutoff_level):
    """
    creates dictionary {bacteria name: unique proteins found}
    from the file: rf"{results_folder}\proteins_per_bacteria.txt"
    """
    dict_list = []
    for raw in raw_list:
        unique_raw_dict = {}
        for group in cutoff_proteins_dictionary.keys():
            unique_raw_dict[group] = []

        df = pd.read_csv(f"{results_folder}\\{raw}\\peptide.tsv", sep="\t")
        protein_col = df['Protein'].values.tolist()
        mapped_proteins_array = df['Mapped Proteins'].values.tolist()
        peptide_col = df['Peptide'].values.tolist()

        for i, protein in enumerate(protein_col):
            appended_protein_as_list = protein.split('_')
            if appended_protein_as_list[0] == "bacteria":
                # check if are_mapped_proteins_contaminant
                are_mapped_proteins_contaminant = True
                if isinstance(mapped_proteins_array[i], str) and mapped_proteins_array[i].lower() != 'nan':
                    mapped_proteins_as_list = mapped_proteins_array[i].split(',')
                    for mapped_protein in mapped_proteins_as_list:
                        # if are_mapped_proteins_contaminant == True protein is unique
                        are_mapped_proteins_contaminant = are_mapped_proteins_contaminant and mapped_protein.split('_')[0].strip() != "bacteria"

                if are_mapped_proteins_contaminant:
                    for bacteria in bacteria_list:
                        if bacteria.species == appended_protein_as_list[2] and bacteria.genus == \
                                appended_protein_as_list[1]:
                            taxonomy_level = get_taxonomy_level_cutoff(bacteria, cutoff_level)
                            unique_raw_dict[taxonomy_level].append(peptide_col[i])

                # what is the taxonomy level and name of protein.
                # check all the proteins under this level. if a mapped protein is not in this group, not unique
                else:
                    taxonomy_level = ''
                    for bacteria in bacteria_list:
                        species_name = appended_protein_as_list[2]
                        if len(appended_protein_as_list) == 5:
                            species_name = appended_protein_as_list[2] + '_' + appended_protein_as_list[3]
                        if bacteria.species == species_name:
                            taxonomy_level = get_taxonomy_level_cutoff(bacteria, cutoff_level)
                            break

                    unique = True
                    mapped_proteins_as_list = mapped_proteins_array[i].split(',')
                    for mapped_protein in mapped_proteins_as_list:
                        if mapped_protein.split('_')[0].strip() == "bacteria":
                            if mapped_protein in cutoff_proteins_dictionary[taxonomy_level]:
                                continue
                            else:
                                unique = False
                                break
                    if unique and taxonomy_level:
                        unique_raw_dict[taxonomy_level].append(peptide_col[i])

        dict_list.append(unique_raw_dict)
    return dict_list


def create_lengths_dicts(bacteria_prots_dict_list):
    len_dict_list = []
    for dict in bacteria_prots_dict_list:
        len_dict = {}
        for key in dict.keys():
            len_dict[key] = len(dict[key])

        len_dict_list.append(len_dict)
    return len_dict_list


def create_bacteria_class_list(organisms):
    list_of_bacterias = []
    for organism in organisms:
        url = fr"https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism_id&format=tsv&query=%28{organism}%29"
        response = requests.get(url)
        if response.status_code == 200:
            taxonomy_complete = False
            # Initialize variables
            phylum = None
            bacteria_class = None
            order = None
            family = None
            loop_number = 1
            organism_ids = response.content.decode()
            while not taxonomy_complete:
                current_row = organism_ids.splitlines()[loop_number]
                if not current_row:
                    Bacteria(phylum, bacteria_class, order, family, organism.split('+')[0], organism.split('+')[1])
                    print(f"bacteria: {organism.split('+')[0]} {organism.split('+')[1]} taxonomy not complete")
                    break
                organism_id = int(current_row.split('\t')[1].strip())

                url = fr"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id={organism_id}&lvl=3&lin=s&keep=1&srchmode=1&unlock&log_op=lineage_toggle"

                try:
                    response = requests.get(url)
                    response.raise_for_status()  # Raise an exception for bad responses (4xx and 5xx status codes)

                    html_content = response.text  # 'html_content' contains the HTML of the web page
                    soup = BeautifulSoup(html_content, 'html.parser')  # Create a BeautifulSoup object

                    # Extracting information
                    if soup.find('a', {'title': 'superkingdom'}).text == "Viruses":
                        loop_number += 1
                        continue
                    if not phylum and soup.find('a', {'title': 'phylum'}):
                        phylum = soup.find('a', {'title': 'phylum'}).text
                    if not bacteria_class and soup.find('a', {'title': 'class'}):
                        bacteria_class = soup.find('a', {'title': 'class'}).text
                    if not order and soup.find('a', {'title': 'order'}):
                        order = soup.find('a', {'title': 'order'}).text
                    if not family and soup.find('a', {'title': 'family'}):
                        family = soup.find('a', {'title': 'family'}).text

                    # Check if all variables are found
                    if phylum is not None and bacteria_class is not None and order is not None and family is not None:
                        taxonomy_complete = True
                        list_of_bacterias.append(Bacteria(phylum, bacteria_class, order, family, organism.split('+')[0],
                                                          organism.split('+')[1]))
                    else:
                        loop_number += 1

                except requests.RequestException as e:
                    print(f"Error downloading HTML content: {e}")

    return list_of_bacterias


def create_cutoff_dict(bacteria_protein_dict, bacteria_list, cutoff_level):
    """
    dict from the example:
     bacteria.phylum ; [list of the mapped proteins]
    """
    dict = {}
    for bacteria in bacteria_list:
        if cutoff_level == "phylum":
            if bacteria.phylum in dict.keys():
                dict[bacteria.phylum].extend(bacteria_protein_dict[bacteria.genus + "_" + bacteria.species])
            else:
                dict[bacteria.phylum] = bacteria_protein_dict[bacteria.genus + "_" + bacteria.species]

        elif cutoff_level == "class":
            if bacteria.bacteria_class in dict.keys():
                dict[bacteria.bacteria_class].extend(bacteria_protein_dict[bacteria.genus + "_" + bacteria.species])
            else:
                dict[bacteria.bacteria_class] = bacteria_protein_dict[bacteria.genus + "_" + bacteria.species]

        elif cutoff_level == "order":
            if bacteria.order in dict.keys():
                dict[bacteria.order].extend(bacteria_protein_dict[bacteria.genus + "_" + bacteria.species])
            else:
                dict[bacteria.order] = bacteria_protein_dict[bacteria.genus + "_" + bacteria.species]

        elif cutoff_level == "family":
            if bacteria.family in dict.keys():
                dict[bacteria.family].extend(bacteria_protein_dict[bacteria.genus + "_" + bacteria.species])
            else:
                dict[bacteria.family] = bacteria_protein_dict[bacteria.genus + "_" + bacteria.species]

        elif cutoff_level == "genus":
            if bacteria.genus in dict.keys():
                dict[bacteria.genus].extend(bacteria_protein_dict[bacteria.genus + "_" + bacteria.species])
            else:
                dict[bacteria.genus] = bacteria_protein_dict[bacteria.genus + "_" + bacteria.species]

    return dict


def calculate_distribution(unique_peptides_dictionary, results_folder, raw):
    df = pd.read_csv(f"{results_folder}\\{raw}\\peptide.tsv", sep="\t")
    bacteria = unique_peptides_dictionary.keys()
    bacteria_intensities = []
    for bac in bacteria:
        intensities_sum = 0
        for peptide in unique_peptides_dictionary[bac]:
            index = df[df['Peptide'] == peptide].index
            if not index.empty:
                intensities_sum += df.loc[index, 'Intensity'].sum()
        bacteria_intensities.append(intensities_sum)

    create_cake_plot_intensities(bacteria, bacteria_intensities, results_folder, raw)


def calculate_distribution_above_threshold(unique_peptides_dictionary, results_folder, raw, cutoff):
    """
        for bacteria with unique number of peptides larger than threshold:
         calculate distribution according to intensity sum
    """
    df = pd.read_csv(f"{results_folder}\\{raw}\\peptide.tsv", sep="\t")
    bacteria = unique_peptides_dictionary.keys()
    bacteria_intensities = []
    for bac in bacteria:
        if len(unique_peptides_dictionary[bac]) > THRESHOLD:
            intensities_sum = 0
            for peptide in unique_peptides_dictionary[bac]:
                index = df[df['Peptide'] == peptide].index
                if not index.empty:
                    intensities_sum += df.loc[index, 'Intensity'].sum()
            bacteria_intensities.append(intensities_sum)
        else:
            bacteria_intensities.append(0)

    # create_cake_plot_intensities(bacteria, bacteria_intensities, results_folder, raw)
    compare_distributions(bacteria, bacteria_intensities, results_folder, raw, cutoff)


def calculate_distribution_intensity_median(unique_peptides_dictionary, results_folder, raw, cutoff, create_cake):
    """
        for bacteria with unique number of peptides larger than threshold:
         calculate distribution according to intensity median per protein
    """
    df = pd.read_csv(f"{results_folder}\\{raw}\\peptide.tsv", sep="\t")
    bacteria = list(unique_peptides_dictionary.keys())
    bacteria_intensities = []
    for bac in bacteria:
        dict_protein_intensity = {}
        if len(unique_peptides_dictionary[bac]) > THRESHOLD:  # number of unique peptides is above the chosen threshold
            for peptide in unique_peptides_dictionary[bac]:
                index = df[df['Peptide'] == peptide].index
                if not index.empty:
                    protein = df.loc[index, 'Protein'].iloc[0]
                    if protein in dict_protein_intensity.keys():
                        dict_protein_intensity[protein] = np.append(dict_protein_intensity[protein], df.loc[index, 'Intensity'].sum())
                    else:
                        dict_protein_intensity[protein] = np.array(df.loc[index, 'Intensity'].sum())

            list_of_medians = np.zeros(0)
            for prot in dict_protein_intensity.keys():
                # list_of_medians is a list of medians of specific proteins according to the order of dict_protein_intensity.keys()
                list_of_medians = np.append(list_of_medians, np.median(dict_protein_intensity[prot]))

            # sum_list_of_medians = np.sum(list_of_medians)
            # bacteria_intensities.append(sum_list_of_medians)
            # mean
            proteins_mean = np.mean(list_of_medians)
            bacteria_intensities.append(proteins_mean)

        else:
            bacteria_intensities.append(0)

    if create_cake:
        create_cake_plot_intensities(bacteria, bacteria_intensities, results_folder, raw)

    # compare_distributions(bacteria, bacteria_intensities, results_folder, raw, cutoff)
    return binary_classification(bacteria, bacteria_intensities, results_folder, raw, cutoff)



def create_cake_plot_intensities(bacteria, bacteria_intensities, results_folder, raw):
    # Filter bacteria and intensities with intensities greater than 0
    filtered_data = [(b, i) for b, i in zip(bacteria, bacteria_intensities) if i > 0]
    sorted_data = sorted(filtered_data, key=lambda x: x[1], reverse=True)
    sorted_bacteria, sorted_intensities = zip(*sorted_data)

    colors = plt.cm.tab20c.colors
    # Define fontsize
    fontsize = 12

    # Create the plot
    plt.figure(figsize=(20, 20))
    plt.pie(sorted_intensities, labels=sorted_bacteria, autopct='%1.1f%%', startangle=140, colors=colors,
            textprops={'fontsize': fontsize})

    # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.axis('equal')

    # Title for the pie chart
    plt.title(f'Bacteria Intensity Distribution {raw[14:-3]}', fontsize=20, y=1.05)
    plt.savefig(fr"{results_folder}\{raw}_cake_plt_intensities.png", format='png', dpi=300)
    plt.show()


def sort_peptides(all_proteomes_folder, organisms, raw_list, cutoff_level):
    bacteria_list = create_bacteria_class_list(organisms)  # bacteria objects list according to taxonomy from the web

    # creates dictionary {bacteria: proteins list}
    bacteria_protein_dict = convert_file_to_bacteria_protein_dictionary(all_proteomes_folder)

    if cutoff_level == "species":
        proteins_dictionary_cutoff = bacteria_protein_dict
    else:
        # create dictionary according to cutoff level
        proteins_dictionary_cutoff = create_cutoff_dict(bacteria_protein_dict, bacteria_list, cutoff_level)

    # creates unique bac:prot dict list for each raw
    unique_peptides_dictionary_list = convert_file_to_unique_protein_dictionary(all_proteomes_folder, raw_list,
                                                                                proteins_dictionary_cutoff,
                                                                                bacteria_list, cutoff_level)
    # creates not unique excel for each raw
    bacteria_peptides_dict_list = proteins_from_bacteria_as_peptides(raw_list, all_proteomes_folder,
                                                                     proteins_dictionary_cutoff)

    sum_TP, sum_FN, sum_FP, sum_TN = 0, 0, 0, 0
    #  creates excel table of unique proteins for each raw
    for i, raw_name in enumerate(raw_list):
        # create_excel_table(unique_peptides_dictionary_list[i], raw_name, all_proteomes_folder, "unique", cutoff_level)
        # calculate_distribution(unique_peptides_dictionary_list[i], all_proteomes_folder, raw_name)
        # # intensity sum calculation

        # add for binary classification
        calculate_distribution_above_threshold(unique_peptides_dictionary_list[i], all_proteomes_folder, raw_name, cutoff_level)
        TP, FN, FP, TN = calculate_distribution_intensity_median(unique_peptides_dictionary_list[i], all_proteomes_folder, raw_name, cutoff_level,create_cake=False)
        sum_TP += TP
        sum_FN += FN
        sum_FP += FP
        sum_TN += TN

    accuracy = (sum_TP + sum_TN) / (sum_TP + sum_TN + sum_FP + sum_TN)
    precision = sum_TP / (sum_TP + sum_FP)
    sensitivity = sum_TP / (sum_TP + sum_FN)

    # Print the counts
    print(f"True Positives: {sum_TP}")
    print(f"False Negatives: {sum_FN}")
    print(f"False Positives: {sum_FP}")
    print(f"True Negatives: {sum_TN}")
    print(f"accuracy: {accuracy}")
    print(f"precision: {precision}")
    print(f"sensitivity: {sensitivity}")

    # lengths_unique_dict = create_lengths_dicts(unique_peptides_dictionary_list)
    # lengths_dict = create_lengths_dicts(bacteria_peptides_dict_list)
    # print()
    # for dict in lengths_dict:
    #     print(dict)
    #
    # print("\nunique:")
    # for unique_dict in lengths_unique_dict:
    #     print(unique_dict)


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

    if run_msfragger == "true":
        # creating raw folders in fasta folder for each raw file
        for raw_name in raw_list:
            create_raw_folder(f"{all_proteomes_folder}", raw_name)

        # running of raw_list in each fasta folder
        for raw_name in raw_list:
            multiple_run_msfragger_runner(fragpipe_directory, all_proteomes_folder, raw_name, fasta_name)

        end_time = time.time()
        running_time = end_time - start_time
        print(f"total running time: {running_time} seconds")

    print(f"cutoff level is {cutoff_level}:")
    sort_peptides(all_proteomes_folder, organisms, raw_list, cutoff_level)

# input example
#     all_proteomes_folder = r"C:\Users\owner\Microbiome\Analyses\intensity_7_3_24"
#     organisms = [
#         'Actinomyces+naeslundii',
#         'Enterobacter+cloacae_complex',
#         'Actinomyces+oris',
#         'Corynebacterium+amycolatum',
#         'Corynebacterium+striatum',
#         'Dermabacter+hominis',
#         'Rothia+mucilaginosa',
#         'Cutibacterium+avidum',
#         'Gemella+haemolysans',
#         'Staphylococcus+epidermidis',
#         'Staphylococcus+hominis',
#         'Enterococcus+faecalis',
#         'Lactobacillus+rhamnosus',
#         'Streptococcus+infantis',
#         'Streptococcus+salivarius',
#         'Clostridium+perfringens',
#         'Negativicoccus+succinicivorans',
#         'Veillonella+atypica',
#         'Veillonella+parvula',
#         'Finegoldia+magna',
#         'Citrobacter+braakii',
#         'Citrobacter+youngae',
#         'Escherichia+coli',
#         'Klebsiella+aerogenes',
#         'Klebsiella+michiganensis',
#         'Klebsiella+oxytoca',
#         'Klebsiella+pneumoniae',
#         'Kluyvera+ascorbata',
#         'Proteus+mirabilis'
#     ]
#  raw_list = [
#         # "2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2",
#         # "2022_05_25_07_PreNatal3_Sample2_S1_SCX_p2",
#         # "2022_05_25_05_PreNatal3_Sample3_S1_SCX_p2",
#         # "2022_05_25_03_PreNatal3_Sample4_S1_SCX_p2",
#         # "2022_05_25_01_PreNatal3_Sample5_S1_SCX_p2",
#         "2021_12_10_11_PreNatal2_Sample2_S1_p2",
#         "2021_12_10_13_PreNatal2_Sample3_S1_p2",
#         "2021_12_10_07_PreNatal1_Sample1_S1_p2",
#         "2021_12_10_09_PreNatal1_Sample2_S1_p2",
#         "2021_12_10_01_PreNatal1_Sample3_S1_p2",
#         "2021_12_10_02_PreNatal2_Sample1_S1_p2"
#     ]
#     cutoff_level = "species"
