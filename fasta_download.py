import os
import requests
from msfragger_runner import current_date_argument


def download_proteomes_table(organism_name, fastas_number):
    """
        download proteomes table of chosen organism
        param: organism_name "Clostridium+perfringens"
    """
    response = requests.get(f"https://rest.uniprot.org/proteomes/stream?format=list&query=%28{organism_name}%29")
    proteome_list = []
    counter = 0
    if response.status_code == 200:
        for row in response.iter_lines(decode_unicode=True):
            proteome_list.append(row.strip())
            # counter += 1
            # if counter == fastas_number:
            #     break
    else:
        print("organism_name not found!")
    return proteome_list


def download_fastas(proteome_list, directory_path):
    """
        for each item in list download fasta file and creates folder with the same name as the fasta file.
        first try to download from:
        https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3A{UP_number}%29%29
        file name:
        uniprotkb_proteome_{UP_number}_{cur_date}
        * cur_date format: year_month_day

        if that doesn't work, try:
        https://rest.uniprot.org/uniparc/stream?format=fasta&query=%28%28upid%3A{UP_number}%29%29
        file name:
        uniparc_upid_{UP_number}_{cur_date}

    """
    folder_names = []
    temp_dir = f"{directory_path}\\fasta"
    create_directory(directory_path)
    cur_date = current_date_argument()
    for proteome in proteome_list:
        create_directory(directory_path)
        response1 = requests.get(f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3A"
                                 f"{proteome}%29%29")
        response2 = requests.get(f"https://rest.uniprot.org/uniparc/stream?format=fasta&query=%28%28upid%3A"
                                 f"{proteome}%29%29")

        if response1.status_code == 200 and response1.content:
            fasta_file_name = f"uniprotkb_proteome_{proteome}_{cur_date}"
            with open(f"{directory_path}\\fasta\\{fasta_file_name}.fasta", 'wb') as fasta_file:
                fasta_file.write(response1.content)
                fasta_file.close()
            os.rename(temp_dir, f"{directory_path}\\{fasta_file_name}")
            folder_names.append(fasta_file_name)

        elif response2.status_code == 200 and response2.content:
            fasta_file_name = f"uniparc_upid_{proteome}_{cur_date}"
            with open(f"{directory_path}\\fasta\\{fasta_file_name}.fasta", 'wb') as fasta_file:
                fasta_file.write(response2.content)
                fasta_file.close()
            os.rename(temp_dir, f"{directory_path}\\{fasta_file_name}")
            folder_names.append(fasta_file_name)
        else:
            os.rmdir(temp_dir)
            print(proteome, "fasta file was not found in uniprot")
    return folder_names


def create_directory(directory_path):
    try:
        os.makedirs(f"{directory_path}\\fasta")
    except OSError as e:
        print(f"Error creating folder: {str(e)}")


def download_all_fastas_of_organism(organism_name, directory_path, fastas_number):
    """
        creates a folder contains the proteomes fasta file.
        return list of names of folders (same as fasta files names)
    """
    proteome_list = download_proteomes_table(organism_name, fastas_number)
    folder_list = download_fastas(proteome_list, directory_path)
    return folder_list


def download_all_fastas_of_organism_random(organism_name, directory_path, fastas_number):
    proteome_list = download_proteomes_table(organism_name, fastas_number)
    jump = len(proteome_list) // fastas_number
    final_list = []
    if len(proteome_list) > fastas_number:
        print("jump:", jump)
        print("range(fastas_number):", range(fastas_number))
        for i in range(fastas_number):
            final_list.append(proteome_list[i*jump])
        folder_list = download_fastas(final_list, directory_path)
    else:
        folder_list = download_fastas(proteome_list, directory_path)
    return folder_list
