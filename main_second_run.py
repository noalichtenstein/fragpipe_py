import fasta_download
import msfragger_runner
import os
import time
import shutil


def multiple_run_msfragger_runner(all_proteomes_folder, proteomes_fasta_list, raw_name):
    """
        return list of folder names of all completed runs
    """
    completed_analysis = []
    #  raw_folder, fasta_folder, work_folder, fasta_file, raw_file = inputs
    for prot in proteomes_fasta_list:
        fasta_folder = f'{all_proteomes_folder}\\{prot}'
        work_folder = f'{all_proteomes_folder}\\{prot}\\{raw_name}'  # folder with the name of the fasta
        msfragger_runner.run_process([all_proteomes_folder, fasta_folder, work_folder, f'{prot}.fasta', raw_name])
        completed_analysis.append(prot)

    return completed_analysis


def multiple_run_msfragger_runner_if_not_empty(all_proteomes_folder, proteomes_fasta_list, raw_name, organism_name):
    completed_analysis = []
    #  raw_folder, fasta_folder, work_folder, fasta_file, raw_file = inputs
    counter = 0
    with open(f"{all_proteomes_folder}\\{organism_name}_completed_runs.txt", 'w') as file:
        for prot in proteomes_fasta_list:
            fasta_folder = f'{all_proteomes_folder}\\{prot}'
            work_folder = f'{all_proteomes_folder}\\{prot}\\{raw_name}'  # folder with the name of the fasta
            if not os.listdir(work_folder):
                msfragger_runner.run_process(
                    [all_proteomes_folder, fasta_folder, work_folder, f'{prot}.fasta', raw_name])
                counter += 1
            completed_analysis.append(prot)
            file.write(f"{prot}_{raw_name}\n")
            print(f"-----------------------------------finished {counter} runs----------------------------------------")

    file.close()


def create_raw_folder(directory_path, raw_name):
    try:
        os.makedirs(f"{directory_path}\\{raw_name}")
    except OSError as e:
        print(f"Error creating folder: {str(e)}")


def convert_list_to_file(folder, organism_name, fasta_list):
    with open(f"{folder}\\{organism_name}.txt", 'w') as file:
        for prot in fasta_list:
            file.write(f"{prot}\n")
        file.close()


def convert_file_to_list(folder, organism_name):
    fasta_list = []
    with open(f"{folder}\\{organism_name}.txt", 'r') as file:
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
        new_folder_name = f"{new_proteomes_folder}\\{prot}"
        with open(f"{new_folder_name}\\{prot}.fasta", 'a') as new_file, open(
                f"{new_proteomes_folder}\\Human_Top_100.fasta", 'r') as human:
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


if __name__ == '__main__':
    start_time = time.time()
    all_proteomes_folder = r"C:\Users\owner\Microbiome\Analyses\Clostridium_perfringens_human_ms_raw2"
    old_proteomes_folder = r"C:\Users\owner\Microbiome\Analyses\Clostridium_perfringens_ms"
    raw_list = [
        # "2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2",
        "2022_05_25_07_PreNatal3_Sample2_S1_SCX_p2",
        # "2022_05_25_05_PreNatal3_Sample3_S1_SCX_p2",
        # "2022_05_25_03_PreNatal3_Sample4_S1_SCX_p2",
        # "2022_05_25_01_PreNatal3_Sample5_S1_SCX_p2"
    ]
    organism_name = "Clostridium+perfringens"

    if old_proteomes_folder:
        proteomes_fasta_list = convert_file_to_list(old_proteomes_folder, organism_name)
        combine_human_and_bacteria_and_move_fastas(old_proteomes_folder, all_proteomes_folder, proteomes_fasta_list)
    else:
        proteomes_fasta_list = fasta_download.download_all_fastas_of_organism(organism_name, all_proteomes_folder, 250)
        combine_human_and_bacteria(all_proteomes_folder, proteomes_fasta_list)

    # creating raw folders in fasta folder for each raw file
    for proteom in proteomes_fasta_list:
        for raw_name in raw_list:
            create_raw_folder(f"{all_proteomes_folder}\\{proteom}", raw_name)

    # running of raw_list in each fasta folder
    for raw_name in raw_list:
        completed_analysis = multiple_run_msfragger_runner_if_not_empty(all_proteomes_folder, proteomes_fasta_list,
                                                                        raw_name, organism_name)
        print(completed_analysis)

    end_time = time.time()
    running_time = end_time - start_time
    print(f"total running time: {running_time} seconds")

    # create_fasta_folder_from_propeties(all_proteomes_folder, organism_name)
