import subprocess
import os
import datetime
import shutil



def current_date_argument():
    # Get the current date as a datetime object
    current_date = datetime.date.today()

    # Format the current date as "year-month-day"
    formatted_date = current_date.strftime("%Y-%m-%d")
    return formatted_date


def add_decoy_and_contaminate(fragpipe_directory, fasta_folder, fasta_file):
    # List of commands to execute
    commands = [
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe'
        ' workspace --init --nocheck --temp '
        'C:\\Users\\owner\\AppData\\Local\\Temp\\be655283-c8e5-415f-8f4e-fa1187f0756a',
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe'
        f' database --custom {fasta_folder}\\{fasta_file} --contam --contamprefix',
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe'
        ' workspace --clean --nocheck'
    ]
    run_commands(commands)
    # creates:
    # work folder:
    # 2023 - 09 - 21 - decoys - contam - Klebsiella_WGLW5.fasta
    # 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_percolator_target_psms.tsv
    # 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_percolator_decoy_psms.tsv
    # 2023 - 10 - 29 - decoys - contam - Klebsiella_WGLW5.fasta


def run_ms_fragger_commands(fragpipe_directory, raw_folder, work_folder, fasta_file, raw_file):
    names = ["1fe734c7-ed97-463f-8db5-f789012e26f9", "db6a1400-18a8-4112-bcac-35f69edac709"]
    for name in names:
        folder_path = fr"C:\Users\owner\AppData\Local\Temp\{name}"
        try:
            os.makedirs(folder_path)
            print(f"Folder '{folder_path}' created successfully.")
        except OSError as e:
            print(f"Error creating folder: {str(e)}")

    cur_date = current_date_argument()

    # List of commands to execute
    commands = [
        # CheckCentroid command
        f'{fragpipe_directory}\\jre\\bin\\java.exe -Xmx20G -cp "{fragpipe_directory}\\lib\\fragpipe-20.0.jar;{fragpipe_directory}\\tools\\batmass-io-1.28.12.jar" com.dmtavt.fragpipe.util.CheckCentroid {raw_folder}\\{raw_file}.raw 18',

        # WorkspaceCleanInit command
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe workspace --clean --nocheck',

        # WorkspaceCleanInit command (again)
        # the folder C:\\Users\\owner\\AppData\\Local\\Temp\\db6a1400-18a8-4112-bcac-35f69edac709 is arbitrary
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe workspace --init --nocheck --temp C:\\Users\\owner\\AppData\\Local\\Temp\\db6a1400-18a8-4112-bcac-35f69edac709',
        # C:\Users\owner\Microbiome\MS_Fragger\fragpipe\tools\philosopher_v5.0.0_windows_amd64\philosopher.exe workspace --init --nocheck --temp C:\Users\owner\AppData\Local\Temp\abfd881e-53ec-44c7-b081-53b0ba28a842

        # MSFragger command
        # requires:
        # fragger.params
        # fragpipe.workflow
        # fragpipe-files.fp-manifest
        f'{fragpipe_directory}\\jre\\bin\\java.exe -jar -Dfile.encoding=UTF-8 -Xmx20G {fragpipe_directory}\\tools\\MSFragger-3.8\\MSFragger-3.8.jar {work_folder}\\fragger.params {raw_folder}\\{raw_file}.raw',
        # the file: 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_uncalibrated is a little different from the one created in fragpipe- there are different values in offset???
        # in work folder:
        # 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2.pin
        # 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2.pepXML
        # 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_uncalibrated

        # MSFragger move pepxml command
        f'{fragpipe_directory}\\jre\\bin\\java.exe -cp {fragpipe_directory}\\lib\\fragpipe-20.0.jar;/C:/Users/owner/Microbiome/MS_Fragger/fragpipe/lib/commons-io-2.11.0.jar com.github.chhh.utils.FileMove --no-err {raw_folder}\\{raw_file}.pepXML {work_folder}\\{raw_file}.pepXML',
        # files does not change

        # MSFragger move pin command
        f'{fragpipe_directory}\\jre\\bin\\java.exe -cp {fragpipe_directory}\\lib\\fragpipe-20.0.jar;/C:/Users/owner/Microbiome/MS_Fragger/fragpipe/lib/commons-io-2.11.0.jar com.github.chhh.utils.FileMove --no-err {raw_folder}\\{raw_file}.pin {work_folder}\\{raw_file}.pin',
        # files does not change

        # try to move 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_uncalibrated
        # todo: added now again
        f'{fragpipe_directory}\\jre\\bin\\java.exe -cp {fragpipe_directory}\\lib\\fragpipe-20.0.jar;/C:/Users/owner/Microbiome/MS_Fragger/fragpipe/lib/commons-io-2.11.0.jar com.github.chhh.utils.FileMove --no-err {raw_folder}\\{raw_file}_uncalibrated.mzML {work_folder}\\{raw_file}_uncalibrated.mzML',

        # Percolator command
        # todo: added now again
        f'{fragpipe_directory}\\tools\\percolator-306\\percolator.exe --only-psms --no-terminate --post-processing-tdc --num-threads 18 --results-psms {raw_file}_percolator_target_psms.tsv --decoy-results-psms {raw_file}_percolator_decoy_psms.tsv --protein-decoy-pattern rev_ {raw_file}.pin',
        # work folder
        # 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_percolator_decoy_psms
        # 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_percolator_target_psms

        # Percolator: Convert to pepxml command
        # todo: added now again
        # think that output path param is neseccery # old verse has C:\\Users\\owner\\Microbiome\\Analyses\\test3_ms_fragger\\2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_uncalibrated.mzML at the end
        f'{fragpipe_directory}\\jre\\bin\\java.exe -cp {fragpipe_directory}\\lib/* com.dmtavt.fragpipe.tools.percolator.PercolatorOutputToPepXML {raw_file}.pin {raw_file} {raw_file}_percolator_target_psms.tsv {raw_file}_percolator_decoy_psms.tsv interact-{raw_file} DDA 0.5 {work_folder}\\{raw_file}_uncalibrated.mzML',
        # work folder
        # interact - 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2.pep.xml

        # Percolator: Delete temp command
        # todo: added now again
        f'{fragpipe_directory}\\jre\\bin\\java.exe -cp {fragpipe_directory}\\lib\\fragpipe-20.0.jar com.github.chhh.utils.FileDelete {work_folder}\\{raw_file}_percolator_target_psms.tsv',
        # work folder
        # deletes 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_percolator_target_psms

        # Percolator: Delete temp command (again)
        # todo: added now again
        f'{fragpipe_directory}\\jre\\bin\\java.exe -cp {fragpipe_directory}\\lib\\fragpipe-20.0.jar com.github.chhh.utils.FileDelete {work_folder}\\{raw_file}_percolator_decoy_psms.tsv',
        # work folder
        # deletes 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_percolator_decoy_psms

        # todo: check new command+create temp folder?:
        # PeptideProphet: Workspace init [Work dir: C:\Users\owner\Microbiome\Analyses\fragpipe_res_8_2_24\fragpipe-2021_12_10_01_PreNatal1_Sample3_S1_p2.pepXML-temp]
        fr'{fragpipe_directory}\tools\philosopher_v5.0.0_windows_amd64\philosopher.exe workspace --init --nocheck --temp C:\Users\owner\AppData\Local\Temp\1fe734c7-ed97-463f-8db5-f789012e26f9',

        # todo: check new command:
        # PeptideProphet [Work dir: C:\Users\owner\Microbiome\Analyses\fragpipe_res_8_2_24\fragpipe-2021_12_10_01_PreNatal1_Sample3_S1_p2.pepXML-temp]
        fr'{fragpipe_directory}\tools\philosopher_v5.0.0_windows_amd64\philosopher.exe peptideprophet --decoyprobs --ppm --accmass --nonparam --expectscore --nontt --nonmc --decoy rev_ --database {work_folder}\{cur_date}-decoys-contam-combined_proteome.fasta.fas ..\{raw_file}.pepXML',

        # todo: check new command:
        # PeptideProphet: Delete temp
        fr'{fragpipe_directory}\jre\bin\java.exe -cp {fragpipe_directory}\lib\fragpipe-20.0.jar com.github.chhh.utils.FileDelete {work_folder}\fragpipe-{raw_file}.pepXML-temp',



        # todo: check new command: check if raw folder
        # TODO: error: Not all given paths exist: C:\Users\owner\Microbiome\Analyses\intensity_python_8_2_24\2021_12_10_01_PreNatal1_Sample3_S1_p2\interact-2021_12_10_01_PreNatal1_Sample3_S1_p2.pep.xml
        # Rewrite pepxml [Work dir: C:\Users\owner\Microbiome\Analyses\fragpipe_res_8_2_24]
        rf"{fragpipe_directory}\jre\bin\java.exe -cp {fragpipe_directory}\lib/* com.dmtavt.fragpipe.util.RewritePepxml {work_folder}\interact-{raw_file}.pep.xml {raw_folder}\{raw_file}\{raw_file}_uncalibrated.mzML",

        # try to move 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_uncalibrated todo: is it necessary? our command
        f'{fragpipe_directory}\\jre\\bin\\java.exe -cp {fragpipe_directory}\\lib\\fragpipe-20.0.jar;/C:/Users/owner/Microbiome/MS_Fragger/fragpipe/lib/commons-io-2.11.0.jar com.github.chhh.utils.FileMove --no-err {raw_folder}\\{raw_file}_uncalibrated.mzML {work_folder}\\{raw_file}_uncalibrated.mzML',

        # PhilosopherDbAnnotate command
        # current folder of python or res fold?
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe database --annotate {work_folder}\\{cur_date}-decoys-contam-{fasta_file}.fas --prefix rev_',
        # C:\Users\owner\AppData\Local\Temp\db6a1400-18a8-4112-bcac-35f69edac709 is arbitrary

        # PhilosopherFilter command
        # requires the file: combined.prot.xml
        # TODO: ERRO[17:24:40] Missing input parameters. missing PeptideProphet pepXML files
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe filter --psm 0.01 --tag rev_ --pepxml {work_folder}',

        # PhilosopherReport command
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe report',
        # work folder:
        # ion
        # peptide
        # protein
        # protein
        # psm

        # WorkspaceClean command
        f'{fragpipe_directory}\\tools\\philosopher_v5.0.0_windows_amd64\\philosopher.exe workspace --clean --nocheck',

        # todo: check new command:
        # IonQuant[Work dir: C:\Users\owner\Microbiome\Analyses\fragpipe_res_8_2_24]
        rf'{fragpipe_directory}\jre\bin\java.exe -Xmx20G -Dlibs.bruker.dir="{fragpipe_directory}\tools\MSFragger-3.8\ext\bruker" -Dlibs.thermo.dir="{fragpipe_directory}\tools\MSFragger-3.8\ext\thermo" -cp "{fragpipe_directory}\tools\jfreechart-1.5.3.jar;{fragpipe_directory}\tools\batmass-io-1.28.12.jar;{fragpipe_directory}\tools\IonQuant-1.9.8.jar" ionquant.IonQuant --threads 18 --ionmobility 0 --minexps 1 --mbr 0 --maxlfq 1 --requantify 1 --mztol 10 --imtol 0.05 --rttol 0.4 --mbrmincorr 0 --mbrrttol 1 --mbrimtol 0.05 --mbrtoprun 10 --ionfdr 0.01 --proteinfdr 1 --peptidefdr 1 --normalization 0 --minisotopes 2 --minscans 3 --writeindex 0 --tp 0 --minfreq 0 --minions 2 --locprob 0.75 --uniqueness 0 --filelist {work_folder}\filelist_ionquant.txt --modlist {work_folder}\modmasses_ionquant.txt'
    ]
    run_commands(commands)

def run_commands(commands):
    # Loop through and execute each command
    for command in commands:
        try:
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Print command output
            print(f'Command: {command}')
            print('Process output (stdout):')
            print(result.stdout)
            print('Process errors (stderr):')
            print(result.stderr)

            if result.returncode == 0:
                print('Command executed successfully.')
            else:
                print(f'Command failed with return code {result.returncode}.')

        except Exception as e:
            print(f'An error occurred: {str(e)}')


def move_files(work_folder, file_name):
    source_path = "C:\\Users\\owner\\Microbiome\\resources\\fragpipe_files\\" + file_name
    shutil.copy(source_path, work_folder)


"""
    This function moves relevant files for the program to the work directory and updates their contents as follows:
    path_fragger -> 'work_directory'\'cur_date'-decoys-contam-'fasta_file'.fas
    path_fragpipe -> 'work_directory'\'cur_date'-decoys-contam-'fasta_file'.fas
    path_files_fragpipe -> 'work_directory'\'raw_file'.raw
    work_dir_fragpipe -> work_directory
"""


def handle_fragpipe_files(raw_folder, work_folder, fasta_file, raw_file):
    files = ["fragger.params", "fragpipe.workflow", "fragpipe-files.fp-manifest", "filelist_ionquant.txt", "modmasses_ionquant.txt"]
    for file_name in files:
        file_path = os.path.join(work_folder, file_name)
        move_files(work_folder, file_name)

        with open(file_path, "r") as cur_file:
            file_contents = cur_file.read()

        # Replace the placeholders with the appropriate values
        file_contents = file_contents.replace("path_fragger",
                                              f'{work_folder}\\{current_date_argument()}-decoys-contam-{fasta_file}.fas')
        file_contents = file_contents.replace("path_fragpipe",
                                              f'{work_folder}\\{current_date_argument()}-decoys-contam-{fasta_file}.fas')
        file_contents = file_contents.replace("work_dir_fragpipe", work_folder)
        file_contents = file_contents.replace("work_dir", raw_folder)
        file_contents = file_contents.replace("path_files_fragpipe", f'{raw_folder}\\{raw_file}.raw')

        # Write the modified contents back to the file
        with open(file_path, "w") as cur_file:
            cur_file.write(file_contents)


def delete_files(work_folder):
    keep_files = ['ion.tsv', 'peptide.tsv', 'protein.fasta', 'protein.tsv', 'psm.tsv']
    files_to_delete = [f for f in os.listdir(work_folder) if f not in keep_files]

    for file_to_delete in files_to_delete:
        file_path = os.path.join(work_folder, file_to_delete)
        os.remove(file_path)


def run_process(inputs):
    fragpipe_directory, raw_folder, fasta_folder, work_folder, fasta_file, raw_file = inputs
    os.chdir(work_folder)
    handle_fragpipe_files(raw_folder, work_folder, fasta_file, raw_file)
    add_decoy_and_contaminate(fragpipe_directory, fasta_folder, fasta_file)
    run_ms_fragger_commands(fragpipe_directory, raw_folder, work_folder, fasta_file, raw_file)
    delete_files(work_folder)
    print("completed full run")

# work_folder:
# 2023-09-21-decoys-contam-Klebsiella_WGLW5.fasta (decoy step)
# 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2_uncalibrated.mzML
# 2023-09-21-decoys-contam-Klebsiella_WGLW5.fasta.fas.1.pepindex
# 2023-09-21-decoys-contam-Klebsiella_WGLW5.fasta.fas.2.pepindex
# whole folder - C:\Users\owner\AppData\Local\Temp\95664e5d-7aa1-44c4-934b-d437d11826f9

# work_folder:
# 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2.pepXML
# 2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2.pin
# combined.prot.xml
# filelist_proteinprophet
# filter
# fragger.params
# fragpipe.workflow
# fragpipe-files.fp-manifest
# interact-2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2.pep.xml
# ion.tsv
# log_2023-09-21_15-37-06
# peptide.tsv
# protein
# protein.tsv
# psm.tsv

# in other location AppData
# C:\Users\owner\AppData\Local\Temp\fragpipe - all folder except workflow
# C:\Users\owner\AppData\Local\Temp\jna-106164915
# C:\Users\owner\AppData\Local\Temp\hsperfdata_owner
