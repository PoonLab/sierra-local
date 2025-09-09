import requests
from pathlib import Path
import os
import argparse

mod_path = Path(os.path.dirname(__file__))

def update_apobec_mutation(target_dir=None):
    """
    Update the APOBEC DRMS file from the Github page
    :param updater_outdir: optional directory to save files
    :return: absolute path to the apobec drms JSON file
    """
    # UPDATE APOBEC DRMS
    print('Downloading the latest APOBEC DRMS File')

    try:
        url = 'https://raw.githubusercontent.com/hivdb/hivfacts/main/data/apobecs/apobec_drms.json'
        filepath = os.path.join(target_dir, "apobec_drms.json")
        request = requests.get(url, allow_redirects=True)
        with open(filepath, 'wb') as file:
            file.write(request.content)
        print("Updated APOBEC DRMs into {}".format(filepath))
        return filepath

    except: # pragma: no cover
        print("Unable to update APOBEC DRMs. Try manually downloading the APOBEC DRM JSON into data/apobec_drms.json")

def update_hivdb(target_dir=None):
    """
    Query the HIVdb Github page for new ASI (algorithm specification interface)
    XML files.
    :param updater_outdir: optional directory to save files
    :return: absolute path to new XML file
    """
    print('Downloading the latest HIVDB XML File')
    try:
        url = requests.get('https://raw.githubusercontent.com/hivdb/hivfacts/main/data/algorithms/HIVDB_latest.xml')
        file = url.text

        filepath = os.path.join(target_dir, file)
        hivdb_latest = 'https://raw.githubusercontent.com/hivdb/hivfacts/main/data/algorithms/{}'.format(file)
        request = requests.get(hivdb_latest, allow_redirects=True)
        with open(filepath, 'wb') as file:
            file.write(request.content)

        print("Updated HIVDB XML into {}".format(filepath))
        return filepath

    except Exception as e: # pragma: no cover
        print("Unable to update HIVDB XML. Try manually downloading the HIVdb ASI2.")
        print(e)
        return None

def update_is_unusual(target_dir=None):
    
    print('Downloading the latest file to determine is unusual')

    try:
        unusual_latest = 'https://raw.githubusercontent.com/hivdb/hivfacts/refs/heads/main/data/aapcnt/rx-all_subtype-all.csv' 
        request = requests.get(unusual_latest)
        filepath = os.path.join(target_dir, "rx-all_subtype-all.csv")
        with open(filepath, 'wb') as file:
            file.write(request.content)

        print(f'Updated is unusual file to {filepath}')
        return filepath

    except:
        print('Could not update file for is unusual (rx-all_subtype-all.csv)\n'
              'Please download manually from https://hivdb.stanford.edu/page/release-notes/#data.files')

def update_sdrms(target_dir=None):
    """
    Query the HIVDB facts github page to find and update SDRM mutations file
    @return: file path of updated file
    """
    print('Downloading the latest file to determine SDRM mutations')
    try:
        latest = 'https://raw.githubusercontent.com/hivdb/hivfacts/main/data/sdrms_hiv1.csv'
        request = requests.get(latest)
        filepath = os.path.join(target_dir, "sdrms_hiv1.csv")
        with open(filepath, 'wb') as file:
            file.write(request.content)

        print(f'Updated SDRM mutations file to {filepath}')
        return filepath

    except:
        print('Could not update file for SDRM Mutations (sdrms_hiv1.csv)\n'
              'Please download manually from https://github.com/hivdb/hivfacts/tree/main/data')

def update_mutation_type(target_dir=None):
    """
    Query the HIVDB facts github page to find and update mutations type file
    @return: file path of updated file
    """
    print('Downloading the latest file to determine mutation type')
    try:
        latest = 'https://raw.githubusercontent.com/hivdb/hivfacts/main/data/mutation-type-pairs_hiv1.csv'
        request = requests.get(latest)
        filepath = os.path.join(target_dir, "mutation-type-pairs_hiv1.csv")
        with open(filepath, 'wb') as file:
            file.write(request.content)

        print(f'Updated mutation type file to {filepath}')
        return filepath

    except:
        print('Could not update file for mutation type (mutation-type-pairs_hiv1.csv)\n'
              'Please download manually from https://github.com/hivdb/hivfacts/tree/main/data')

def update_apobec(target_dir=None):
    """
    Query the HIVDB facts github page to find and update apobec file
    @return: file path of updated file
    """
    print('Downloading the latest file to determine apobec')
    try:
        latest = 'https://raw.githubusercontent.com/hivdb/hivfacts/main/data/apobecs/apobecs.csv'
        request = requests.get(latest)
        filepath = os.path.join(target_dir, "apobecs.csv")
        with open(filepath, 'wb') as file:
            file.write(request.content)

        print(f'Updated apobecs file to {filepath}')
        return filepath

    except:
        print('Could not update file for apobecs (apobecs.csv)\n'
              'Please download manually from https://github.com/hivdb/hivfacts/tree/main/data')

def update_reference_fasta(target_dir=None):
    """
    update reference fasta file for subtyper script
    """
    print("Downloading the latest subtype reference fasta file")
    try:
        latest = "https://cms.hivdb.org/prod/downloads/hiv-genotyper/genotype-references.fasta"
        request = requests.get(latest)
        filepath = os.path.join(target_dir, 'genotype-references.fasta')
        with open(filepath, 'wb') as file:
            file.write(request.content)
        
        print(f'Updated reference fasta to {filepath}')
    except:
        print("Couldn't update subtyper reference fasta, please get manually at: https://hivdb.stanford.edu/page/hiv-subtyper/")
        
def update_genotype_properties(target_dir=None):
    """
    update genotype property file for subtyper script
    """
    print("Downloading the latest subtype genotype property File")
    try:
        latest = 'https://cms.hivdb.org/prod/downloads/hiv-genotyper/genotype-properties.tsv'
        request = requests.get(latest)
        filepath = os.path.join(target_dir, 'genotype-properties.csv')
        with open(filepath, 'wb') as file:
            file.write(request.content)
            
        print(f'Updated reference fasta to {filepath}')
    except:
        print("Couldn't update subtyper genotype property file, please get manually at: https://hivdb.stanford.edu/page/hiv-subtyper/")
    
def main(updater_outdir=None): # pragma: no cover
    update_hivdb(updater_outdir)
    update_apobec(updater_outdir)
    update_is_unusual(updater_outdir)
    update_sdrms(updater_outdir)
    update_mutation_type(updater_outdir)
    update_apobec_mutation(updater_outdir)
    update_genotype_properties(updater_outdir)
    update_reference_fasta(updater_outdir)

if __name__ == '__main__':
    # Add argument parsing for when running updater.py directly
    parser = argparse.ArgumentParser(description='Update HIVdb data files')
    parser.add_argument('-updater_outdir', default=None,
                        help='<optional> Path to folder to store updated files from updater (default: sierralocal/data folder))')
    args = parser.parse_args()

    if args.updater_outdir:
        target_dir = args.updater_outdir
    else:
        target_dir = os.path.join(mod_path, "data")
    
    # Create directory if it doesn't exist
    os.makedirs(target_dir, exist_ok=True)

    main(updater_outdir=target_dir)