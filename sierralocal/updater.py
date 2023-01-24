import requests

from pathlib import Path
import re
import os

# this needs to be modified to point Python to your local chromedriver binary
mod_path = Path(os.path.dirname(__file__))


def update_APOBEC():
    """
    Update the APOBEC DRMS file from the Github page
    :return: absolute path to the apobec drms JSON file
    """
    # UPDATE APOBEC DRMS
    print('Downloading the latest APOBEC DRMS File')

    try:
        url = 'https://raw.githubusercontent.com/hivdb/hivfacts/main/data/apobecs/apobec_drms.json'
        filepath = os.path.join(mod_path, "data", "apobec_drms.json")
        request = requests.get(url, allow_redirects=True)
        with open(filepath, 'wb') as file:
            file.write(request.content)
        print("Updated APOBEC DRMs into {}".format(filepath))
    except:
        print("Unable to update APOBEC DRMs. Try manually downloading the APOBEC DRM JSON into data/apobec_drms.json")

    return filepath


def update_HIVDB():
    """
    Query the HIVdb Github page for new ASI (algorithm specification interface)
    XML files.
    :return: absolute path to new XML file
    """
    print('Downloading the latest HIVDB XML File')
    try:
        url = requests.get('https://raw.githubusercontent.com/hivdb/hivfacts/main/data/algorithms/HIVDB_latest.xml')
        file = url.text

        filepath = os.path.join(mod_path, "data", file)
        hivdb_latest = 'https://raw.githubusercontent.com/hivdb/hivfacts/main/data/algorithms/{}'.format(file)
        request = requests.get(hivdb_latest, allow_redirects=True)
        with open(filepath, 'wb') as file:
            file.write(request.content)

        print("Updated HIVDB XML into {}".format(filepath))

    except Exception as e:
        print("Unable to update HIVDB XML. Try manually downloading the HIVdb ASI2.")
        print(e)
        return None

    return filepath


def main():
    update_HIVDB()
    update_APOBEC()


if __name__ == '__main__':
    main()
