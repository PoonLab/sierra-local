import requests

from selenium import webdriver
from selenium.webdriver.chrome.options import Options

from pathlib import Path
import re
import os


BASE_URL = 'https://hivdb.stanford.edu'

# this needs to be modified to point Python to your local chromedriver binary
mod_path = Path(os.path.dirname(__file__))
driver_path = mod_path.parent / 'bin/chromedriver'

options = Options()
options.add_argument('--headless')
browser = webdriver.Chrome(executable_path=str(driver_path), chrome_options=options)


def updateAPOBEC():
    # UPDATE APOBEC DRMS
    try:
        release_notes = 'https://hivdb.stanford.edu/page/release-notes/'
        browser.get(release_notes)
        html = browser.page_source

        apobec_url = BASE_URL + re.search(u"\/assets\/media\/apobec\-drms.*?tsv",html).group(0)
        r = requests.get(apobec_url, allow_redirects=True)
        apobec_path = mod_path/'data'/'apobec.tsv'
        open(str(apobec_path), 'wb').write(r.content)
        print("Updated APOBEC DRMs from", apobec_url, "into {}".format(apobec_path))
    except:
        print("Unable to update APOBEC DRMs. Try manually downloading the APOBEC DRM TSV into data/apobec.tsv")


def update_HIVDB():
    """
    Query the HIVdb website for new ASI (algorithm specification interface)
    XML files.
    :return: absolute path to new XML file
    """
    URL = 'https://hivdb.stanford.edu/page/algorithm-updates'
    try:
        #response = urllib.request.urlopen(URL)
        browser.get(URL)
        html = browser.page_source

        # search HTML contents for links to XML files, extract link for most recent version
        matches = re.findall(u"(/assets/media/HIVDB_[0-9a-z.]+\.xml)", html)
        versions = map(lambda x: re.search(u"_([0-9]+\.[0-9.]+)\.", x).group(1), matches)
        intermed = [(v, i) for i, v in enumerate(versions)]
        intermed.sort(reverse=True)  # descending order
        xml_url = BASE_URL + matches[intermed[0][1]]

        # download the XML file
        r = requests.get(xml_url, allow_redirects=True)
        xml_filename = mod_path/'data'/os.path.basename(xml_url)
        with open(str(xml_filename), 'wb') as file:
            file.write(r.content)

        print("Updated HIVDB XML from {} into {}".format(xml_url, xml_filename))

    except Exception as e:
        print("Unable to update HIVDB XML. Try manually downloading the HIVdb ASI2.")
        print(e)
        return None

    return(xml_filename)

def main():
    update_HIVDB()
    updateAPOBEC()

if __name__ == '__main__':
    main()
