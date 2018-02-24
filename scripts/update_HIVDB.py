import urllib.request
import re
import requests
import os

base_URL = 'https://hivdb.stanford.edu'
URL = 'https://hivdb.stanford.edu/page/algorithm-updates'
os.chdir('..')
os.chdir('./data')

# UPDATE ALGORITHM
try:
    response = urllib.request.urlopen(URL)
    html = response.read().decode('utf-8')
    xml_url = base_URL + re.search(u"/assets.*?HIVDB_.*?xml",html).group(0)
    r = requests.get(xml_url, allow_redirects=True)
    open('HIVDB.xml', 'wb').write(r.content)
    print("Updated HIVDB XML from",xml_url,"into {}\HIVDB.xml".format(os.getcwd()))
except:
    print("Unable to update HIVDB XML")

# UPDATE APOBEC DRMS
try:
    release_notes = 'https://hivdb.stanford.edu/page/release-notes/'
    response = urllib.request.urlopen(release_notes)
    html = response.read().decode('utf-8')
    apobec_url = base_URL + re.search(u"\/assets\/media\/apobec\-drms.*?tsv",html).group(0)
    r = requests.get(apobec_url, allow_redirects=True)
    open('apobec.tsv', 'wb').write(r.content)
    print("Updated APOBEC DRMs from", apobec_url, "into {}\\apobec.tsv".format(os.getcwd()))
except:
    print("Unable to update APOBEC DRMs")
