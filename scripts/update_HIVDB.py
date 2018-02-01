import urllib2
import re
import requests
import os

base_URL = 'https://hivdb.stanford.edu'
URL = 'https://hivdb.stanford.edu/page/algorithm-updates'
os.chdir('..')
os.chdir('./data')

# UPDATE ALGORITHM
response = urllib2.urlopen(URL)
html = response.read()

xml_url = base_URL + re.search(u"/assets.*?HIVDB_.*?xml",html).group(0)
print(xml_url)
r = requests.get(xml_url, allow_redirects=True)

open('HIVDB.xml', 'wb').write(r.content)

# UPDATE APOBEC DRMS
release_notes = 'https://hivdb.stanford.edu/page/release-notes/'
response = urllib2.urlopen(release_notes)
html = response.read()

apobec_url = base_URL + re.search(u"\/assets\/media\/apobec\-drms.*?tsv",html).group(0)
print(apobec_url)
r = requests.get(apobec_url, allow_redirects=True)

open('apobec.tsv', 'wb').write(r.content)
