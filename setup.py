import sys
from distutils.core import setup
from setuptools.command.install import install
from distutils import log
import os
import pathlib
from setuptools import find_packages

if sys.version_info.major < 3:
    sys.exit('Sorry, sierra-local requires Python 3.x or above')

# adapted from
# https://stackoverflow.com/questions/5932804/set-file-permissions-in-setup-py-file/25761434
class OverrideInstall(install):
    def run(self):
        uid, gid = 0, 0  # root user
        mode = 0o755
        set_data_dir = False
        install.run(self)
        for filepath in self.get_outputs():
            path = os.path.dirname(filepath)
            if path.endswith('data') and not set_data_dir:
                log.info('Changing permissions of %s to %s' % (path, oct(0o777)))
                os.chmod(path, 0o777)
                set_data_dir = True


setup(
    name='sierralocal',
    description='Local execution of HIVdb algorithm',
    long_description="Lightweight implementation of the Stanford HIVdb "
                     "resistance genotyping algorithm (Sierra web "
                     "service) for local execution.",
    packages=find_packages(),
    install_requires=[
    'cython >=0.29.35',
    'post-align @ https://github.com/hivdb/post-align/archive/8e2ee118261987208c17add6cef5c1270e325a4c.zip',
    'more-itertools>=9.1.0',
    'orjson>=3.9.1',
    'types-setuptools>=67.8.0.0'
    ],
    version='0.4.1',
    author=['Jasper Ho', 'Garway Ng', 'Gopi Gugan', 'Mathias Renaud', 'William Wang', 'Art Poon'],
    author_email=['jasperchho@gmail.com', 'apoon42@uwo.ca'],
    url='https://github.com/PoonLab/sierra-local',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent'
    ],
    scripts=['bin/sierralocal'],
    options={'build_scripts': {'executable': '/usr/bin/env python3'}},
    package_data={
            'sierralocal': [
                'bin/nucamino-*',
                'data/genotype-properties.*.csv',
                'data/genotype-references.*.fasta',
                'data/*Prevalences.tsv',
                'data/*-comments.csv',
                'data/HIVDB_9.4.xml',
                'data/apobec-drms.221b0330.tsv',
                'data/alignment-config_hiv1.json',
                'data/apobec_drms.json',
                'data/apobecs.csv',
                'data/mutation-type-pairs_hiv1.csv',
                'data/sdrms_hiv1.csv',
                'data/rx-all_subtype-all.csv'
            ]},
    cmdclass={'install': OverrideInstall
              }
)
