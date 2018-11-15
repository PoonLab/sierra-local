import sys
from distutils.core import setup
#from setuptools import setup
from setuptools.command.install import install
from distutils import log
import os

if sys.version_info.major < 3:
    sys.exit('Sorry, sierra-local requires Python 3.x')

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
	name = 'sierralocal',
    description = 'Local execution of HIVdb algorithm',
    long_description="Lightweight implementation of the Stanford HIVdb "
                     "resistance genotyping algorithm (Sierra web "
                     "service) for local execution.",

	packages = ['sierralocal'],
	version = '0.1.2',
	author = 'Jasper Ho',
	author_email = 'jasperchho@gmail.com',
	url = 'https://github.com/PoonLab/sierra-local',
	classifiers = [
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
            'data/*-comments.csv'
    ]},
    cmdclass={'install': OverrideInstall},
    #include_package_data=True
    #entry_points={
    #    'console_scripts': [
    #        'sierralocal=sierralocal:main',
    #    ],
    #}
)
