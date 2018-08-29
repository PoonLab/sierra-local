import sys
from distutils.core import setup

if sys.version_info.major < 3:
    sys.exit('Sorry, sierra-local requires Python 3.x')

setup(
	name = 'sierralocal',
    description = 'Local execution of HIVdb algorithm',
    long_description="Lightweight implementation of the Stanford HIVdb "
                     "resistance genotyping algorithm (Sierra web "
                     "service) for local execution.",

	packages = ['sierralocal'],
	version = '0.0.1',
	author = 'Jasper Ho',
	author_email = 'jasperchho@gmail.com',
	url = 'https://github.com/PoonLab/sierra-local',
	classifiers = [
		'Programming Language :: Python :: 3',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Operating System :: OS Independent'
	],

	scripts=['bin/sierralocal'],
    package_data={'sierralocal': [
        'data/*Prevalences.tsv',
        'data/*-comments.csv'
    ]},
    #include_package_data=True
    #entry_points={
    #    'console_scripts': [
    #        'sierralocal=sierralocal:main',
    #    ],
    #}
)
