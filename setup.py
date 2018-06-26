from distutils.core import setup

setup(
	name = 'sierralocal',
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
	include_package_data=True
)