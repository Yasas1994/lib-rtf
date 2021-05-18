from librtf.version import __version__

from setuptools import setup, find_packages


url = 'https://github.com/Yasas1994/lib-rtf'

with open('README.md') as f:
    long_description = f.read()

setup(
    name='librtf',
    version=__version__,

    description='return time based phylogeny',
    long_description=long_description,
    long_description_content_type='text/markdown',

    url=url,
    download_url=url + '/tarball/' + __version__,
    author='Rajitha Yasas Wijesekara',
    author_email='yasas.wijesekara@uni-greifswald.de',
    keywords = 'Bioinformatics, Phylogeny, Sequence analysis',
    license='MIT',
    packages=find_packages(),

    install_requires=['numpy>=1.20.1', 'biopython>=1.78', 'dendropy>=4.5.2',
                      'scikit-learn>=0.24.1', 'pandas>=1.2.4'],
    include_package_data=True,
    python_requires='>=3.7'
)
