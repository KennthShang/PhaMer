from setuptools import setup

with open('README.md') as f:
    long_description = f.read()

setup(name='phamer',
      version='0.1',
      description='Phage identification',
      url='https://github.com/KennthShang/PhaMer',
      author='SHANG Jiayu',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author_email='jyshang2-c@my.cityu.edu.hk',
      license='GPLv3',
      packages=['preprocessing.py', 'PhaMer.py'],
      package_data={'PhaMer': ['data/database.fa.bz2',
                               'data/contigs.csv',
                               'data/pc2wordsid.dict',
                               'data/pcs.csv',
                               'data/profiles.csv',
                               'data/proteins.csv',
                                  ]},
      scripts=['bin/vcontact2', 'bin/vcontact2_gene2genome'],
      install_requires=[
        'networkx>=2.2',
        'numpy>=1.20.1',
        'pandas>=1.0.5',
      ]
      )
