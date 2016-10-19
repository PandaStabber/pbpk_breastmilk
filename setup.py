from setuptools import setup, find_packages

setup(
    name='pbpk_breastmilk',
    version='0.1',
    author='Eli Goldberg and Tenzing Gyalpo',
    author_email='elisgoldberg@gmail.com',
    packages=find_packages(),
    description='A generic, one-box, congener transport model describing mother to child pollutant transport via '
                'breastmilk',
    long_description=open('README.MD').read(),
    install_requires=[
        'args==0.1.0',
        'clint==0.5.1',
        'cycler==0.10.0',
        'matplotlib==1.5.3',
        'numpy==1.11.0',
        'pandas==0.18.1',
        'pkginfo==1.3.2',
        'pyparsing==2.1.10',
        'python-dateutil==2.5.3',
        'pytz==2016.7',
        'requests==2.11.1',
        'requests-toolbelt==0.7.0',
        'scipy==0.17.1',
        'six==1.10.0',
        'twine==1.8.1',
    ])
