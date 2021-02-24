# NPMINE 

<p align="center">
  <img src="https://github.com/computational-chemical-biology/npmine/blob/master/img/npmine.png" alt="logo"/>
</p>

NPMINE package performs semi-automated retrieval of chemical and phylogenetic information from natural products literature.

## Installation

Install conda

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh

```
   
Create a dedicated conda environment and activate

```
conda env create -f environment.yml
conda activate nplibrary
pip install git+https://gitlab.com/rsilvabioinfo/chemwalker.git
```
 
## Installation with docker


```
 docker pull ridasilva/npmine:latest
```

Alternatively one can buil the docker image

```
docker build . -t npmine 
```

If a new library needs to be installed, don't forget to update the environment.yml file 

```
conda env export | grep -v "^prefix: " > environment.yml 
```

## Using npmine docker

Go to your `npmine_library/notebooks` directory, and execute the commands

```
docker run -it -p 8888:8888 -v "$PWD":/home/npmine ridasilva/npmine:latest
jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root
```

## Getting started
 
- Load python terminal and load packages
- `Install docker and pull OSRA <https://hub.docker.com/r/cyclica/osra>`_
- `Install gnfinder <https://github.com/gnames/gnfinder>`_
- `Install oscarpdf2json command line tool <https://bitbucket.org/mjw99/chemextractor/src/master/>`_

## Usage

```
from npmine.retrieve_chemical_entities import retrieve_chemical_entities
help(retrieve_chemical_entities) 
```

## Contributing

To contribute, fork the repository and make a "Pull Request" with your edits. For major changes, please open an issue first to discuss what you would like to change

## License

This project is licensed under the BSD 3-Clause License - see the LICENSE file for details

## Acknowledgements

- `OSRA (Optical Structure Recognition Application) <https://cactus.nci.nih.gov/osra/#9>`_ - **Igor Filippov** - 2007, SAIC-Frederick, Frederick National Laboratory for Cancer Research, NIH, DHHS, Frederick, MD 

- `Image file of the latest OSRA build, based on open-babel and osra <https://hub.docker.com/r/cyclica/osra>`_ - **Leonard Morayniss**

- `gnfinder (Global Names Finder) <https://github.com/gnames/gnfinder>`_ - **gnames**

- `OSCAR4 (Open-Source Chemistry Analysis Routines) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3205045/>`_ -  **David M. Jessop, Sam E. Adams, Egon L. Willighagen, Lezan Hawizy, Peter Murray-Rust** 

- `oscarpdf2json <https://bitbucket.org/mjw99/chemextractor/src/master/>`_ - **Mark Williamson**

