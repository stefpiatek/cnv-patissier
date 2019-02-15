# CNV-patissier

Orchestrates your Copy Number Variant (CNV) bakeoff, part of research project for Clinical Bioinformatics Scientist Training Programme. 

## Overview

This project is aimed to collect data on the ability of CNV callers (targeted-capture Next Generation Sequencing) to detect clinically relevant CNVs. All samples should have a gold-standard known-CNV status for the gene in question. This project is aimed at bioinformaticians working genetics laboratories but please let me know via the issues if there are any problems in installation or running, regardless. 

## Installation 

Example setup for the system and the repository is shown below, using either python's virtual environment or conda environments if you'd prefer. If you are unsure, it is probably simpler to use the python virtual environment

### System requirements

At least Python3.6 and Docker (at least engine 1.10) are required for this project, and this has only been developed for unix systems. 

1. Example Python installation for [Ubuntu](https://docs.python-guide.org/starting/install3/linux/)

    Example code at time of writing shown, please use link if this is not recent. 

    ```
    # Check python3 version
    python3 --version
    ```

    If Python version is less than 3.6

    ```
    # Add deadsnakes PPA
    sudo apt-get install software-properties-common
    sudo add-apt-repository ppa:deadsnakes/ppa
    # Update apt
    sudo apt-get update
    # install python3.6
    sudo apt-get install python3.6 
    ```

2. Example Docker installation for [Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/) with x86_64 architecture. CentOS, Debian and Fedora, along with other achitectures also available from the link.

    Example code at time of writing shown, please use link if this is not recent. 

    ```
    # Remove older versions of docker. It's fine if apt-get reports that none of the packages are installed
    sudo apt-get remove docker docker-engine docker.io containerd runc
    # Update apt
    sudo apt-get update
    # Install packages to allow apt to use repository over HTTPS
    sudo apt-get install \
        apt-transport-https \
        ca-certificates \
        curl \
    software-properties-common
    # Add Docker's official GPG key
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    # Verify key by searching for the last 8 characters of fingerprint
    sudo apt-key fingerprint 0EBFCD88
    # Add stable repository for amd64 architecture
    sudo add-apt-repository \
        "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
        $(lsb_release -cs) \
        stable"
    # Update apt
    sudo apt-get update
    # Install Docker ce
    sudo apt-get install docker-ce
    ```

### Installation choice 1: Python's virtual enviroment library, [venv](https://docs.python.org/3/library/venv.html)

1. Clone repository 

    ```
    git clone https://github.com/stefpiatek/cnv-patissier.git
    ```

2. With **Python 3.6**, set up virtual environment
    
    ```
    cd cnv-patissier
    # Create the virtual environment
    python3.6 -m venv .venv 
    # Enable the virtual environment
    source .venv/bin/activate
    ```

3. Install requirements

    ```
    pip install -r requirements.txt
    ```

### Installation choice 2: [Anaconda](https://conda.io/en/latest/index.html)

1. Clone repository 

    ```
    git clone https://github.com/stefpiatek/cnv-patissier.git
    ```

2. Create conda environment

    ```
    cd cnv-patissier
    # Create conda environment
    conda create -n cnv-patissier python=3.6 anaconda
    # Update python
    conda update python
    # Enable the cona environment
    source activate cnv-patissier
    ```

3. Install requirements

    ```
    pip install -r requirements.txt
    ```

## Setup and usage

### Overview of steps taken by CNV-patissier

- CNV-patissier is run per capture and scans the gene sample sheets and runs all CNV callers on every sample per gene. 
    - The output from each caller is created in `cnv-patissier/output/<capture>/<date-time-of-run>/<cnv-caller>/`
    - The known CNV status and called CNVs are saved in a sqlite database in `cnv-patissier/output/` 
    - A successful run settings file is written to `cnv-patissier/successful-run-settings/<capture>/<cnv-caller>/<gene>.toml`
    - Logs are written to `cnv-patissier/logs/`
- Each caller checks if there as been a successful run for that gene, if there has an no settings have changed (i.e. sample paths) then it moves onto the next. 
- If there hasn't been a successful run  the caller is run on that gene.
- If the settings have changed, the previous output will be deleted and the caller will be rerun. If you want to force a rerun, just delete the releveant successful run settings file. 

### Setup of a capture

1. Create a directory with the name of your capture in this example I will use `ICR_example`, and then created `bed` and `sample-sheets` sub-directories. If in doubt, look at the `input/ICR_example` directory in cnv-patissier.

    ```
    cd cnv-patissier
    mkdir input/ICR_example
    mkdir input/ICR_example/bed
    mkdir input/ICR_example/sample-sheets
    ```

2. Copy your *sorted* capture bed file as the <capture name>.bed e.g. `cnv-patissier/input/ICR_example/bed/ICR_example.bed`

    - The bed file should have the chromosome in the same format as you reference genome (i.e. "chr1" or "1")
    - Please follow the [BED format](http://genome.ucsc.edu/FAQ/FAQformat) with chrom, start, end and name delimited by tabs. No header.
    - The name column should be the gene name. If in doubt, please look at the example data in input for ICR_example. 


3. For each gene in the capture where you have known CNV-status using a gold-standard, create a tab-delimited sample sheet. e.g. `cnv-patissier/input/ICR_example/sample-sheets/BRCA1.txt` and `cnv-patissier/input/ICR_example/sample-sheets/BRCA2.txt`

    - Make sure to create a tab delimited file
    - The name of the file should match the name column of the bed file from step 2
    - Only create a sample sheet where you have at least 30 samples which are known not to have CNV, and have samples with known CNVs. You do not need to have a sample sheet for every gene in your capture. 
    - The column names are: 
        - sample_id: unique name for the sample
        - sample_path: full path to the bam file
        - result_type :  samples that are known to have no CNVs can be either `normal-panel` or `normal`. Samples which have a CNV are `positive`
            - There should be at least 30 `normal-panel` samples, as many `positive` samples and a similar number of `normal` samples
        - Data for positive CNV samples:
            - cnv_call: if dupliation `DUP`, if deletion `DEL`. If you really have no way of knowing, please put `unknown`
            - chromosome: check prefix matches bed file, and reference genome. Can be left blank
            - start: most 3' position of CNV from gold-standard detection, can be left blank
            - end : most 5' position of CNV from gold-standard detection, can be left blank

4. Create `settings.py` file in the base directory of cnv-patissier

    ```
    cd cnv-patissier
    cp example_settings.py settings.py
    # edit the values in the `cnv_pat_settings` dictionary of `settings.py` for your setup
    ```

### Running CNV-patissier

```
# in the root directory of cnv-patissier, with your environment activated 
# e.g.
python cnv-patissier.py ICR_example
```


## Testing

To run the tests

```
# in the root directory of cnv-patissier, with your environment activated 
python -m pytest tests/
```
