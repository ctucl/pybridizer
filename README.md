
![pybridizer logo](https://github.com/WhoIsJack/pybridizer/raw/main/_media/pybridizer_logo.png#gh-light-mode-only)
![pybridizer logo](https://github.com/WhoIsJack/pybridizer/raw/main/_media/pybridizer_logo_lite.png#gh-dark-mode-only)


Pybridizer is a python repository to generate HCR v3.0 probes for in situ hybridization of mRNA. 

The module allows for quick and easy design of probe pairs for the Hybridization Chain Reaction approach (Choi et al. Development 2018.)


### Installation: BLAST+ (Pre-requisite)

- Install BLAST+ by downloading executable from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    - Windows users: download the file ending with `win64.exe` and follow the instructions
    - Windows users: as of early 2026, BLAST+ **cannot deal with whitespace in paths!**
        - Install in a folder without blank spaces, e.g. `C:\Users\yourname\NCBI`

* Windows users: Add required environment variables
    - Open “Edit environment variables for your account”
    - Ensure the install added the `...\NCBI\blast-2.17.0+\bin` (or similar) path to the “Path” variable
    - Add variable `BLASTDB` and set it to track a `db` folder like `...\NCBI\blast-2.17.0+\db`
    - Add variable `BLASTDB_LMDB_MAP_SIZE` and set it to `1000000`

- Download refseq data for the desired model organism
    - Search NCBI for `danio rerio reference genome` (or equivalent)
    - Click `Download`, select `Refseq only`, tick only the `Transcripts (FASTA)` sequences
    - Unzip/extract the `rna.fna` data file that should be inside a subfolder like `ncbi_dataset\data\GCF_000002035.6`
    - Ensure the extracted file is located in a folder without spaces, e.g. `C:\Users\yourname\NCBI\fasta_src\`
        - Also give the file a clear name with all relevant information, e.g. `drerio_ncbi-GRCz12tu_refseq_transcripts.fna`

* Open cmd and run: `makeblastdb -in <fna-file-path> -parse_seqids -out <db-file-path> -dbtype nucl`
    - Given the examples above, `<fna-file-path>` and `<db-file-path>` would be:
        - `C:\Users\yourname\NCBI\fasta_src\drerio_ncbi-GRCz12tu_refseq_transcripts.fna`
        - `C:\Users\yourname\NCBI\blast-2.17.0+\db\drerio_ncbi-GRCz12tu_refseq_transcripts`


### Installation: pybridizer

Requires python (e.g. the [miniforge conda distribution](https://conda-forge.org/download/)) and jupyter notebook.

First create a new virtual environment using conda:

```
conda create --name pyrbidize python ipykernel setuptools wheel
conda activate pybridize
```

Optionally, add this environment as an available kernel for jupyter notebook:

```
python -m ipykernel install --user --name=pybridize
```

Use pip to install the HCR probe design tool (from PyPI):

```
pip install pybridizer
```

Developers may wish to install in editable mode from within the cloned repo to be able to make changes to the code:

```
pip install -e .
```


### Usage

See `examples\pybridizer_basic_usage.ipynb` for a step-by-step guide for designing probes.

You can also visualize the alignment of the generated probes to the target sequence and store the alignment data in a FASTA file.

![Oligo pair alignment plot example](https://github.com/WhoIsJack/pybridizer/raw/main/_media/oligp_pair_alignment_plot_example.png)