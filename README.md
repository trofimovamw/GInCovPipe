# nCovPopDyn pipeline

This pipeline was created as a simple tool to study the change of nucleotide diversity over time in a collection of sequences.

## Input
As an input the pipeline requires a file containing sequences and a one with a reference consenus sequence.

For the sequences it is important that they contain a sequencing-, or better, sample-date. The date must have the format **%YYYY-%mm-%dd**
and has to be either part of the sequence-name or provided in an additional tsv-file.
- If the date is part of the sequence-name, then the name should look like this: **'some_name | %YYYY-%mm-%dd'**.   
- If the date is provided in an additional file, then the file must have a column **'strain'**, containing the sequence name, and a column **'date'**, containing the date.

## Output
The pipeline creates a folder **'results'**, containing all (intermediate) outputs, with the following structure:
```
    ├── results                                 # Main results folder
    │   ├── analysis                            # Nucleotide diversity plots and table
    │   ├── bam                                 # sorted and indexed bam files
    │   ├── bins                                # binning results
    │       ├── cal_week                        # binned by calendar week
    │       ├── eq_days_10                      # binned by equal days                       
    │       ├── eq_size_100                     # binned by equal number of sequences
    │       ├── fuzzy_days_100                  # binned by equal number of sequences (fuzzy)
    │               ├── bin_*.bam               # binned sequences as BAM
    │               ├── counts_*.tsv            # count matrix                       
    │               ├── header_*.tsv            # header files (seq. name & date)
    |               ├── range_*.tsv             # range of dates of the corresponding bin
    │   ├── plots                               # Results plots (and tables)
    │   └── raw                                 # Preprocessed files
    │   
    └── ...
```
The final diversity analysis (plots and tables) can be found in the subfolder **results/plots** and the smoothed spline trajectory can be found in subfolder **results/splines**.


## How to run this pipeline - A small instruction

This is a small guide on how to run the pipeline. If you follow this instruction, it will be easier to understand what went wrong in case of any trouble.

### 1. Prerequisites
To run this pipeline, some tools have to be installed. While some are necessary (Snakemake), others are optional (Conda/Miniconda).
However, we recommend to follow all steps, since we cannot guarantee functionality otherwise.

#### 1.1 Install Conda/Miniconda - if you haven't yet

Conda will manage the dependencies of our pipeline. Instructions can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install


#### 1.2 Create the working environment

Create a new environment where the pipeline will be executed, for example like this:

```
conda create --name ncov_pipeline
```

Then to activate this environment, type:

```
conda activate ncov_pipeline
```

#### 1.3 Install Snakemake

Snakemake is the workflow management system we use. Install it like this:

```
conda install snakemake
```

### 2. Initialize the pipeline

As input the pipeline requires names of the sequence file and the reference genome, and binning parameters.
These variables are stored in [`config.yaml`](./config.yaml) and used as wildcards to create and link files with each other or as parameters for the binning.

#### 2.1 Raw sequences
The pipeline requires a file containing sequences, with the date in the sequence-name in GISAID format (date in format "%Y-%m-%d" at the end of header after a vertical bar).

For sequence files containing the date within the sequence-name, copy the file path and paste it into the variable **samples** of [`config.yaml`](./config.yaml):

  ```
  samples: "path/to/sequences/data"
  ```
If the headers in the sequence file do not contain the date, you can add it to headers using a script provided (TODO!).

#### 2.2 Reported cases data file

To compare estimated population dynamics with reported active cases, include a table in the folder [`reported_cases`](./reported_cases). Also provide the following parameters in the corresponding config field, for example like this:

  ```
  reported_cases: ["reported_cases.csv","\t","date","active_cases","%m/%d/%y"]
  ```
  
where the first element of the list is the file name with format extension, the second element is the delimiter type in this file, followed by date column name, active cases column name and a format the date is stored in.

#### 2.3 Reference consensus sequence
Copy and paste the file path of reference/consensus sequence into the variable **consensus** of [`config.yaml`](./config.yaml).

  ```
  consensus: "path/to/consensus/sequence"
  ```

#### 2.4 Binning parameters
You also have to set the parameters for some of the binning methods in [`config.yaml`](./config.yaml).
You can set the number of sequences per bin, and the number of days.
Parameters should be given as a list. Additionally, minimal bin size and maximal days span should be
provided.

```
number_per_bin: [20, 30]
days_per_bin: [7, 10, 30]
min_bin_size: 15
max_days_span: 21
```



### 3. Run

To run the pipeline, go to Snakemake directory (where the Snakefile is) and activate the conda environment you created in step 1.2. Then enter the following command to execute the pipeline:


```
snakemake --use-conda --snakefile Snakefile --cores 2
```

The ---use-conda parameter allows Snakemake to install packages that are listed in the environment file [`env.yml`](./env/env.yml). The --cores parameter defines how many CPUs the pipeline will use.

Output of each pipeline step can be found in folder **results**.
