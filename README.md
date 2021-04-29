# *G*enome-based *In*cidence Estimation of *Cov*id-19 *Pipe*line

This pipeline was created as an easy-to-use tool to study the change of nucleotide diversity over time in a collection of sequences.

## Input
As an input the pipeline requires a file containing sequences and a file with a reference consenus sequence.

For the sequences it is important that they contain a sequencing-, or better, sample-date. The date must have the format **%YYYY-%mm-%dd**
and has to be either part of the sequence-name or provided in an additional tsv-file.
- If the date is part of the sequence-name, then the name should look like this: **'some_name | %YYYY-%mm-%dd'**.   
- If the date is provided in an additional file, add the date to corresponding FASTA headers.

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
conda create --name GInCovPipe
```

Then to activate this environment, type:

```
conda activate GInCovPipe
```

#### 1.3 Install Snakemake

Snakemake is the workflow management system we use. Install it like this:

```
conda install snakemake
```

### 2. Initialize the pipeline

As input the pipeline requires names of the sequence file, the reference genome, and binning parameters.
These variables are stored in [`config.yaml`](./config.yaml) and used as wildcards to create and link files with each other or as parameters for the binning. For more information about the YAML markup format refer to documentation: https://yaml.org

#### 2.1 Raw sequences
The pipeline requires a file containing sequences, with the date in the sequence-name in GISAID format (date in format "%Y-%m-%d" at the end of header after a vertical bar).

For sequence files containing the date within the sequence-name, copy the file path and paste it into the variable **samples** of [`config.yaml`](./config.yaml):

  ```
  samples: "path/to/sequences/data"
  ```
If the headers in the sequence file do not contain the date, you can add it to headers using a script provided (TODO!).

#### 2.2 Reported cases data file

To compare estimated population dynamics with reported active cases, include a table in the folder [`reported_cases`](./reported_cases). Also provide the following parameters in the corresponding config field like this:

  ```
  reported_cases: ["reported_cases.csv","\t","date","active_cases","%m/%d/%y"]
  ```

where the first element of the list is the file name with format extension, the second element is the delimiter type in this file, date column name, active cases column name, and a format the date is stored in.

If no reported cases data is provided, leave the fields empty like this:

  ```
  reported_cases: []
  ```

#### 2.3 Reference consensus sequence
Copy and paste the file path of reference/consensus sequence into the variable **consensus** of [`config.yaml`](./config.yaml).

  ```
  consensus: "path/to/consensus/sequence"
  ```

#### 2.4 Binning parameters
You also have to set the parameters for some of the binning methods in [`config.yaml`](./config.yaml).
You can set the number of sequences per bin, and the number of days.
Parameters can be given as an array. Additionally, minimal bin size and maximal days span should be
provided.

```
number_per_bin: [20, 30]
days_per_bin: [7, 10, 30]
min_bin_size: 15
max_days_span: 21
```

If parameter **number_per_bin** is an empty list, a default mode with predefined fractions of reads (2%, 5%, 7%) is used. ALternatively, all arrays can be given in the configuration file as a list, like this:

```
number_per_bin: 
    - 20
    - 30
```

#### 2.5 Reproduction number prediction

The workflow can calculate and plot the prediction of effective reproduction number. If this prediction is desired, specify it via a Boolean variable in the configuration file, for example like this:

```
R0: y
```

If no prediction is wanted, specify it in the configuration file like this:

```
R0: n
```

Other options for specifying this parameter also work. For examples see https://yaml.org/type/bool.html

### 3. Run

To run the pipeline, go to the pipeline directory (where the Snakefile is) and activate the conda environment you created in step 1.2. Then enter the following command to execute the pipeline:

```
snakemake --use-conda --snakefile Snakefile --configfile path/to/config.yaml -j -d path/to/workdir
```

The ---use-conda parameter allows Snakemake to install packages that are listed in the environment file [`env.yml`](./env/env.yml). With parameter --configfile you can give the configuration file [`config.yml`], described above. The -j parameter determines the number of available CPU cores to use in the pipeline. Optionally you can provide the number of cores, e.g. -j 4. With parameter -d you can set the work directory, i.e. where the results of the pipeline are written to.

## Output
The pipeline creates a folder **'results'**, containing all (intermediate) outputs, with the following structure:
```
    ├── results                                 # Main results folder
    │   ├── bam                                 # sorted and indexed bam files
    │   ├── bins                                # binning results
    │       ├── cal_week                        # binned by calendar week
    │       ├── eq_days_10                      # binned by equal days                       
    │       ├── eq_size_100                     # binned by equal number of sequences
    │       └── fuzzy_days_100                  # binned by equal number of sequences (fuzzy)
    │               ├── bin_*.bam               # binned sequences as BAM
    │               ├── bin_*.bai               # index files                       
    │               ├── header_*.tsv            # header files (seq. name & date)
    |               ├── range_*.tsv             # range of dates of the corresponding bin
    |               └── list_of_files.tsv       # list of file names in the binning mode
    │   ├── bins_results                        # Individual binning results plots (and tables)
    │   ├── meta                                # Meta information about all used sequences (name and collection date)
    │   ├── interpolation                       # Plots and tables for final interpolated trajectory
    |       ├── interpolation.csv               # table with interpolated population size estimates
    │       ├── estimates.csv                   # table wih raw population size estimates                       
    │       ├── interpolation.pdf               # plot of interpolated population size
    |       └── wdots_interpolation.pdf         # plot of interpolated population size with dot size scaled by bin size
    │   ├── r0                                  # Reproduction number estimate
    |       ├── r0.csv                          # table with daily reproduction number estimates
    │       └── r0.pdf                          # plot of daily reproduction number estimates                       
    │   └── raw                                 # Preprocessed files
    │   
    └── ...
```
The final diversity analysis (plots and tables) can be found in the subfolder **results/bins_results** and the smoothed trajectory can be found in subfolder **results/interpolation**. The final output (depending on the setup of the pipeline) contains the following tables and plots:

- In folder *interpolation*:
    - *estimates.csv*: table containing estimates of population size for all binning strategies
    - *interpolation.csv*: table containing interpolated final trajectory of population size; the trajectory is calculated by  combining all binning strategies
    - *interpolation.pdf*: plot of interpolated trajectory, optionally overlayed with reported cases trajectory (if the corresponding table was given)
    - *wdots_interpolation.pdf*: plot of interpolated trajectory with point estimate dots scaled by corresponding sub-sample size, optionally overlayed with reported cases trajectory (if the corresponding table was given)
- In folder *r0* (if option to calculate reproduction number is chosen):
    -  *r0.csv*: table containing daily reproductive number estimates; calculated from the interpolated trajectory
    -  *r0.pdf*: plot of daily reproductive number estimates with confidence interval 

