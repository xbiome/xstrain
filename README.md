# StrainPanDA -- A strain analysis pipeline based on pangenome

StrainPanDA is a tool that deconvolute pangenome coverage into strain composition and strain gene profile.

## Installation

### Local installation

1. install nextflow

```sh
curl -s https://get.nextflow.io | bash
```

2. install PanPhlAn (for mapping short reads to pan-genome databases)

https://github.com/SegataLab/panphlan/wiki#download-the-panphlan-software

3. install [R](https://www.r-project.org/) and [strainpandar](src/strainpandar): decomposing pangenome into strains

Install required packages `dplyr`, `foreach`, `MASS`, `NMF`

```sh
tar -czf strainpandar.tar.gz src/strainpandar
R CMD INSTALL strainpandar.tar.gz
```

### Using the strainpanda container (recommended)

We provide two docker files to build the two images required by the nextflow pipeline.

```
docker build -t strainpanda-mapping:dev . -f docker/Dockerfile_mapping
docker build -t strainpanda-strainpandar:dev . -f docker/Dockerfile_strainpandar
```

### Pre-built pangenome databases

TODO: url for downloading

## Run analysis

### Run full analysis
Download test data to a foler e.g. `reads`

```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/005/SRR5813295/SRR5813295_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/005/SRR5813295/SRR5813295_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/006/SRR5813296/SRR5813296_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/006/SRR5813296/SRR5813296_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/007/SRR5813297/SRR5813297_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/007/SRR5813297/SRR5813297_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/008/SRR5813298/SRR5813298_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/008/SRR5813298/SRR5813298_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/009/SRR5813299/SRR5813299_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR581/009/SRR5813299/SRR5813299_2.fastq.gz
```

Create a species list

```sh
echo "Faecalibacterium prausnitzii" > species_list.txt
```

#### With docker

```sh
PATH_TO_REPO/main.nf -profile docker --ref_path PATH_TO_REFERENCE --path reads/ --ref_list species_list.txt
```

#### Local

```sh
PATH_TO_REPO/main.nf --ref_path PATH_TO_REFERENCE --path reads/ --ref_list species_list.txt
```

### Run strainpandar

Assuming you have a count matrix generated, we provide a [Rscript wrapper](bin/run_strainpandar.r) to run strainpandar. A sample count matrix can be found [here](data/Faecalibacterium-prausnitzii-202009.counts.csv).


```sh
Rscript bin/run_strainpandar.r -c data/Faecalibacterium-prausnitzii-202009.counts.csv -r data/refs/Faecalibacterium-prausnitzii-202009 -o work -t 8 -m 8 -n 0
```

The parameters of the wrapper can be found from the help message:

```sh
> Rscript bin/run_strainpandar.r -h
A wrapper script to perform strain decomposition using strainpandar package.
Usage: bin/run_strainpandar.r [-[-help|h]] [-[-counts|c] <character>] [-[-reference|r] <character>] [-[-output|o] [<character>]] [-[-threads|t] [<integer>]] [-[-max_rank|m] [<integer>]] [-[-rank|n] [<integer>]]
    -h|--help         Show this help message
    -c|--counts       Gene-sample count matrix (CSV file) obtained from mapping reads to a reference pangenome [required]
    -r|--reference    Pangenome database path [required]
    -o|--output       Output prefix [default: ./strainpandar]
    -t|--threads      Number of threads to run in parallele [default: 1]
    -m|--max_rank     Max number of strains expected [default: 8]
    -n|--rank         Number of strains expected. If 0, try to select from 1 to `max_rank`. If not 0, overwrite `max_rank`. [default: 0]
```


## Outputs

Ouput files from the above run:

```sh
strainpanda_out/
├── Faecalibacterium-prausnitzii-202009.counts.csv  ## Merged count matrix (gene family by sample)
├── Faecalibacterium-prausnitzii-202009_mapping     ## Sample specific count files
│   ├── SRR5813295_Faecalibacterium-prausnitzii-202009.csv.bz2
│   ├── SRR5813296_Faecalibacterium-prausnitzii-202009.csv.bz2
│   ├── SRR5813297_Faecalibacterium-prausnitzii-202009.csv.bz2
│   ├── SRR5813298_Faecalibacterium-prausnitzii-202009.csv.bz2
│   └── SRR5813299_Faecalibacterium-prausnitzii-202009.csv.bz2
├── Faecalibacterium-prausnitzii-202009_strainpandar_out ## strainpandar outputs
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.anno_strain_sample.pdf ## annotation to the closest reference
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.genefamily_strain.csv ## gene family-strain matrix
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.genefamily_strain.pdf ## heatmap visualization
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.rds ## R object contains strainpandar results
│   ├── Faecalibacterium-prausnitzii-202009.strainpanda.strain_sample.csv ## strain-sample matrix
│   └── Faecalibacterium-prausnitzii-202009.strainpanda.strain_sample.pdf ## barplot visualization
└── pipeline_info  ## pipeline run statistics
    ├── strainpanda_DAG.svg
    ├── strainpanda_report.html
    ├── strainpanda_timeline.html
    └── strainpanda_trace.txt
```