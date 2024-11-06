# ATTNETION!!!
- The program has been updated and moved to [PhaBOX 2](https://github.com/KennthShang/PhaBOX), which is more user-friendly. In the new version, PhaMer is generalized to all kinds of viruses, more than just phages. Hope you will enjoy it. This folder will be no longer maintained. 


- Our web server for viruses-related tasks (including virus identification, taxonomy classification, lifestyle prediction, host prediction, and protein annotation) is available! You can visit [Web Server](http://phage.ee.cityu.edu.hk/) to use the GUI. We also provided more detailed intermediate files and visualization for further analysis.

![PhaMer](logo.jpg)

PhaMer is a python library for identifying bacteriophages from metagenomic data. PhaMer is based on a Transorfer model and rely on protein-based vocabulary to convert DNA sequences into sentences. 


1. The paper is accepted by Briefings in Bioinformatics. If you use PhaMer for your research, please use the following citation: 
```
Jiayu Shang, Xubo Tang, Ruocheng Guo, Yanni Sun, Accurate identification of bacteriophages from metagenomic data using Transformer, Briefings in Bioinformatics, 2022;, bbac258, https://doi.org/10.1093/bib/bbac258
```
2. Thanks for the help from [@sjaenick](https://github.com/sjaenick), the program is much smooth in this currect version.


# Overview
The main function of PhaMer is to identify phage-like contigs from metagenomic data. The input of the program should be fasta files and the output will be a csv file showing the predictions. Since it is a Deep learning model, if you have GPU units on your PC, we recommand you to use them to save your time. 

If you have any trouble installing or using PhaMer, please let us know by emailing us (jyshang2-c@my.cityu.edu.hk).


## Required Dependencies
* Python 3.x
* Numpy
* Pandas
* Pytorch>1.8.0
* [Diamond](https://github.com/bbuchfink/diamond)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [git-lfs](http://arfc.github.io/manual/guides/git-lfs)


If you want to use the gpu to accelerate the program:
* cuda
* Pytorch-gpu

* For cpu version pytorch: `conda install pytorch torchvision torchaudio cpuonly -c pytorch`
* For gpu version pytorch: Search [pytorch](https://pytorch.org/) to find the correct cuda version according to your computer

### Quick install
*Note*: we suggest you to install all the package using conda (both miniconda and [Anaconda](https://anaconda.org/) are ok).

After cloning this respository, you can use anaconda to install the **PhaMer.yaml**. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: `conda env create -f PhaMer.yaml -n phamer`


### Prepare the database and environment
Due to the limited size of the GitHub, we zip the database. Before using PhaMer, you need to unpack them using the following commands.

1. When you use PhaMer at the first time
```
cd PhaMer/
conda env create -f PhaMer.yaml -n phamer
conda activate phamer
cd database/
bzip2 -d database.fa.bz2

# initial files
pip install gdown
gdown  --id 11QsvDJQwRBO-LWAXrA07N5WpxHcQJv8Y
cd ..
```
*Note:* 
* Because the parameter is larger than 100M, please make sure you have downloaded transformer.pth correctly.
* if you cannot download the *transformer.pth* from the command lines above, please use the [Google Drive link](https://drive.google.com/file/d/1d_6DGhN4Q-NZPKHEIo4yD4InLkP2U1rI/view?usp=sharing) to download it and place it in the *database/* folder. 
* You can use command `md5sum -c check.md5` in *database/* folder to check whether the parameter file is correct (thanks for the suggestion from [@nick-youngblut](https://github.com/nick-youngblut).



2. If the example can be run without any bugs, you only need to activate your 'phamer' environment before using PhaMer.
```
conda activate phamer
```


## Usage

```
python PhaMer_preprocessing.py [--contigs INPUT_FA] [--len MINIMUM_LEN] [--midfolder DIR] [--proteins PROTEIN_FA] [--threads NUM] [--dbdir DR]
python PhaMer.py [--out OUTPUT_CSV] [--reject THRESHOLD] [--midfolder DIR] [--threads NUM] [--dbdir DR]
```

**Options**


      --contigs INPUT_FA
                            input fasta file
      --len MINIMUM_LEN
                            predict only for sequence >= len bp (default 3000)
      --proteins PROTEIN_FA
                            An optional protein file. If you have already annotated your contigs, you can use them as the inputs. 
                            Otherwise, PhaMer will run prodigal to translate your contigs.
      --threads NUM
                            Number of threads to run PhaMer (default 8)
      --dbdir DR
                            An optional path to store the database directory (default database/)
      --out OUTPUT_CSV
                            The output csv file (prediction)
      --reject THRESHOLD
                            Threshold to reject prophage. The higher the value, the more prophage will be rejected (default 0.3)
      --midfolder DIR
                            Folder to store the intermediate files (default phamer/)

**Example**

Prediction on the example file:

    python PhaMer_preprocessing.py --contigs test_contigs.fa
    python PhaMer.py --out example_prediction.csv

The prediction will be written in *example_prediction.csv*. The CSV file has three columns: contigs names, prediction, and prediction score. The test_contig.fasta contain a phage genome, so the output is phage.
    
### References
The paper is accepted by Briefings in Bioinformatics and you can find it via: [BIB version](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbac258/6620872)

The arXiv version can also be found via: [arXiv version](http://arxiv.org/abs/2201.04778)

### Contact
If you have any questions, please email us: jyshang2-c@my.cityu.edu.hk


