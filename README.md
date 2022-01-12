![PhaMer](logo.jpg)

PhaMer is a python library for identifying bacteriophages from metagenomic data. PhaMer is based on a Transorfer model and rely on protein-based vocabulary to convert DNA sequences into sentences. 

# Overview
The main function of PhaMer is to identify phage-like contigs from metagenomic data. The input of the program should be fasta files and the output will be a csv file showing the predictions. Since it is a Deep learning model, you can either re-train your own model or use our pre-train model for prediction. However, if you do not have GPU units on your PC, we recommand you to use the pre-trained model to save your time. We will continuously update the database for your convenience.

If you have any trouble installing or using PhaMer, please let us know by opening an issue on GitHub or emailing us (jyshang2-c@my.cityu.edu.hk).


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

### An easiler way to install
*Note*: we suggest you to install all the package using conda (both miniconda and [Anaconda](https://anaconda.org/) are ok).

After cloning this respository, you can use anaconda to install the **PhaMer.yaml**. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: `conda env create -f PhaMer.yaml`

* For cpu version pytorch: `conda install pytorch torchvision torchaudio cpuonly -c pytorch`
* For gpu version pytorch: Search [pytorch](https://pytorch.org/) to find the correct cuda version according to your computer

### Prepare the database
Due to the limited size of the GitHub, we zip the database. Before using PhaMer, you need to unpack them using the following commands.

```
cd PhaMer/database
bzip2 -d database.fa.bz2
cd ..
```

## Usage
