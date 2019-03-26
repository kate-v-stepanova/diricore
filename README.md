#  Diricore pipeline installation
This is a diricore version from E.Stepanova. For the other versions there is other documentation (in docs directory).

```
git clone https://github.com/kate-v-stepanova/diricore.git
cd diricore
conda create -n diricore python=2.7
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
pip install -r requirements.txt
```

In `diricore/programs` there are installed packages required by diricore.  In ideal case everything will work smoothly just like that. If not, install the packages manually:

1. Samtools: version 0.1.19

`https://sourceforge.net/projects/samtools/files/samtools/0.1.19/`

2. cutadapt 

```
pip install cutadapt
```

3. Bowtie: version 2.0.6

```
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.6/
unzip bowtie2-2.0.6*.zip
cp -R bowtie2-2.0.6 diricore_pipeline/programs
```

4. One of the scripts requires [seqtk](https://github.com/lh3/seqtk)

```
git clone https://github.com/lh3/seqtk
cd seqtk
make
```
