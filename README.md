# MetaPalette #

## What is MetaPalette? ##
CommonKmers is a k-mer based bacterial community reconstruction technique that utilizes sparsity promoting ideas from the field of compressed sensing to reconstruct the composition of a bacterial community. This method allows for strain-level abundance estimation, and can quantify the evolutionary distance between organisms in the sample and in the training database (thereby allowing for successful classification even with incomplete training data).


## How Do I Install CommonKmers? ##
###Build from source###
You will need the Kmer counting tool Jellyfish to be installed. Please see [the Jellyfish installation page](http://www.genome.umd.edu/jellyfish.html) for installation directions. Briefly, this can be installed using:

```bash
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz
tar -xf jellyfish-2.2.3.tar.gz
cd jellyfish-2.2.3
./configure
make
```
The binary will then be located in ``jellyfish-2.2.3/bin/``.

You will need to download [this data repository](http://www.math.oregonstate.edu/~koslickd/CommonKmersData.tar.gz), and then extract using ``tar -xf CommonKmersData.tar.gz``. This folder contains all the default training data. This can be accomplished with:

```bash
curl http://www.math.oregonstate.edu/~koslickd/CommonKmersData.tar.gz > CommonKmersData.tar.gz
tar -xf CommonKmersData.tar.gz
```

Please refer to [the Julia installation page](http://julialang.org/downloads/) to install Julia.
You will need to add the HDF5 and ArgParse packages. These can be added using `Pkg.add("HDF5")` and `Pkg.add("ArgParse")`. In Ubuntu, this can be accomplished with something like:

```bash
apt-get install -y software-properties-common python-software-properties
add-apt-repository -y ppa:staticfloat/juliareleases
add-apt-repository -y ppa:staticfloat/julia-deps
apt-get -y update || echo "ok" 
apt-get install -y julia
apt-get install -y hdf5-tools
julia -e 'Pkg.add("HDF5"); Pkg.add("ArgParse");'
```

You will also need to compile the ``query_per_sequence`` code using a command such as:
```bash
g++ -I /jellyfish/jellyfish-2.2.3/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.3/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.3/.libs query_per_sequence.cc sequence_mers.hpp -o query_per_sequence
```

###Or Use Docker###
A Dockerfile is included in this repository. See the [Docker homepage](https://www.docker.com/) for more information.

You will need to download [this data repository](http://www.math.oregonstate.edu/~koslickd/CommonKmersData.tar.gz), and then extract using ``tar -xf CommonKmersData.tar.gz``. This folder contains all the default training data.

You can either pull the docker image from DockerHub using
```bash
docker pull dkoslicki/commonkmers
```

Or you can build the docker image from the Dockerfile by cloning the repository, starting Docker, and then in the ``CommonKmers/Docker`` folder, using the command:
```bash 
docker build -t username/imagename .
```


## Running the program ##
####From the command line####
To classify a sample using the Julia version, use the ``ClassifyFull.jl`` command located in ``CommonKmers/src/Julia``. An example of running the program in the sensitive mode using 48 threads and a minimum quality score (for kmers to be counted) of C (phred33 ascii code 35) is given by

```julia
julia -p 48 Classify.jl -d /path/to/CommonKmersData/ -o /path/to/output/file.profile -i /path/to/input/file.fastq -Q C -k sensitive -j /path/to/./jellyfish -q /path/to/./query_per_sequence --normalize
```

FASTQ and FASTA files are acceptable input. Note that if FASTA files are used, no error correction will be done (which can lead to poor results).

The optional flag ``--save_y`` will save the normalized common kmer counts into a file called ``/path/to/ouput/file.fastq-y30.txt``. After running the script with this flag, the script can then be re-run using the optional flag ``--re_run`` with a different ``--kind`` specified (this will significantly speed this and any other subsequent runs). Note that if you change the quality score ``-Q``, you must ``--save_y`` before you can use the ``--re_run`` flag again.

The optional flag ``--save_x`` will save the reconstruction to a text file. This is necessary for plotting features.

The optional flag ``--normalize`` will normalize the output profile to sum to 1 (so it will appear that 100% of the sample has been classified). This is similar to the default options of MetAPhlAn.


####Using Docker####
To run the tool from docker, mount the appropriate folders and run using the following command:
```bash
docker run --rm -e "QUALITY=C" -e "DCKR_THREADS=48" -v /path/to/CommonKmersData:/dckr/mnt/camiref/CommonKmersData:ro -v /path/to/Output:/dckr/mnt/output:rw -v /path/to/Input:/dckr/mnt/input:ro -t username/imagename [type]
```
In the input folder must be a collection of gzipped FASTQ (or FASTA) files, as well as a file (called ``sample.fq.gz.list`` (given by the docker image environmental variable ``$CONT_FASTQ_FILE_LISTING``) listing the files on which to run the tool.
Here ``[type]`` is one of ``default, sensitive, specific``.
The ``--rm`` flag deletes temporary files after exit (otherwise they might persist in ``/var/lib/docker/volumes`` or the like).
If the environmental variable ``QUALITY`` is not passed to docker (via ``-e QUALITY=<ascii character>``), a default value of "C" will be used. 



## Output format ##
The output format complies with the [CAMI format](https://github.com/CAMI-challenge/contest_information/blob/master/file_formats/CAMI_TP_specification.mkd).
The docker complies with the [Bioboxes profiling format 0.9](https://github.com/bioboxes/rfc/tree/master/data-format).

## Recommendations ##
I recommend using a quality score roughly equal to the average first quartile quality score in the file. This can be found with the following commands:
```bash
#Convert non ACTGN characters to N
awk '{if(NR%4==2){gsub(/[^ACGT]/,"N");print $0}else{print $0}}' input.fq > input_ACTGN.fq 
#Compute average first quartile quality score
/Fastx/bin/./fastx_quality_stats -i input_ACTGN.fq -Q33 | cut -f7 | sed -n '1!p' | awk '{a+=$1} END{print a/NR}' | awk '{printf "%.0f",$1}'
```
The FastX toolbox can be downloaded [here](http://hannonlab.cshl.edu/fastx_toolkit/).


## Custom Training Databases ##
If you wish to use a custom training database, the following steps must be performed:

0. Install Bcalm
1. Create a directory to contain the training data (called ``CommonKmerTrainingData`` below).
2. Create an acceptable taxonomy for the training genomes, and place it in the ``CommonKmerTrainingData`` folder.
3. Create a file consisting of the full paths of the training genomes, and save this to a file (for example, ``FileNames.txt``).
4. Compile the code contained in ``CommonKmers/src/CountInFile/``.
5. Run the script ``Train.jl``.

Alternatively, you can use Docker (though an acceptable taxonomy still needs to be created).

####Install Bcalm####
To install Bcalm, do something like the following:
```bash
wget https://github.com/Malfoy/bcalm/archive/1.tar.gz && \
 tar -zxf 1.tar.gz && \
 cd bcalm-1 && \
 make && \
 cp bcalm /usr/local/bin
```
Note that the [current Bcalm git repository](https://github.com/Malfoy/bcalm) does not compile correctly, so you must install from the release.

####Creating custom taxonomy####
For each genome in ``FileNames.txt`` (and in the same order), a taxonomy file must be created. This file MUST be a newline delimitated file with each line having the following format:
```bash
<organismName>\t<TaxID>\t<TaxPath>
```

``<organismName>`` must be a unique identifier for each genome.

``<TaxID>`` must be a unique TaxID for each genome

``<TaxPath>`` must be a pipe delimitated list that gives the taxonomy of the given organism. The format is: 

```
k__<KingdomTaxID>_<KingdomName>|p__<PhylumTaxID>_<PhylumName>|c__<ClassTaxID>_<ClassName>|o__<OrderTaxID>_<OrderName>|f__<FamilyTaxID>_<FamilyName>|g__<GenusTaxID>_<GenusName>|s__<SpeciesTaxID>_<SpeciesName>|t__<StrainTaxID>_<StrainName>
```

The taxonomy is only required at the kingdom level, with lower levels being optional. Missing ranks can be included with ``||``.

An example line is as follows:

```
1184607_Austwickia_chelonae_NBRC_105200	1184607	k__2_Bacteria|p__201174_Actinobacteria|c__1760_Actinobacteria|o__2037_Actinomycetales|f__85018_Dermatophilaceae|g__1184606_Austwickia|s__100225_Austwickia_chelonae|t__1184607_Austwickia_chelonae_NBRC_105200
```

For your convenience, the script ``CommonKmers/src/Taxonomy/generate_taxonomy_taxid.py`` generates such a taxonomy using the NCBI taxonomy. This file must be placed in the ``CommonKmerTrainingData`` folder and MUST be called ``Taxonomy.txt``.

####Compile the ``count_in_file`` code####
The ``/CommonKmers/src/CountInFile/count_in_file.cc`` code can be compiled using a command like:

```bash
g++ -I /jellyfish/jellyfish-2.2.3/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.3/.libs -l jellyfish-2.0 -l pthread -Wl,--rpath=/jellyfish/jellyfish-2.2.3/.libs count_in_file.cc -o count_in_file
```

####Run the script ``Train.jl``####
The script ``Train.jl`` can be called using a command such as:
```bash
julia -p 48 Train.jl -i FullFileNames.txt -o /path/to/output/CommonKmerTrainingData/ -b /path/to/./bcalm -r /path/to/fast/IO/device/ -j /path/to/jellyfish -c /path/to./count_in_file -s 500 -t 20
```

The option ``-s`` specifies how many training genomes at a time are held in memory. Increasing/decreasing this increases/decreases the amount of RAM used.

The option ``-t`` specifies how many jellyfish instances are created. Too many will cause disk thrashing (default is 20).

The option ``-r`` specifics the location of a temporary folder on a fast IO device. Unfortunately Bcalm uses a considerable amount of file IO, and so a fast storage device is required. Note that you can create a RAM disk using a command like:

```bash
mkdir /tmp/ramdisk; chmod 777 /tmp/ramdisk
sudo mount -t tmpfs -o size=100G tmpfs /tmp/ramdisk/
```

To unmount the RAM disk, use a command like:
```bash
sudo umount -v /tmp/ramdisk
```

Note that the time required to complete the training step can be considerable (depending on hardware available). Using 48 cores and 256GB of RAM, training on ~7,000 genomes can take upwards of a week.

####Using Docker####

Alternatively, after creating the acceptable taxonomy, Docker can be used to form the training data. You will need to have access to a folder containing all the uncompressed training fasta/fastq files. You will also need to create a file (name it ``sample.fna.list``) that contains all the base names of the training fasta/fastq files, and put this in the same folder.

Docker can then be called with
```bash
docker run --rm --privileged -e "DCKR_THREADS=48"  -e "RAM_DISK_SIZE=100G" -v /path/to/input/data:/dckr/mnt/input:ro -v /path/to/output/folder:/dckr/mnt/output:rw -t username/imagename train
```
The flag ``--privileged`` is required since docker will then be allowed to automatically create the RAM disk. Note the default RAM disk size is 10G (I suggest using around half the available RAM).

####Run the ``Classify.jl`` script####
You can now run the ``Classify.jl`` script as before, but this time utilizing the directory ``CommonKmerTrainingData`` for the option ``-d``.
## Contact ##
For issues with this software, contact david.koslicki@math.oregonstate.edu

## License ##
This project is released under the GPL-3 License. Please view the [LICENSE](LICENSE)
file for more details.


## Contributors ##
+ David Koslicki (all main source code, unless otherwise noted)
+ Daniel Alonso Alemany (early version of the software)
+ Daniel Falush
+ Nam Nguyen
