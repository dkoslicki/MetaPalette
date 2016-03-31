#This docker file will create a docker container for the CommonKmers method

FROM cami/profiling

MAINTAINER David Koslicki version: 0.1

#Install required packages
RUN apt-get -y update && \
 apt-get install -y wget && \
 apt-get install -y g++ && \
 apt-get install -y build-essential && \
 apt-get install -y git && \
 apt-get install -y gzip

#Install python, biopython
RUN apt-get -y update && \ 
 apt-get install -y python && \
 apt-get install -y python-numpy python-scipy python-dev python-pip && \
 easy_install -f http://biopython.org/DIST/ biopython && \
 apt-get install -y libatlas-base-dev gfortran && \ 
 pip install --upgrade scipy

#Install Jellyfish
RUN apt-get -y update && \ 
 mkdir jellyfish && \ 
 cd jellyfish && \ 
 wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz && tar -xzf jellyfish-2.2.3.tar.gz && \ 
 cd /jellyfish/jellyfish-2.2.3/ && \ 
 ./configure && \ 
 make && \
 cp /jellyfish/jellyfish-2.2.3/bin/jellyfish /usr/local/bin

#Install ete3
RUN apt-get install -y python-numpy python-qt4 python-lxml python-six python-matplotlib && \
 #easy_install -U ete3 && \
 easy_install -U ete2

#Install X server for ETE
RUN apt-get install -y xvfb x11-xkb-utils xfonts-100dpi xfonts-75dpi xfonts-scalable xfonts-cyrillic x11-apps
#preface all python commands with xvfb-run 

#Install h5py
RUN apt-get install -y python-h5py

#Get the MetaPalette code, Install query per sequence and count_in_file
RUN git clone https://github.com/dkoslicki/MetaPalette.git && \
 cd /MetaPalette/src/QueryPerSeq/ && \
 g++ -I /jellyfish/jellyfish-2.2.3/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.3/.libs -Wl,--rpath=/jellyfish/jellyfish-2.2.3/.libs query_per_sequence.cc sequence_mers.hpp -l jellyfish-2.0 -l pthread -o query_per_sequence && \
 cp /MetaPalette/src/QueryPerSeq/query_per_sequence /usr/local/bin && \
 ls /MetaPalette/src/Python/ | xargs -I{} cp /MetaPalette/src/Python/{} /usr/local/sbin/ && \
 cd /MetaPalette/src/CountInFile/ && \
 g++ -I /jellyfish/jellyfish-2.2.3/include -std=c++0x -Wall -O3 -L /jellyfish/jellyfish-2.2.3/.libs -Wl,--rpath=/jellyfish/jellyfish-2.2.3/.libs count_in_file.cc -l jellyfish-2.0 -l pthread -o count_in_file && \
 cp /MetaPalette/src/CountInFile/count_in_file /usr/local/bin

#Install Bcalm
RUN wget https://github.com/Malfoy/bcalm/archive/1.tar.gz && \
 tar -zxf 1.tar.gz && \
 cd BCALM-1 && \
 sed -i 's/-march=native//g' makefile && \
 make && \
 cp bcalm /usr/local/bin

#Make training env variable
ENV TRAINING_FILES_LIST $DCKR_MNT/input/FileNames.txt
ENV INPUT_LIST_OF_FILE_NAMES $DCKR_MNT/input/InputFileNames.txt
ENV TRAINED_DATA $DCKR_MNT/MetaPaletteData

#Create tasks
RUN cp /MetaPalette/Docker/default /dckr/etc/tasks.d/default && \
 cp /MetaPalette/Docker/sensitive /dckr/etc/tasks.d/sensitive && \
 cp /MetaPalette/Docker/specific /dckr/etc/tasks.d/specific && \
 cp /MetaPalette/Docker/train /dckr/etc/tasks.d/train && \
 echo -e "\n \n NOTE: Be sure to download the training data at http://files.cgrb.oregonstate.edu/Koslicki_Lab/MetaPalette/ and uncompress it. The resulting directory should be passed to docker at the location CONT_DATABASES_DIR via: -v /path/to/uncompressed/Data:/dckr/mnt/camiref/CommonKmersData:ro"
