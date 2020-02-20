FROM ubuntu:18.04
#FROM openjdk:8
#MAINTAINER Gert wohlgemuth <wohlgemuth@ucdavis.edu>

#installing build tools
RUN \
	apt-get update && \
	apt-get install -y build-essential automake checkinstall git cmake subversion openjdk-8-jdk wget

#fetching source codes
RUN \
	cd /tmp && \
    git clone https://github.com/metamolecular/gocr-patched.git && \
#	svn checkout svn://svn.code.sf.net/p/openbabel/code/openbabel/trunk /tmp/openbabel && \
	svn checkout svn://svn.code.sf.net/p/osra/code/tags/2.0.1 /tmp/osra

#installing dependencies for osra
#libgraphicsmagick3 
#E: Package 'libgraphicsmagick3' has no installation candidate
RUN \
	apt-get install -y libtclap-dev libpotrace0  libpotrace-dev  libocrad-dev libgraphicsmagick++1-dev libgraphicsmagick++1-dev libgraphicsmagick++3 && \
	apt-get install -y libeigen3-dev libgraphicsmagick1-dev libnetpbm10-dev libpoppler-dev libpoppler-cpp-dev libleptonica-dev wget tesseract-ocr tesseract-ocr-eng

#apt-get install -y openbabel libopenbabel-dev
RUN \
    echo 'deb http://us.archive.ubuntu.com/ubuntu/ trusty universe' >> /etc/apt/sources.list && \
    echo 'deb http://security.ubuntu.com/ubuntu xenial-security main' >> /etc/apt/sources.list && \
    apt-get update

RUN apt-get install -y openbabel libopenbabel-dev libgraphicsmagick3


RUN \
    wget -O /tmp/openbabel.tgz https://sourceforge.net/projects/osra/files/openbabel-patched/openbabel-2.3.2-tr1-memory.tgz && \
    cd /tmp/ && \
    tar -xvf openbabel.tgz && \
    cd openbabel-2.3.2-tr1-memory && \
    mkdir build  && \
    cd build  && \
    cmake ..  && \
    make -j2  && \
    make install

#patching gocr
RUN \
	cd /tmp/gocr-patched && \
	./configure && \
	make libs && \
	make all install

#installing osra
RUN \
	cd /tmp/osra && \
        wget https://gist.githubusercontent.com/mcs07/7b722cfafe8bad81aa69/raw/8886e5feb2f97183b3117b6a84e46c19c05a807b/osra-adaptiveThreshold.diff && \
        git apply osra-adaptiveThreshold.diff && \
	./configure --with-tesseract && \
	make all && \
	make install && \
	echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib >> ~/.bashrc

#useful converter
RUN \
	apt-get install -y imagemagick

#install gnfinder 
RUN \
    cd /bin && \
    git clone https://github.com/gnames/gnfinder.git && \
    cd gnfinder && \
    make install

#install chemextractor/oscarpdf2json
RUN \
    cd /bin && \
    git clone https://bitbucket.org/mjw99/chemextractor.git && \
    cd chemextractor ; mvn clean package ; sudo dpkg -i ./target/*.deb

#install npmine
RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O /tmp/miniconda.sh  && \
    echo 'e1045ee415162f944b6aebfe560b8fee */tmp/miniconda.sh' | md5sum -c - && \
    bash /tmp/miniconda.sh -f -b -p /opt/conda && \
    /opt/conda/bin/conda install --yes -c conda-forge \
    /opt/conda/bin/pip install --upgrade pip && \
    rm /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml
RUN echo "source activate nplibrary" > ~/.bashrc
ENV PATH /opt/conda/envs/nplibrary/bin:$PATH

#cleanup
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/*

CMD ["bash"]
