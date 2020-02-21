FROM ubuntu:14.04 

MAINTAINER Ricardo R. da Silva <ridasilva@usp.br>

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    lzip \
    curl \
    wget \
    libtclap-dev \
    libpotrace0  \
    libpotrace-dev  \
    libocrad-dev \
    openbabel \
    libopenbabel-dev  \
    libgraphicsmagick++-dev \
    libnetpbm10-dev \
&&  apt-get clean \
&&  rm -rf /var/lib/lists/*
ENV CXXFLAGS -pthread 

RUN mkdir -p /tmp/OSRA \
    && curl http://ftp.gnu.org/gnu/ocrad/ocrad-0.25.tar.lz \
    | tar --lzip -xC /tmp/OSRA \
    && cd /tmp/OSRA/ocrad-0.25 \
    && ./configure \
    && make \
    && make install

RUN mkdir -p /tmp/OSRA &&  \
    cd /tmp/OSRA &&  \
    wget http://cactus.nci.nih.gov/osra/gocr-0.50pre-patched.tgz && \
    tar -xvf gocr-0.50pre-patched.tgz && \
    cd gocr-0.50pre-patched && \
    ./configure &&  \
    make libs && \
    make all install

RUN mkdir -p /tmp/OSRA \
    && cd /tmp/OSRA \
    && wget http://downloads.sourceforge.net/project/osra/osra/2.0.1/osra-2.0.1.tgz \
    && tar -xvf /tmp/OSRA/osra-2.0.1.tgz \
    && cd /tmp/OSRA/osra-2.0.1 \
    && ./configure \
    && make all \
    && make install
RUN cd && rm -rf /tmp/OSRA

RUN \
    apt-get install git default-jdk -y

#RUN \
#    apt-get install unzip -y
#
#RUN \
#     cd /opt && \
#     wget https://downloads.apache.org/maven/maven-3/3.6.3/binaries/apache-maven-3.6.3-bin.zip && \
#     unzip apache-maven-3.6.3-bin.zip && \
#     ln -s /opt/apache-maven-3.6.3 /opt/maven && \
#     echo 'export JAVA_HOME=/usr/lib/jvm/default-java' > ~/.bashrc && \
#     echo 'export M2_HOME=/opt/maven' > ~/.bashrc && \
#     echo 'export MAVEN_HOME=/opt/maven' > ~/.bashrc && \
#     echo 'export PATH=${M2_HOME}/bin:${PATH}' > ~/.bashrc && \
#     . ~/.bashrc
#
#RUN \
#     echo $MAVEN_HOME && \
#     ls $JAVA_HOME && \
#     /opt/maven/bin/mvn -version
#
#RUN \
#    cd /bin && \
#    git clone https://bitbucket.org/mjw99/chemextractor.git && \
#    cd chemextractor ; /opt/maven/bin/mvn clean package ; dpkg -i ./target/*.deb
#

COPY chemextractor_1.0~SNAPSHOT_all.deb /tmp/
RUN dpkg -i /tmp/chemextractor_1.0~SNAPSHOT_all.deb

#install pdftotext
RUN apt-get install -y poppler-utils

#install gnfinder 
RUN \
	wget https://github.com/gnames/gnfinder/releases/download/v0.9.1/gnfinder-v0.9.1-linux.tar.gz -P /tmp && \
	tar xf /tmp/gnfinder-v0.9.1-linux.tar.gz -C /bin

#install npmine
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh  && \
    bash /tmp/miniconda.sh -f -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

#COPY environment.yml /tmp/environment.yml

#RUN conda env create -f /tmp/environment.yml
RUN conda create -n nplibrary python=3 -y
RUN echo "source activate nplibrary" > ~/.bashrc
ENV PATH /opt/conda/envs/nplibrary/bin:$PATH

RUN conda init bash
RUN bash ~/.bashrc
RUN conda install -n nplibrary -c rdkit rdkit -y
RUN pip install pandas requests bs4 configparser ipykernel jupyter
RUN pip install seaborn biopython 
RUN pip install snakemake 
#RUN pip install git+https://gitlab.com/rsilvabioinfo/npmine_library
RUN python -m ipykernel install --user --name nplibrary --display-name nplibrary
COPY . .
RUN python setup.py install

RUN mkdir -p /home/npmine/
WORKDIR /home/npmine/
#COPY ../notebooks notebooks
EXPOSE 8888

#cleanup
RUN rm -rf /var/lib/apt/lists/*
RUN rm -rf /tmp/*

CMD ["bash"]
