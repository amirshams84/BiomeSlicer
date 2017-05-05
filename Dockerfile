FROM centos:latest
MAINTAINER Amir Shams <amir.shams84@gmail.com>
##############################################################
# Dockerfile Version:   1.0
# Software:             centos7
# Software Version:     7.0
# Software Website:     https://www.centos.org/
# Description:          Centos7
##############################################################

RUN yum -y update ;
RUN yum clean all ;
RUN yum install -y epel-release ;
RUN yum install -y ansible git gcc gcc-c++ make net-tools sudo which wget file patch libtool texinfo tar unzip bzip2 bzip2-devel ; 
RUN yum install -y openssl openssl-devel readline readline-devel sqlite-devel tk-devel zlib zlib-devel ncurses-devel python-pip mc ;
RUN yum install -y tbb psmisc python-devel ;
RUN yum clean all ;
##############################################################
# Software:             PIP INSTALL PACKAGES
# Software Version:     1.0
# Software Website:     -
# Description:          required javascript library
##############################################################

RUN pip install --upgrade pip
RUN pip install numpy
RUN pip install scipy
RUN pip install plotly
RUN pip install pandas
RUN pip install xlrd
RUN pip install openpyxl
RUN pip install xlwt
RUN pip install XlsxWriter
RUN pip install lxml
RUN pip install zip
RUN pip install biom-format
##############################################################
# Software:             Regular Directories
# Software Version:     1.0
# Software Website:     -
# Description:          required javascript library
##############################################################

RUN mkdir /INPUTDIR /EXECDIR /OUTPUTDIR /TESTDIR /INDEXDIR
RUN chmod -R 0755 /INPUTDIR /EXECDIR /OUTPUTDIR /TESTDIR /INDEXDIR
##############################################################
# Dockerfile Version:   1.0
# Software:             mothur
# Software Version:     1.39
# Software Website:     www.mothur.org
# Description:          mothur 
##############################################################
 
RUN wget https://github.com/mothur/mothur/releases/download/v1.39.5/Mothur.linux_64.zip -P /EXECDIR
RUN unzip /EXECDIR/Mothur.linux_64.zip -d /EXECDIR
RUN rm -rf /EXECDIR/Mothur.linux_64.zip /EXECDIR/__MACOSX
RUN chmod -R 0755 /EXECDIR/mothur
##############################################################
# Dockerfile Version:   1.0
# Software:             biomeslicer
# Software Version:     1.0
# Software Website:     .
# Description:          biomeeslicer python script
##############################################################
RUN wget https://github.com/amirshams84/biomeslicer/blob/master/biomeslicer.pyc?raw=true -O /biomeslicer.pyc
RUN chmod -R 0755 /biomeslicer.pyc

