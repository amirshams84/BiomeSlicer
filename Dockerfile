FROM amirshams/centos7:4.0

MAINTAINER Amir Shams <amir.shams84@gmail.com>
ENV ROOT=/
ENV CURRENT_PATH=.

##############################################################
# Software:             PIP INSTALL PACKAGES
# Software Version:     1.0
# Software Website:     -
# Description:          required javascript library
##############################################################
pip install numpy
pip install scipy
pip install plotly
pip install pandas
pip install biom-format
##############################################################
# Software:             javascript
# Software Version:     1.0
# Software Website:     -
# Description:          required javascript library
##############################################################
RUN git clone https://github.com/amirshams84/javascript $CURRENT_PATH/javascript
RUN chmod -R 0755 $CURRENT_PATH/javascript/ ;

##############################################################
# Software:             exec
# Software Version:     1.0
# Software Website:     -
# Description:          required execution files
##############################################################
RUN git clone https://github.com/amirshams84/exec $CURRENT_PATH/BIOME_SLICER_EXECDIR
RUN chmod -R 0755 $CURRENT_PATH/BIOME_SLICER_EXECDIR/ ;

##############################################################
# Software:             TEST_DATA
# Software Version:     1.0
# Software Website:     -
# Description:          required test files
##############################################################
RUN git clone https://github.com/amirshams84/test_data $CURRENT_PATH/BIOME_SLICER_TESTDIR
RUN chmod -R 0755 $CURRENT_PATH/BIOME_SLICER_TESTDIR/ ;

VOLUME $CURRENT_PATH/BIOME_SLICER_OUTPUTDIR

RUN wget https://raw.githubusercontent.com/amirshams84/16S_Data_Parser/master/biome_slicer.py -P $CURRENT_PATH/

CMD ["bin/bash"]
