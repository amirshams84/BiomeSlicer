FROM amirshams/centos7:1.0

MAINTAINER Amir Shams <amir.shams84@gmail.com>
ENV ROOT=/
ENV CURRENT_PATH=.
CMD ["/bin/bash"]
##############################################################
# Software:             PIP INSTALL PACKAGES
# Software Version:     1.0
# Software Website:     -
# Description:          required javascript library
##############################################################
RUN pip install numpy
RUN pip install scipy
RUN pip install plotly
RUN pip install pandas
RUN pip install biom-format
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
RUN git clone https://github.com/amirshams84/exec $CURRENT_PATH/BIOM_SLICER_EXECDIR
RUN chmod -R 0755 $CURRENT_PATH/BIOM_SLICER_EXECDIR/ ;

##############################################################
# Software:             TEST_DATA
# Software Version:     1.0
# Software Website:     -
# Description:          required test files
##############################################################
RUN git clone https://github.com/amirshams84/test_data $CURRENT_PATH/BIOM_SLICER_TESTDIR
RUN chmod -R 0755 $CURRENT_PATH/BIOM_SLICER_TESTDIR/ ;

VOLUME $CURRENT_PATH/BIOM_SLICER_OUTPUTDIR

RUN wget https://raw.githubusercontent.com/amirshams84/16S_Data_Parser/master/biom_slicer.py -P $CURRENT_PATH/

CMD ["bin/bash"]
