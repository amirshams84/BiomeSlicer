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
RUN git clone https://github.com/amirshams84/exec $CURRENT_PATH/exec
RUN chmod -R 0755 $CURRENT_PATH/exec/ ;

##############################################################
# Software:             TEST_DATA
# Software Version:     1.0
# Software Website:     -
# Description:          required test files
##############################################################
RUN git clone https://github.com/amirshams84/test_data $CURRENT_PATH/test_data
RUN chmod -R 0755 $CURRENT_PATH/test_data/ ;

VOLUME $CURRENT_PATH/microbiome_slicer_results_folder

RUN wget https://raw.githubusercontent.com/amirshams84/16S_Data_Parser/master/microbiome_slicer.py -P $CURRENT_PATH/

CMD ["bin/bash"]
