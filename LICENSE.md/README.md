# LFIT_GUI
A GUI-enabled lens modeling tool

Dependencies
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
lfit_gui is written in python, and uses some standard python packages such as numpy, scipy, matplotlib, etc. Some non-standard python packages are also required. The easiest way of installing those packages (if you haven't) will be through pip install. 

1. lmfit (a non-linear fitter)
https://pypi.python.org/pypi/lmfit/0.8.0

Note that lfit_gui uses an older version of lmfit (0.8.0), which is unfortunately not compatible with its later releases. So if you were to install lmfit via pip install, you would need to specify the version number by “pip install lmfit==0.8.0”

2. emcee (a mcmc tool)
http://dan.iel.fm/emcee/current/

3. cosmolopy
http://roban.github.io/CosmoloPy/

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


Installation
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
1. Unpack the tarball lfit_gui.tar in a directory of your choice (hereby refer to INSTALL_DIR).

mv lfit_gui.tar INSTALL_DIR
cd INSTALL_DIR
tar -xvf lfit_gui.tar

A directory named lfit_gui will be created. You are free to remove the tar file. 

2. Set up the environmental variables (use bash as an example, if using c shell, modify accordingly).

export LFIT_GUI_DIR=INSTALL_DIR/lfit_gui/
export PYTHONPATH=$LFIT_GUI_DIR:$PYTHONPATH

3. Create aliases for the executables 
alias lfit_gui='INSTALL_DIR/lfit_gui/lfit_gui'
alias lfit_script='INSTALL_DIR/lfit_gui/lfit_script'
alias lfit_generate_biz='INSTALL_DIR/lfit_gui/lfit_generate_biz'
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


Set up data structure
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
lfit_gui is originally designed to work with HST data, and therefore requires the data to be organized in a specific manner. 

1. Create HST_DIR
You'll need a root directory for all the data (hereby refer to HST_DIR)

mkdir HST_DIR
export HST_DIR=HST_DIR
export HST_DATAROOT=HST_DIR

2. Make cycle/program/visit/ subdirectories 
