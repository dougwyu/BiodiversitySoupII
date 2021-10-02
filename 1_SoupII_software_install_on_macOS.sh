#!/bin/bash
set -e
set -u
set -o pipefail

####################################################################################################
# Software installation
# tested on macOS Catalina (zsh)
####################################################################################################
# The script uses open-source software and has been tested on macOS Catalina with zsh shell
# Much of this software can be installed using Homebrew, which is a package manager for macOS.
# The whole installation process can require several hours, depending on internet speed and how much you have previously installed

# Install command line tools for macOS
# installs compilers and other Unix tools for programming
# if you have Xcode installed, you do not need to install CLT (but CLT is much smaller than Xcode)
xcode-select --install

# Install Homebrew for macOS
# go to http://brew.sh and follow the instructions for installing Homebrew on macOS

# Install these software packages using Homebrew
brew update; brew upgrade; brew cleanup # get the latest version homebrew
brew tap brewsci/bio # a "tap" is a repository of "installation formulas" of specialist software, here, bioinformatics
brew install git
brew install coreutils
brew install wget
brew install rename
brew install gnu-sed # (gsed == GNU version of sed == Linux version of sed)
brew install grep # gnu-grep
brew install gawk # gnu-awk
brew install parallel # GNU parallel
brew install seqkit # https://github.com/shenwei356/seqkit
brew install vsearch # https://github.com/torognes/vsearch
brew install seqtk # https://github.com/lh3/seqtk
brew install fastp # http://github.com/OpenGene/fastp
brew install python@2 # python2
brew cleanup # remove unneeded files

# Install env_parallel with these two commands; orig instructions at https://stackoverflow.com/a/27559537
echo '. `which env_parallel.zsh`' >> $HOME/.zshenv
source $HOME/.zshenv

# Install miniconda3 for python3, numpy, and matplotlib
    # Installer package at:
         # https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg
    # Instructions at:
        # https://conda.io/projects/conda/en/latest/user-guide/install/index.html
# after installing, confirm that you are using the correct installation of conda
which conda # should return somthing like /Users/HOME/opt/miniconda3/bin/conda
# A common problem with python installations is previous python installations, which you should remove
# There is help at:  https://conda.io/projects/conda/en/latest/user-guide/install/index.html

# Until Begum is upgraded to python3, the above steps should have installed both python2 and python3, which you can check as follows:
which python2 # /usr/local/bin/python
python2 --version # Python 2.7.17
which python3 # /Users/Negorashi2011/opt/miniconda3/bin/python3
python3 --version # Python 3.8.3

# Install numpy and matplotlib
conda install numpy
conda install matplotlib


# Install usearch11, free binary is 32-bit, available at http://drive5.com/usearch/download.html
# download usearch binary after registering at drive5.com
mv usearch11.0.667_i86osx32 /usr/local/bin # adjust usearch11 pathname as needed
cd /usr/local/bin
ln -s usearch11.0.667_i86osx32 usearch11  # make symbolic link (symlink)
chmod +x usearch11


# For the following software, create a dedicated folder at ~/src/
mkdir ~/src/

# Install Begum
# http://github.com/shyamsg/Begum
cd ~/src
git clone https://github.com/shyamsg/Begum.git
# no building needed because these are python2 scripts

# Install DAMe, for one function: convertToUSearch.py
cd ~/Desktop
git clone https://github.com/shyamsg/DAMe.git
mv ~/Desktop/DAMe /usr/local/bin # on ubuntu will need sudo mv
# no building needed because these are python3 scripts


# The following GUI programs are installed in your Applications folder

# Install Atom:  Combined text editor and terminal
# download binary from https://atom.io # (currently 1.52.0)
    # In Atom, go to Atom/Preferences (or type cmd, (command comma)). This opens the Settings tab.  Click on the "Install" button on the left and type 'terminus' in the 'Search packages' window, and when the package appears, click on the blue install button.
    # Terminus installs a terminal inside Atom. You can open terminal windows by clicking on the + sign in the lower left of the Atom window. (you might need to restart Atom)

    # This allows you to send commands from an Atom text window to a terminal window (just like in RStudio), using 'ctrl-enter'.  In some cases, 'ctrl-enter' does not work. This is because the 'keymap' is missing.
         # Go to Atom/Keymap... This opens a new tab entitled keymap.cson
         # add this text to the end of the keymap.cson file:
'atom-workspace atom-text-editor:not([mini])': 'ctrl-enter': 'terminus:insert-selected-text'
        # Fix adapted from:  https://github.com/platformio/platformio-atom-ide-terminal/issues/67

# Install R:  Statistical analysis
# install from https://cran.r-project.org

# Install RStudio:  GUI to R
# install free RStudio Desktop from https://rstudio.org

# Download MIDORI reference file from:  http://www.reference-midori.info/download.php#
# navigate here on the webpage and download:
# /20180221/MIDORI_UNIQUE_20180221/SINTAX_format/MIDORI_UNIQUE_20180221_COI_SINTAX.fasta.zip');

# Install R packages within R
# Launch RStudio and run these commands in R. This step can take hours the first time
install.packages(c("vegan", "car", "data.table", "here", "tidyverse", "readxl", "cowplot", "lubridate", "patchwork", "sjmisc", "broom", "lme4", "MuMin", "ggeasy", "RColorBrewer", "envDocument", "conflicted","BiocManager","remotes", dependencies = TRUE))

remotes::install_github("grunwaldlab/metacoder")

library(BiocManager)
BiocManager::install(c("GenomicRanges", "Biobase", "IRanges", "AnnotationDbi", "phyloseq")) # install Bioconductor packages

BiocManager::valid() # identify packages that are out of date or unexpected versions. Creates a command for updating out of date bioconductor packages
