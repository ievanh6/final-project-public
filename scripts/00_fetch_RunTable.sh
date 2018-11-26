#!/bin/bash
# 00_fetch_RunTable.sh --- 
# 
# Filename: 00_fetch_RunTable.sh
# Description: Download a RunInfoTable given an NCBI SRA BioProject identifier 
# Author: Naupaka Zimmerman
# Maintainer: same
# Created: Mon Nov 27 21:21:18 2017 (-0800)
# Version: 
# Package-Requires: ()
# Last-Updated: November 27, 2017
#           By: NBZ
#     Update #: 0
# URL: 
# Doc URL: 
# Keywords: 
# Compatibility: bash, linux
# 
# 

# Commentary: 
# 
# run as, eg: `bash 00_fetch_RunTable.sh PRJNA138459`
# 
# 

# Change Log:
# 
# 
# 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
# 
# 

# Code:

# Install NCBI e-utils if needed
if ! [ -x "$(command -v esearch)" ]; then
  echo 'Error: git is not installed. Now installing...'
  cd ~
  /bin/bash
  perl -MNet::FTP -e \
       '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
     $ftp->login; $ftp->binary;
     $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
  gunzip -c edirect.tar.gz | tar xf -
  rm edirect.tar.gz
  builtin exit
  export PATH=$PATH:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"
  ./edirect/setup.sh  
fi

# Download SraRunTable from NCBI given BioProject number
# Idea from: https://www.biostars.org/p/260190/
esearch -db sra -query $1 | efetch -format runinfo

# 
# 00_fetch_RunTable.sh ends here
