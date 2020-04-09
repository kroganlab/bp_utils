#! /bin/sh

# A shell script to run SAINTexpress from within a docker image. Docker is
# required to run this script successfully.  The Docker image bpolacco/saint
# is hosted on docker hub, so docker will need to download the image on first
# running of this file. 
#
# Use this to run SAINT-express on either intensity or spectral count data.
#
# For file access inside the image, the current working directory gets
#  mounted in the docker image, so paths to SAINT input files (baits,preys,interactions)
# must be contained within the current working directory, within subfolders is fine.
#
# Output is list.txt in current working directory, with no option to change.  That seems 
# to be a limitation of SAINTexpress.
#
# NOTE: If it appears to run successfully but no list.txt file is created, it is possible
# that SAINTexpress failed with a segmentation fault. I found this happened once when
# baits file included replicates (column 1) that were not in the interactions file.
# You may need to enter the container with a shell and troubleshoot there:
# 
# docker run --rm -v`pwd`:/wd -w /wd -it  bpolacco/saint:0.1 /bin/bash

# another source of errors:  only 1 control run with SAINTExpress-int will
# fail by throwing a const char*



usage ()
{

 echo "usage: ./saintDocker.sh int|spc interactions.txt preys.txt baits.txt"
 echo "   all paths must be descendents of current directory"
 echo '   i.e. `pwd`/$2 must be a valid path name'
 echo "output: list.txt in current directory"

}



if [ "$#" -eq 4 ]; then
 case $1 in
    int)
        saintCommand="SAINTexpress-int"
        ;;
    spc)
        saintCommand="SAINTexpress-spc"
        ;;
    *)
        echo "Error: First argument must be one of int, spc"
        usage
        exit 1
        ;;
 esac
 docker run --rm -v`pwd`:/wd -w /wd bpolacco/saint:0.1 $saintCommand $2 $3 $4
else
 usage
 exit 1
fi
