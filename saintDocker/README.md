
# saintDocker

Run SAINTexpress PPI scorer on any machine with docker.  (This has only been tested on a few Apple notebook computers.). This directory contains the files necessary to create the docker image, and a shell script to help you run saintExpress from the image. Because the docker image is already created and shared on DockerHub, you can get by with just the shell script, or none of the files if you want to manage `docker run` yourself on the command line. 

## Usage:

```
#best is to start in a directory that contains your bait,prey,interactions file output from artMS...

cd msint
../../../saintDocker/saintDocker.sh int  my_interactions.txt my_preys.txt my_baits.txt
```
If docker is installed and this is the first time you've run this, it may take a minute to download the docker image.

Output is lists.txt in current directory
