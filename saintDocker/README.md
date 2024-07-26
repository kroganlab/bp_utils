
# saintDocker

Run SAINTexpress PPI scorer on any machine with docker.  (This has only been tested on a few Apple notebook computers.). This directory contains the files necessary to create the docker image, and a shell script to help you run saintExpress from the image. 

To run Docker you will first build the image locally, then use the saintDocker.sh script to help run SAINT within the docker iamge.


## Usage:

### option 1) build the docker image
```
# build the docker image within the  saintDocker directory that you downloaded from github
cd saintDocker
docker build -t bpolacco/saint . 
# this will run for a few minutes
```

### option 2) (easiest) use pre-built docker image from docker hub
There is a prebuilt image available at docker hub.  You can download it, and then change its name (tag) with:

```
docker pull martingordon808/saint_bp

docker image tag bpolacco/saint  martingordon808/saint_bp
```


### verify the build (optional)
If it worked you should see the new image in your list of docker images from the command `docker images`:
```
$ docker images
REPOSITORY                                TAG       IMAGE ID       CREATED          SIZE
bpolacco/saint                            latest    a012accb19cb   10 minutes ago   1.2GB
...
```

To  further test if it worked you can open a BASH shell within a disposable container and verify the tools are available.  The following ouput is from running these commands (but once in the shell feel free to look around with any BASH/linux commands):

```
docker run -it --rm -v`pwd`:/wd -w /wd bpolacco/saint  /bin/bash
SAINTexpress-spc
SAINTexpress-int
exit
```



```
$ docker run -it --rm -v`pwd`:/wd -w /wd bpolacco/saint  /bin/bash
root@6c8077997957:/wd# SAINTexpress-spc
Input files not supplied, using defaults.
Input files are: inter.dat, prey.dat, bait.dat
Interaction file: "inter.dat"
Prey file: "prey.dat"
Bait file: "bait.dat"
GO file: ""
Parsing prey file prey.dat ...terminate called after throwing an instance of 'std::runtime_error'
  what():  invalid delimiter
Aborted
root@6c8077997957:/wd# SAINTexpress-int
Input files not supplied, using defaults.
Input files are: inter.dat, prey.dat, bait.dat
Interaction file: "inter.dat"
Prey file: "prey.dat"
Bait file: "bait.dat"
GO file: ""
Parsing prey file prey.dat ...terminate called after throwing an instance of 'std::runtime_error'
  what():  invalid delimiter
Aborted
root@6c8077997957:/wd# exit
exit
$
```


### run SAINTexpress using shell script
Best is to start in a directory that contains your bait,prey,interactions file output from artMS, and then use full path to the shell script in the saintDocker folder.

```
cd msint
../../../saintDocker/saintDocker.sh int  my_interactions.txt my_preys.txt my_baits.txt
```

Output is lists.txt in current directory
