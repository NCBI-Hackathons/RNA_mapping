![](https://github.com/NCBI-Hackathons/RNA_mapping/blob/master/images/RNAmappermatchup.png)

* Do you want to learn how to align RNA-sequences to a reference genome? 

* Do you want to bypass the process of installing alignment software and get a preview into a panel of different alignment strategies?

* Do you want to compare how popular mapping algorithms perform on your data?

#### Then this is the tutorial for you!  [Please proceed to our wiki](https://github.com/NCBI-Hackathons/RNA_mapping/wiki)

#### *Just interested in comparing two or more .bam files?  Check out [bamdiff](https://github.com/NCBI-Hackathons/RNA_mapping/tree/master/build/bamdiff)*

#### Docker
A Docker container is available for this project.  For a detailed description of docker, please refer to this [overview.](https://www.docker.com/what-docker)

##### Install Docker
Follow [instructions](https://www.docker.com/docker-mac) to install Docker for your environment.

##### Get a pre-built image from DockerHub and run the server
```
docker pull stevetsa/rna_mapping
docker run -it stevetsa/rna_mapping

At the root prompt inside the container - 
/RNA_mapping/build/doAll.sh
```
