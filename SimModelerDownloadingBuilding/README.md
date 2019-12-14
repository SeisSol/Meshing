# Downloading SimModeler

Go to [the support area of Simmetrix web site](http://www.simmetrix.com/index.php/support/support-downloads) and download the binaries of the code.

# Customizing SimModeler for SeisSol
Add the files of [this folder](https://github.com/SeisSol/Meshing/tree/master/SimModelerDownloadingBuilding/SimModelerCustomization) to the main directory of SimModeler. These files allow defining properly boundary conditions and exporting the mesh in the proper format.


# Downloading SimModeler modeling suite

Go to [the support area of Simmetrix web site](http://www.simmetrix.com/index.php/support/support-downloads) to see the available versions of SimModeler modeling suite. Then use this [script](https://github.com/SeisSol/Meshing/tree/master/SimModelerDownloadingBuilding/downloadSimLib.sh) to download all components of the library. Please make sure that all components that you download are included in your licence agreeement.

# Building SimModeler modeling suite
First untar all the files, using for instance:
```bash
for filename in *.tgz
do
  tar zxf $filename
done
```
Then adapt the Makefile in  `code/PartitionWrapper` to your plateform (basically update variables PQUAL (e.g. -mpich2), CC (e.g. mpicc) and CXX (e.g. mpiCC)). As this Makefile usually do not change from one version to the other, it can be wise to copy it to the root folder of your SimModeler installation path, and use this customized file to build any version of the SimModeler library. In this case, the library is then built using:
`sh compileSimLib.sh <release string ex) 7.1-110613 >` where compileSimLib.sh is this file:


```bash
#!/bin/sh
cp Makefile $1/code/PartitionWrapper/
cd  $1/code/PartitionWrapper/
make PARALLEL=mpich2
cp lib/libSimPartitionWrapper-*.a ../../lib/x64_rhel6_gcc44/
cd ../../..
```
