# Downloading SimModeler

Go to [the support area of Simmetrix web site](http://www.simmetrix.com/index.php/support/support-downloads) and download the binaries of the code.

# Customizing SimModeler for SeisSol
Add the files of [this folder](https://github.com/SeisSol/Meshing/tree/master/SimModelerDownloadingBuilding/SimModelerCustomization) to the main directory of SimModeler. These files allow defining properly boundary conditions and exporting the mesh in the proper format.


# Downloading SimModeler modeling suite

Go to [the support area of Simmetrix web site](http://www.simmetrix.com/index.php/support/support-downloads) to see the available versions of SimModeler modeling suite. Then use this [script](https://github.com/SeisSol/Meshing/tree/master/SimModelerDownloadingBuilding/downloadSimLib.sh) to download all components of the library. Please make sure that all components that you download are included in your licence agreeement.

Finally, untar all the files, using for instance:
```bash
for filename in *.tgz
do
  tar zxf $filename
done
```
The libraries are precompiled, so no need to build them!
