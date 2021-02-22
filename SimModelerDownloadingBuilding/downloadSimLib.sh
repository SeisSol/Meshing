#!/bin/bash -ex

EXPECTED_ARGS=4

if [ $# -ne $EXPECTED_ARGS ];then
  echo $#
  echo "Usage: `basename $0` <username> <password> <downloadcode (0 or 1): 0: library 1:documentation><release string ex) 7.1-110613 > "
  exit 1
fi;

username=$1
password=$2
downloadcode=$3
release=$4

#is the release a dev release?
str=$4
i=$((${#str}-3))
last3char="${str:$i:3}"

if [ "$last3char" = "dev" ]; then 
    fold='DM' 
else
    fold='M'  
fi;

# in vi add newlines
# %s/<\/a>/<\/a>\r/g
# to source.php( source code of support download site ) then run the command 
# grep documentation.*zip ~/Downloads/source.php | sed 's/.*documentation\/\([A-Za-z]*.zip\).*/\1/g'
# to extract the zip file names


if [ "$downloadcode" = "1" ]; then 
    documentation=( GeomSim.zip GeomSimAcis.zip GeomSimDiscrete.zip FieldSim.zip GeomSimAbstract.zip MeshSimCrack.zip MeshSimAdapt.zip MeshSimAdv.zip MeshSim.zip MeshSimCrack.zip ParallelMeshSimAdapt.zip ParallelMeshSim.zip GeomSimProe.zip GeomSimParasolid.zip GeomSimSolidWorks.zip )
    echo $documentation
else
    documentation=()
fi;

for doc in ${documentation[@]}; do 
  wget --user=${username} --password=${password} http://www.simmetrix.com/application/release/${fold}/${release}/documentation/${doc}
done

# run the command
# grep linux64.tgz ~/Downloads/source.php | sed 's/.*release\/\([a-z]*-linux64.tgz\).*/\1/g'
# to extract the tarball names

if [ "$downloadcode" = "0" ]; then 
    components=( gmcore-linux64.tgz aciskrnl-linux64.tgz discrete-linux64.tgz fdcore-linux64.tgz gmabstract-linux64.tgz gmadv-linux64.tgz msadapt-linux64.tgz msadv-linux64.tgz mscore-linux64.tgz msparalleladapt-linux64.tgz msparallelmesh-linux64.tgz pskrnl-linux64.tgz )
else
    components=()
fi;

for comp in ${components[@]}; do 
  wget --user=${username} --password=${password} http://www.simmetrix.com/application/release/${fold}/${release}/release/${comp}
done

