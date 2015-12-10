#!/bin/bash
# script for execution of CLImAT-HET algorithm
# 10/12/2015 by Zhenhua, yzh163@mail.ustc.edu.cn

echo 

function usage() {
	echo -e "$0 <MCRDIR> <dataDir> <outputDir> <plotDir> <configFile>\n"
}

if [ $# -ne 5 ]; then
	echo "5 parameters are required!"
	echo -n "Usage: " 
	usage
	exit 0
fi

MCRDIR=$1
dataDir=$2
outputDir=$3
plotDir=$4
configFile=$5

MCRDIR=${MCRDIR%/}
dataDir=${dataDir%/}
outputDir=${outputDir%/}
plotDir=${plotDir%/}

[ ! -e $dataDir ] && echo "Error: Data directory $dataDir does not exist!" && exit 0
[ ! -e $configFile ] && echo "Error: Configuration file $configFile does not exist!" && exit 0

[ ! -e $outputDir ] && echo "Create the directory $outputDir!" && mkdir -p $outputDir
[ ! -e $plotDir ] && echo "Create the directory $plotDir!" && mkdir -p $plotDir

echo "Setting up environment variables"
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRDIR}/runtime/glnxa64:${MCRDIR}/bin/glnxa64:${MCRDIR}/sys/os/glnxa64:${MCRDIR}/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:${MCRDIR}/sys/java/jre/glnxa64/jre/lib/amd64/server:${MCRDIR}/sys/java/jre/glnxa64/jre/lib/amd64
#LD_LIBRARY_PATH=.:${MCRDIR}/runtime/glnxa64
#LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRDIR}/bin/glnxa64
#LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRDIR}/sys/os/glnxa64
XAPPLRESDIR=${MCRDIR}/X11/app-defaults
export XAPPLRESDIR
export LD_LIBRARY_PATH

eval ./CLImATHET $dataDir $outputDir $configFile
eval ./CLImATHET_plot $dataDir $outputDir $plotDir

exit 0
