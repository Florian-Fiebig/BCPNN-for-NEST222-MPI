#!/bin/sh
#Input: Take module-DATE as input

# compile with the GNU compiler
module swap PrgEnv-cray PrgEnv-gnu
module add python/2.7.13 site-python/2.7
module add nest/2.2.2-py27
#For the correct directory
CURR=$(pwd)
echo "current dir: $CURR"

#Start time watch
START=$(date +%s)

#Get number of processors on the system
#noProcs=$(grep -c 'model name' /proc/cpuinfo) 
noProcs=1

#Source directory
srcDir="$CURR/$1"

#Bootstrap directory, used temoprally. Removed at end of script.
bootstrapDir="$CURR/bootstrap-$1"

#Delete old bootstrap and build directories
echo "Removing old build and bootstrap directories"
rm -r "$CURR/build-$1"
rm -r "$CURR/bootstrap-$1"

#Build directory
buildDir="$CURR/build-$1"

echo "Source dir: $srcDir"
echo "Bootstrap dir: $bootstrapDir"
echo "Build dir: $buildDir"

#Copy source to bootstrap directory
#sudo cp -r $srcDir $bootstrapDir
cp -r $srcDir $bootstrapDir

echo $(pwd)
#Go into bootstrap dir and run bootstrap
cd $bootstrapDir
"$bootstrapDir/bootstrap.sh"

#Move out
cd ..

#Make new build directory, configure and run make, make install and make installcheck
mkdir $buildDir
echo "Entering $buildDir"
cd $buildDir

export NEST_INSTALL_DIR=/pdc/vol/nest/2.2.2-py27
$bootstrapDir"/configure" --with-nest=${NEST_INSTALL_DIR}/bin/nest-config --prefix=$CURR/ 2>&1 | tee "$CURR/mymodule-configure.log"
make -j $noProcs 2>&1 | tee "$CURR/mymodule-make.log"
make -j $noProcs install  

#Stop time watch
END=$(date +%s)
DIFF=$(( $END - $START ))

# Move out
cd ..

#Display script execution time
echo "It took $DIFF seconds"
