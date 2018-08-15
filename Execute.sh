#!/usr/bin/bash

startt=`date +%s`

compile_sources=true
#compile_sources=false
generate_data=true
#generate_data=false
#fast_marching=true
fast_marching=false
migration=true
#migration=false
slicing=true
#slicing=false
plotting=true
#plotting=false
#qsub_batch=true
qsub_batch=false

echo
echo
echo '              +  +++++  +++++ + +++ ++ +++ + +++++  +++++  +'
echo '            ++  ++++                                  ++++  ++'
echo '          +++  +++                                      +++  +++'
echo '        ++++  ++                                          ++  ++++'
echo '      +++++  +        RUNNING THE FULLY 3D MIGRATION        +  +++++'
echo '        ++++  ++                                          ++  ++++'
echo '          +++  +++                                      +++  +++'
echo '            ++  ++++                                  ++++  ++'
echo '              +  +++++  +++++ + +++ ++ +++ + +++++  +++++  +'
echo

echo
echo
echo
echo '    1/ Compiling the sources .................... '$compile_sources
echo '    2/ Generating the data files ................ '$generate_data
echo '    3/ Performing the fast marching ............. '$fast_marching
echo '    4/ Performing the migration ................. '$migration
echo '    5/ Slicing throught the results ............. '$slicing
echo '    6/ Plotting the scattering maps ............. '$plotting
echo
echo '    >> Doing all this on the server ............. '$qsub_batch

echo
echo
echo
echo '    The following files need to be edited for different runs:'
echo '         - Input/GenData.in         to switch modes'
echo '         - Run/fm3d/execute_fm3d.sh for number of sources and receivers'
echo '         - Input/GenData.in         for the synthetic geometries'
echo '         - Run/raysum/mod           for the synthetic sharp model'
echo '         - Source/fm3d/vdefs.f90    for the 2D smooth model'
echo '         - Source/fm3d/idefs.f90    for the lower limit of the model'
echo '         - Run/fm3d/propgrid.in     for the dimensions of the computation grid'
echo '         - Input/mk3d.in            for the dimensions of the migration grid'
echo '         - Plot/execute_plots.sh    to change the output names'

now=`date`
echo
echo
echo
echo
echo '                      ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~'
echo '                      '$now
echo '                      ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~'
echo
echo
echo


if [ "$compile_sources" = true ] ; then
  echo
  echo '                     +++++++++++++++++++++++++++++++'
  echo '                  +++++++++++++++++++++++++++++++++++++ '
  echo '                ++++                                 ++++'
  echo '               +++                                     +++'
  echo '               ++          Compiling the codes          ++'
  echo '               +++                                     +++'
  echo '                ++++                                 ++++'
  echo '                  +++++++++++++++++++++++++++++++++++++'
  echo '                     +++++++++++++++++++++++++++++++'
  echo
  
  cd Source/
  make gendata
  make slice
  make mk3d
  make move
  make clean

  cd fm3d/
  make mvel
  make mint
  make fm3d
  make move
  make clean

  cd ../../
fi
  

if [ "$generate_data" = true ] ; then
  echo
  echo '                  ///\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\\\'
  echo '                  ||                                  ||'
  echo '                  ||         Generating data          ||'
  echo '                  ||                                  ||'
  echo '                  \\\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\///'
  echo

  cd Run/
  ./gendata

  cd ../
fi


if [ "$fast_marching" = true ] ; then
  echo
  echo '               o o o o o o o o o o o o o o o o o o o o o o'
  echo '              o                                           o'
  echo '             o                                             o'
  echo '            o           Running the Fast Marching           o'
  echo '             o                                             o'
  echo '              o                                           o'
  echo '               o o o o o o o o o o o o o o o o o o o o o o'
  echo

  cd Run/fm3d/
  ./mint
  ./mvel
  if [ "$qsub_batch" = true ] ; then
    qsub batch_fm3d
  else
    ./execute_fm3d.sh
  fi

  cd ../../
fi
  

if [ "$migration" = true ] ; then
  echo
  echo '                 + - - - - - - - - - - - - - - - - - - +'
  echo '                 |                                     |'
  echo '                 |         Migrating the data          |'
  echo '                 |                                     |'
  echo '                 + - - - - - - - - - - - - - - - - - - +'
  echo
  
  cd Run/
  if [ "$qsub_batch" = true ] ; then
    qsub batch_mk3d
  else
    ./mk3d
  fi

  cd ../
fi
    
  
if [ "$slicing" = true ] ; then
  latA="$(awk 'NR==09 {print $1}' Input/Slice.in)"
  lonA="$(awk 'NR==10 {print $1}' Input/Slice.in)"
  latB="$(awk 'NR==13 {print $1}' Input/Slice.in)"
  lonB="$(awk 'NR==14 {print $1}' Input/Slice.in)"
  echo
  echo '            +----------------------------------------------+'
  echo '            | +------------------------------------------+ |'
  echo '            | |                                          | |'
  echo '            | |     Slice from point A('$latA','$lonA')      | |'
  echo '            | |             to point B('$latB','$lonB')      | |'
  echo '            | |                                          | |'
  echo '            | +------------------------------------------+ |'
  echo '            +----------------------------------------------+'
  echo
  
  cd Run/
  ./slice

  cd ../
fi

  
if [ "$plotting" = true ] ; then
  echo
  echo '                  + + ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ + +'
  echo '                   \ \ ~~ ~~ ~~ ~~ ~~~ ~~ ~~ ~~ ~~ / /'
  echo '                    \ \                           / /'
  echo '                     \ \        Starting         / / '
  echo '                     / /        the plots!       \ \ '
  echo '                    / /                           \ \'
  echo '                   / / ~~ ~~ ~~ ~~ ~~~ ~~ ~~ ~~ ~~ \ \'
  echo '                  + + ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ + +'
  echo
  
  cd Plot/
  ./execute_plots.sh
  #./Depth_map.gmt

  cd ../
fi


now=`date`
echo
echo
echo
echo
echo '                      ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~'
echo '                      '$now
echo '                      ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~'


endd=`date +%s`
runtime=$((endd-startt))
runhou=$((  runtime/3600 ))
runmin=$(( (runtime - (runhou*3600)) / 60 ))
runsec=$((  runtime - (runhou*3600) - (runmin*60) ))


echo
echo
echo '                    Total runtime:'
echo '                    '$runtime' seconds'
echo
echo '                    Total runtime:'
echo '                    '$runhou' hours, '$runmin' minutes and '$runsec' seconds'
echo
echo
