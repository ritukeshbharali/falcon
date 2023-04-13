#!/bin/bash

VAR1="Elapsed"
num=$(find . -mindepth 1 -type d | wc -l)
i=0
passd=0
faild=0

for d in */ ; do
    cd $d
    i=$((i+1))
    echo "******************************************************* "
    echo " Running problem $i of $num"
    echo " $(head -n 1 problem.pro)"
    echo "******************************************************* "
    ../../build/falcon-opt problem.pro
    line=($(tail -n -2 problem.log))
    line1=($(tail -n -4 problem.log))
    
    if [ "$line" = "$VAR1"  ] || [ "$line1" = "$VAR1"  ]; then
    	passd=$((passd+1))
    	echo "*-----------------------------------------------------* "
    	echo "                       Passed                           "
    	echo "*-----------------------------------------------------* "
    	echo " "
    else
    	faild=$((faild+1))
    	echo "*-----------------------------------------------------* "
    	echo "                       Failed                           "
    	echo "*-----------------------------------------------------* "
    	echo " "
    fi

    cd ..
done

echo "******************************************************* "
echo " Ran    $i of $num problems"
echo " Passed $passd "
echo " Failed $faild "
echo "******************************************************* "