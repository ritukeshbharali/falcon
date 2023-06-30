#!/bin/bash

echo "******************************"  | tee    tests.log
echo "   Running tests for falcon"     | tee -a tests.log
echo "******************************"  | tee -a tests.log
echo " "                               | tee -a tests.log

VAR1="Elapsed"
num=$(find . -mindepth 1 -type d | wc -l)
i=0
passd=0
faild=0

for d in */ ; do
    cd $d
    i=$((i+1))
    echo "------------------------------------------------------- " | tee -a ../tests.log
    echo " Running problem $i of $num"                              | tee -a ../tests.log
    echo " "                                                        | tee -a ../tests.log
    echo " $(sed -n '2,7p' < problem.pro)"                          | tee -a ../tests.log
    echo " "                                                        | tee -a ../tests.log    
    
    ../../build/falcon-opt problem.pro
    line=($(tail -n -2 problem.log))
    line1=($(tail -n -4 problem.log))
    
    if [ "$line" = "$VAR1"  ] || [ "$line1" = "$VAR1"  ]; then
    	passd=$((passd+1))
    	echo " STATUS: Passed                           " | tee -a ../tests.log
    	echo " "                                          | tee -a ../tests.log
    else
    	faild=$((faild+1))
    	echo " STATUS: Failed                           " | tee -a ../tests.log
    	echo " "                                          | tee -a ../tests.log  
    fi

    cd ..
done

echo " "                                                        | tee -a ./tests.log    
echo "******************************************************* " | tee -a ./tests.log
echo " SUMMARY:  "                                              | tee -a ./tests.log
echo " Ran    $i of $num problems"                              | tee -a ./tests.log
echo " Passed $passd "                                          | tee -a ./tests.log
echo " Failed $faild "                                          | tee -a ./tests.log
echo "******************************************************* " | tee -a ./tests.log
