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
skipd=0

# Function to check if CUDA is available
check_cuda() {
    if ! command -v nvidia-smi &> /dev/null; then
        return 1
    fi
    return 0
}

for d in */ ; do

    # Remove trailing slash to get the folder name
    dir_name=${d%/}
    
    # Check if the folder name starts with "cu" and if CUDA is required
    if [[ "$dir_name" == cu* ]]; then
        # Only run the test if CUDA is available
        if ! check_cuda; then
            echo "Skipping directory $dir_name because CUDA is not available." | tee -a tests.log
            skipd=$((skipd+1))
            continue
        fi
    fi
    
    cd "$d"
    rm -f *.log *.dat *.pvd *.vtu
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
echo " Ran     $i of $num problems"                              | tee -a ./tests.log
echo " Passed  $passd "                                          | tee -a ./tests.log
echo " Failed  $faild "                                          | tee -a ./tests.log
echo " Skipped $skipd "                                          | tee -a ./tests.log
echo "******************************************************* " | tee -a ./tests.log
