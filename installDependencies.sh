#!/bin/bash

# Check if the operating system is Linux

if [[ "$(uname)" == "Linux" ]]; then
    echo "Found Linux OS"
    echo "Installation with proceed ..."
        printf "\n\n"
else
    echo "Script works only on Linux"
    exit 1;
fi

# Ask user if it is okay to install some libraries
# in the root directory (/usr/lib) and (/usr/include)

read -p "Dependencies will be installed to /usr directory. Do you want to proceed? (Y/N): " choice

# Check if the user wants to proceed
if [[ $choice == "y" ]]; then
    
    echo "Installing dependencies: build-essential g++ zlib1g-dev libreadline-dev freeglut3-dev"
    sudo apt-get install -y git build-essential g++ zlib1g-dev libreadline-dev freeglut3-dev
    
    printf "\n\n"
    echo "Downloading Jem and Jive, version 3.0"
    cd ..
    git clone https://github.com/ritukeshbharali/jemjive-3.0.git
    
    read -p "Do you want to install the MPI version? (Y/N): " choice
    
    if [[ $choice == "y" ]]; then
        
        sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev
        
        cd ./jemjive-3.0/jem-3.0 && ./configure --cxx=mpic++ && make -j5
        printf "\n\n"
        echo "Built MPI jem library"
        echo $(pwd)    
        export JEMDIR=$(pwd)
        
    else
    
        cd ./jemjive-3.0/jem-3.0 && ./configure && make -j5
        printf "\n\n"
        echo "Built Multithreaded jem library"
        echo $(pwd)    
        export JEMDIR=$(pwd)
        
   fi
            
    cd ../jive-3.0 && ./configure && make -j5
    echo $(pwd)
    export JIVEDIR=$(pwd)
    cd ../../
    
    printf "\n\n"
    echo "Export these environment variables for future simulations:"
    echo "export JEMDIR='$JEMDIR'"
    echo "export JIVEDIR='$JIVEDIR'"
    
else

   exit 1;    
    
fi    
    
    
    
