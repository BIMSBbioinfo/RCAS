#!/bin/bash

#the recipe follows the instruction in
#https://github.com/jgm/pandoc/blob/master/INSTALL
#sections: Custom install & Creating a relocatable binary

#the binary is compiled on Ubuntu 14.04 trusty
#using ghc 7.10.1, cabal 1.22.6.0, hsb2hs-0.3.1

#to obtain ghc 7.10.1 and cabal 1.22 on Ubuntu 14.04:
#sudo add-apt-repository ppa:hvr/ghc
#sudo apt-get update
#sudo apt-get install ghc-7.10.1
#sudo apt-get install cabal-install-1.22
#export PATH=/opt/ghc/7.10.1/bin:$PATH
#export PATH=/opt/cabal/1.22/bin:$PATH

#to obtain hsb2hs 0.3.1:
#https://hackage.haskell.org/package/hsb2hs-0.3.1/hsb2hs-0.3.1.tar.gz
#cabal install hsb2hs.cabal
#export PATH=$HOME/.cabal/bin:$PATH

cabal update
cabal install --only-dependencies

#Creating a relocatable binary
#The resulting binary can be run from any directory
#and is completely self-contained.
cabal install hsb2hs  # required for --flags="embed_data_files"
cabal configure --flags="embed_data_files" #embed all data files into the binary
cabal build

mkdir -p $PREFIX/bin

cp $SRC_DIR/dist/build/pandoc/pandoc $PREFIX/bin/pandoc
