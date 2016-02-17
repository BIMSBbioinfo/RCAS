#!/bin/sh

#transform py2 to py3
#so there is no conflict, during conda installation,
#between meme-chip and other py3 based programs in RCAS pipeline.
#https://github.com/BIMSBbioinfo/RCAS
#though it is an ad-hoc feature, it should not affect
#the general use of meme-suite.
2to3 -w scripts/*py
2to3 -w scripts/*py.in

./configure --prefix=$PREFIX --with-url="http://meme-suite.org"

make install
