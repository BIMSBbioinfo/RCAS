#!/bin/sh

#transform py2 to py3
#so there is no conflict beween
2to3 -w scripts/*py
2to3 -w scripts/*py.in

./configure --prefix=`pwd`/meme --with-url="http://meme-suite.org"

make install

cp -R ./meme/bin/ $PREFIX
