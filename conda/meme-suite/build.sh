#!/bin/sh

./configure --prefix=`pwd`/meme --with-url="http://meme-suite.org"

make install

cp -R ./meme/bin/ $PREFIX
