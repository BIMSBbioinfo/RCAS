#!/bin/bash

autoreconf -vif

./configure MEME_CHIP=meme-chip   --prefix=$PREFIX

make install
