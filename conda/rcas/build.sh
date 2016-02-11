#!/bin/bash

autoreconf -vif

./configure --prefix=$PREFIX

make install
