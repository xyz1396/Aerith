#!/bin/bash
case $1 in
"clean")
	rm -r ~/.R
	;;
"cleanProject")
	rm -r .Rproj.user/*
	;;
*)
	mkdir ~/.R
	touch ~/.R/Makevars
	echo "CXXFLAGS= -O3 -std=c++14 -Wall -march=native" >~/.R/Makevars
	echo "CXX14FLAGS = -O3 -std=c++14 -Wall -march=native" >>~/.R/Makevars
	echo "CXX11FLAGS = -O3 -std=c++11 -Wall -march=native" >>~/.R/Makevars
	echo "CXX17FLAGS = -O3 -std=c++17 -Wall -march=native" >>~/.R/Makevars
	echo "MAKEFLAGS = -j4" >>~/.R/Makevars
	;;
esac
