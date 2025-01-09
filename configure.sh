#!/bin/bash
case $1 in
"clean")
	rm -r ~/.R/Makevars
	;;
"cleanProject")
	rm -r .Rproj.user/*
	;;
"copy2Git")
	rsync -a --exclude 'rmd' . '/mnt/d/work/202301/Aerith for github'
	;;
"debug")
	mkdir ~/.R
	rm -r ~/.R/Makevars
	touch ~/.R/Makevars
	echo "MAKEFLAGS = -j8" >>~/.R/Makevars
	;;
*)
	mkdir ~/.R
	rm -r ~/.R/Makevars
	touch ~/.R/Makevars
	echo "CXXFLAGS= -O3 -std=c++14 -Wall -march=native" >~/.R/Makevars
	echo "CXX14FLAGS = -O3 -std=c++14 -Wall -march=native" >>~/.R/Makevars
	echo "CXX11FLAGS = -O3 -std=c++11 -Wall -march=native" >>~/.R/Makevars
	echo "CXX17FLAGS = -O3 -std=c++17 -Wall -march=native" >>~/.R/Makevars
	echo "MAKEFLAGS = -j8" >>~/.R/Makevars
	;;
esac
