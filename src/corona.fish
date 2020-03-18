#! /usr/bin/fish
cd ~/data/COVID-19/
git pull
cd ~/prog/corona/src/
julia ./corona.jl
evince ../figs
