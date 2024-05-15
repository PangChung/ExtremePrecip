#!/bin/bash
cd ~/data/temperature/

tar -czvf temperature1950_1980.tar.gz temperature/

rm temperature/*
for i in {1981..2000}
do
bash unzip_file.sh $i
done

tar -czvf tempearture1981_2000.tar.gz temperature/
rm temperature/*
for i in {2001..2021}
do
bash unzip_file.sh $i
done
tar -czvf tempearture2001_2021.tar.gz temperature/
rm temperature/*

