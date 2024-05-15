#!/bin/bash

cd data/temperature/

if [ -f "2t_day_era5-land_$1.nc" ]
then
	echo "file $1 exists"
else
	unzip -jq temp.zip "temperature/2t_day_era5-land_$1.nc" -d "."
	echo "Warning: file $1 DOESNOT exist"
fi

