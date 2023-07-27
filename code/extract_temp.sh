#!/bin/bash
cd ~/Desktop/ExtremePrecip/
START=$1;
END=$2

##conda activate RVecchia
echo $START $END
for ((i=$START; i<=$END; i++))
do
	if [ ! -f "/temperature_$i.RData" ]
	then 
	echo $i 
	bash code/unzip_file.sh $i

	cd ~/Desktop/ExtremePrecip/

	Rscript code/extract_temperature.R "year=$i" 

	rm ~/Desktop/Temperatures/2t_day_era5-land_${i}.nc
	fi
done