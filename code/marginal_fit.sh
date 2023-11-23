#!/bin/bash
cd ~/Desktop/ExtremePrecip/
for i in {1..8}
do Rscript code/marginal_fit.R "idx=$i"
done

echo "Your program has been finished" | mail -s "Your program has been finished" -a "From:zhongprw@gmail.com" -s smtp="smtp.gmail.com:587" -s smtp-use-starttls -s smtp-auth=login -s smtp-auth-user="zhongprw@gmail.com" -s "Marginal Fit" peng.zhong@unsw.edu.au 