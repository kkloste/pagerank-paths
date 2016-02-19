#!/bin/bash

# to call this shell script, use this, from directory [project]/experiments/timing_experiments:
# ./rho-scaling-exps.sh
clear
echo 'Rho scaling experiment starting'

# change permissions to ensure this script can run everything it's supposed to
#chmod 744 *.sh
#chmod 744 *.m

nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling('itdk0304', 100, 50);" > rs-td.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'dblp', 100, 50);" > rs-db.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'fb-one', 100, 50);" > rs-one.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'fbA', 100, 50);" > rs-fbA9.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'youtube', 100, 50);" > rs-yout.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'hollywood-2009', 100, 50);" > rs-hly.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'twitter-2010', 100, 50);" > rs-tw.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'friendster', 100, 50);" > rs-fr.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "rho_scaling( 'ljournal-2008', 100, 50);" > rs-lj.txt &

echo 'finished calling all experiments'