#!/bin/bash

# to call this shell script, use this, from directory [project]/experiments/timing_experiments:
# ./grid_experiments.sh
clear
echo 'grid_experiments'

# change permissions to ensure this script can run everything it's supposed to
# chmod 744 *.sh
# chmod 744 *.m

echo 'grid_experiments runs'

nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'itdk0304' , 100);" > c-11-29-itdk.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'dblp' , 100);" > c-11-29-dblp.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'fb-one' , 100);" > c-11-29-one.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'fbA' , 100);" > c-11-29-fbA.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'youtube' , 100);" > c-11-29you.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'ljournal-2008' , 100);" > c-11-29-lj.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'twitter-2010' , 100 );" > ~/ppr-all-eps/experiments/c-11-29-tw.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'hollywood-2009' , 100);" > ~/ppr-all-eps/experiments/c-11-29-hly.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'friendster' , 100);" > ~/ppr-all-eps/experiments/c-11-29-fr.txt &

echo 'finished calling all experiments'