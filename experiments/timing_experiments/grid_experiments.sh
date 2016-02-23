#!/bin/bash

# to call this shell script, use this, from directory [project]/experiments/timing_experiments:
# ./grid_experiments.sh
clear
echo 'grid_experiments'

# change permissions to ensure this script can run everything it's supposed to
# chmod 744 *.sh
# chmod 744 *.m

echo 'grid_experiments runs'

nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'itdk0304' , 100);" > grid-itdk.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'dblp' , 100);" > grid-dblp.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'fb-one' , 100);" > grid-one.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'fbA' , 100);" > grid-fbA.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'youtube' , 100);" > grid-you.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'ljournal-2008' , 100);" > grid-lj.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'twitter-2010' , 100 );" > grid-tw.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'hollywood-2009' , 100);" > grid-hly.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_new( 'friendster' , 100);" > grid-fr.txt &

echo 'finished calling all experiments'