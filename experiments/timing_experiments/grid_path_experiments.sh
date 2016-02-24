#!/bin/bash

# to call this shell script, use this, from directory [project]/experiments/timing_experiments:
# ./grid_experiments.sh
clear
echo 'grid_path_experiments'

# change permissions to ensure this script can run everything it's supposed to
# chmod 744 *.sh
# chmod 744 *.m

echo 'grid_path_experiments runs'

nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'itdk0304' , 100);" > gridp-itdk.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'dblp' , 100);" > gridp-dblp.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'fb-one' , 100);" > gridp-one.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'fbA' , 100);" > gridp-fbA.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'youtube' , 100);" > gridp-you.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'ljournal-2008' , 100);" > gridp-lj.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'twitter-2010' , 100 );" > gridp-tw.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'hollywood-2009' , 100);" > gridp-hly.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "grid_v_grow_path( 'friendster' , 100);" > gridp-fr.txt &

echo 'finished calling all experiments'