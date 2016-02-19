#!/bin/bash

# to call this shell script, use this, from directory [project]/experiments/timing_experiments:
# ./path_time_exps.sh
clear
echo 'path time experiment'

# change permissions to ensure this script can run everything it's supposed to
# chmod 744 *.sh
# chmod 744 *.m

echo 'netscience runs'

nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'itdk0304' , 1e-5, 100, 0.9);" > pat-td.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'dblp' , 1e-5, 100, 0.9);" > pat-db.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'hollywood-2009' , 1e-5, 100, 0.9);" > pat-hly.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'fb-one' , 1e-5, 100, 0.9);" > pat-one.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'fbA' , 1e-5, 100, 0.9);" > pat-fbA9.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'youtube' , 1e-5, 100, 0.9);" > pat-you9.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'twitter-2010' , 1e-5, 100, 0.9);" > pat-tw.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'friendster' , 1e-5, 100, 0.9);" > pat-fr.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time_new( 'ljournal-2008' , 1e-5, 100, 0.9);" > pat-lj.txt &

echo 'finished calling all experiments'