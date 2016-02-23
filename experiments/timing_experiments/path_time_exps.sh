#!/bin/bash

# to call this shell script, use this, from directory [project]/experiments/timing_experiments:
# ./path_time_exps.sh
clear
echo 'path time experiment'

# change permissions to ensure this script can run everything it's supposed to
# chmod 744 *.sh
# chmod 744 *.m

nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'itdk0304' , 1e-5, 100, 0);" > pat-td.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'dblp' , 1e-5, 100, 0);" > pat-db.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'hollywood-2009' , 1e-5, 100, 0);" > pat-hly.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'fb-one' , 1e-5, 100, 0);" > pat-one.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'fbA' , 1e-5, 100, 0);" > pat-fbA9.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'youtube' , 1e-5, 100, 0);" > pat-you9.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'twitter-2010' , 1e-5, 100, 0);" > pat-tw.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'friendster' , 1e-5, 100, 0);" > pat-fr.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'ljournal-2008' , 1e-5, 100, 0);" > pat-lj.txt &

# RHO = 0.9
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'itdk0304' , 1e-5, 100, 0.9);" > pat-td.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'dblp' , 1e-5, 100, 0.9);" > pat-db.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'hollywood-2009' , 1e-5, 100, 0.9);" > pat-hly.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'fb-one' , 1e-5, 100, 0.9);" > pat-one.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'fbA' , 1e-5, 100, 0.9);" > pat-fbA9.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'youtube' , 1e-5, 100, 0.9);" > pat-you9.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'twitter-2010' , 1e-5, 100, 0.9);" > pat-tw.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'friendster' , 1e-5, 100, 0.9);" > pat-fr.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "do_path_time( 'ljournal-2008' , 1e-5, 100, 0.9);" > pat-lj.txt &

echo 'finished calling all experiments'