#!/bin/bash

# to call this shell script, use this, from directory [project]/experiments/timing_experiments:
# ./path_v_fullgrow_exps.sh
clear
echo 'ppr-path v pprgrow experiment2'

# change permissions to ensure this script can run everything it's supposed to
# chmod 744 *.sh
# chmod 744 *.m


nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'itdk0304' , 1e-5, 100, 0);" > pa-td.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'dblp' , 1e-5, 100, 0);" > pa-db.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'fb-one' , 1e-5, 100, 0);" > pa-one.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'fbA' , 1e-5, 100, 0);" > pa-fbA9.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'youtube' , 1e-5, 100, 0);" > pa-yout.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'hollywood-2009' , 1e-5, 100, 0);" > pa-hly.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'twitter-2010' , 1e-5, 100, 0);" > pa-tw.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'friendster' , 1e-5, 100, 0);" > pa-fr.txt &
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'ljournal-2008' , 1e-5, 100, 0);" > pa-lj.txt &


# RHO = 0.9

#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'itdk0304' , 1e-5, 100, 0.9);" > pa-td.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'dblp' , 1e-5, 100, 0.9);" > pa-db.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'fb-one' , 1e-5, 100, 0.9);" > pa-one.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'fbA' , 1e-5, 100, 0.9);" > pa-fbA9.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'youtube' , 1e-5, 100, 0.9);" > pa-yout.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'hollywood-2009' , 1e-5, 100, 0.9);" > pa-hly.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'twitter-2010' , 1e-5, 100, 0.9);" > pa-tw.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'friendster' , 1e-5, 100, 0.9);" > pa-fr.txt &
#nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r "path_v_fullgrow( 'ljournal-2008' , 1e-5, 100, 0.9);" > pa-lj.txt &

echo 'finished calling all experiments'