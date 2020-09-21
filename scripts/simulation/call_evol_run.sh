#!/bin/bash

WORK_DIR=""
RESULT_DIR=""

# Setting the INTRO_REPL2 = 0, i.e. sequence appears and immediately disappears for extreme numbers of introductions
mkdir -p  $RESULT_DIR

#p_repl=(0.52)
#p_repl_2=(0.52 0.49)
p_repl=(1.05)
p_repl_2=(1.05 0.95)
p_mut=(0.0001)
L=200
N=100
T_FINAL=100
N_SIM=10
# max bin size
N_BIN=12
num_intro=(0 5 10 50 100)
#num_intro=(50 100 500 1000)
T_SWITCH_ORIG=(true false)

echo "Running simulation"

i=1

until [ $i -gt $N_SIM ]
do
	for pm in ${p_mut[@]}
	do
   		#echo $pm
   	 	for pr in ${p_repl[@]}
    		do
                for pr2 in ${p_repl_2[@]}
                do
                    for ni in ${num_intro[@]}
                    do
                        for sw in ${T_SWITCH_ORIG[@]}
                        do
                            #echo $pr
                            INTRO_REPL2=$p_repl_2
                            echo "python3 $WORK_DIR/evol_run.py $RESULT_DIR $pr $pm $L $N $T_FINAL $i $N_SIM $N_BIN $pr2 $ni $INTRO_REPL2 $sw"
                            $HOME/homebrew/bin/python3 $WORK_DIR/evol_run_modular.py $RESULT_DIR $pr $pm $L $N $T_FINAL $i $N_SIM $N_BIN $pr2 $ni $INTRO_REPL2 $sw
                        done
                    done
                done
    		done
    	echo i: $i
 		((i=i+1))
	done
done
