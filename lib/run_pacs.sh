#!/bin/bash


if [ -z "$1" ]
then
   echo :Usage ./$( basename $0 ) FIRST_CYCLE LAST_CYCLE NUM_OF_CANDICATES
   exit
fi

WDIR=$(pwd)
RUNID="$(cat /dev/urandom | base64 | tr -dc "A-Z" | fold -w 5 | head -n 1)"
FIRST_CYCLE=$1
LAST_CYCLE=$2
NUM_OF_CANDICATES=$3

source $WDIR/lib/preference.sh

for curr_cycle in $(seq $FIRST_CYCLE $LAST_CYCLE) ; do
  
  #If the first cycle is zero, initialization is performed
  #When you already run initialization, the FIRST_CYCLE should be 1
  if [ $curr_cycle -eq 0 ]
  then
    $WDIR/lib/run_initial.sh initial $NUM_OF_CANDICATES $RUNID &
    wait
    continue
  fi

  #Make directory for the current cycle
  mkdir $WDIR/cyc${curr_cycle}
  cd $WDIR/cyc${curr_cycle}

  #Extract the ranked structures as .gro file. 
  while IFS=, read pre_rank pre_cycle pre_candi pre_frame pre_score ; do
    mkdir $WDIR/cyc${curr_cycle}/candi${pre_rank}
    cd $WDIR/cyc${curr_cycle}/candi${pre_rank}

    #前回の .xtc ファイルから .gro を抽出する
    $GMX trjconv -f $WDIR/cyc$((curr_cycle-1))/candi${pre_candi}/sample.xtc -s $WDIR/input/em.gro -o ./cyc${curr_cycle}_candi${pre_rank}.gro -dump ${pre_frame} << eof1
0
eof1
    
  done < $WDIR/cyc$((curr_cycle-1))/rank.csv
  
  #入力ファイルの作成
  for candi in $(seq 1 $NUM_OF_CANDICATES) ; do
    cd $WDIR/cyc${curr_cycle}/candi${candi}
    #初期の sample.mdpを再利用すれば良い
  done


  #MDの実行 と スコアの計算
  cd $WDIR/cyc${curr_cycle}

  for i in $(seq 1 $NUM_OF_CANDICATES) ; do
    cd $WDIR/cyc${curr_cycle}/candi${i}
    name=cyc${curr_cycle}_candi${i}

    $GMX grompp -maxwarn -1 -f $WDIR/input/sample.mdp -c ./${name}.gro -p $WDIR/input/initial.top -o ./sample.tpr
    $GMX_MPI mdrun -deffnm sample -ntomp ${OMP}
  done
  
  echo "All candicates are ready."

  #ランクを付ける →　条件を満たせば終了
  for i in $(seq 1 $NUM_OF_CANDICATES) ; do
    cd $WDIR/cyc${curr_cycle}/candi${i}
    python3 $WDIR/lib/$SCORING_SCRIPT
  done

  #ランクを付ける
  cd $WDIR/cyc${curr_cycle}
  python3 $WDIR/lib/merge_score.py ${curr_cycle} $NUM_OF_CANDICATES 
  

done

echo "PaCS-MD has done."
