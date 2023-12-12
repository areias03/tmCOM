#!/bin/zsh


echo Starting all samples

echo Writing results to ../data/results/

conda activate tmcom

python sample_runner.py 1 &

python sample_runner.py 1 -li &

python sample_runner.py 1 --ec &

python sample_runner.py 1 -li --ec &

python sample_runner.py 2 &

python sample_runner.py 2 -li &

python sample_runner.py 2 --ec &

python sample_runner.py 2 -li --ec &
