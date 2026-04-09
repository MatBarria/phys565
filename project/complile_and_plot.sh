make clean
make -j4
python3 script/run_tuples.py 
cd python 
python3 draw_histograms_noratio.py 
cd ..
