DATASETS="Superchris,Stella,Barat,Buchanan,Mitt"
TIMEWINDOW=0.250
INPLEN=10
RESULTSDIR="results/"
RESULTSDIR_CSV="results_csv/"

echo "Processing pipine for time window of $TIMEWINDOW seconds"

### Do dimensionality reduction
python -u scripts/dim_reduction.py --timewindow $TIMEWINDOW --inpLen $INPLEN --datasets $DATASETS | tee $RESULTSDIR/dim_reduction_results.txt

### Do dimensionality reduction
python scripts/classification.py --timewindow $TIMEWINDOW --inpLen $INPLEN --num_shuffles 100 --datasets $DATASETS

### Plot encodings and separation
python scripts/plot_results.py --timewindow $TIMEWINDOW --inpLen $INPLEN --results_dir $RESULTSDIR/ --datasets $DATASETS

### Plot encodings and separation
python scripts/convert_result_files.py --timewindow $TIMEWINDOW --inpLen $INPLEN --datasets $DATASETS --out_dir $RESULTSDIR_CSV
