### Instructions
Fetch the data from Zenodo:
```
for i in {01..22}; do wget "https://zenodo.org/records/6899613/files/simulations_${i}.zip?download=1" -O "simulations_${i}.zip"; done
```
Unpack the data into a subdirectory `simulations`:
```
mkdir -p ./simulations
for i in {1..20}; do batch="simulations_$(printf "%02d" $i)"; unzip "$batch.zip" && mv ./"$batch"/* ./simulations && rmdir ./"$batch"; done
```
Get [stableCoalescence_cladeAnalysis.py](https://github.com/nizzaneela/multi-introduction/blob/Resample-of-cladeAnalysis_stableCoalescence.py/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) and [cladeAnalysis.ipynb](https://github.com/nizzaneela/multi-introduction/blob/Resample-of-cladeAnalysis_stableCoalescence.py/notebooks/cladeAnalysis.ipynb) from this branch:
```
curl -L -o cladeAnalysis.ipynb "https://raw.githubusercontent.com/nizzaneela/multi-introduction/Resample-of-cladeAnalysis_stableCoalescence.py/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py"
curl -L -o cladeAnalysis.ipynb "https://raw.githubusercontent.com/nizzaneela/multi-introduction/Resample-of-cladeAnalysis_stableCoalescence.py/notebooks/cladeAnalysis.ipynb"
```
Launch and run the notebook.

### Supplementary code for:

J. E. Pekar, A. Magee, E. Parker, N. Moshiri, K. Izhikevich, J. L. Havens, K. Gangavarapu, L. M. Malpica Serrano, A. Crits-Christoph, N. L. Matteson, M. Zeller, J. I. Levy, J. C. Wang, S. Hughes, J. Lee, H. Park, M.-S. Park, K. Ching Zi Yan, R. T. Pin Lin, M. N. Mat Isa, Y. M. Noor, T. I. Vasylyeva, R. F. Garry, E. C. Holmes, A. Rambaut, M. A. Suchard, K. G. Andersen, M. Worobey, J. O. Wertheim, "The molecular epidemiology of multiple zoonoses of SARS-CoV-2".
