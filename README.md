# ZeroBind
This is the implementation of ZeroBind: A protein-specific zero-shot predictor with subgraph matching for drug-target interactions.
## Installation
ZeroBind is built on Python3, we recommend using a virtual conda environment as enviroment management for the installation of ZeroBind and its dependencies. The virtual environment can be created as follows:
```bash
conda create -n your_environment python==3.9
conda activate your_environment
```
Download the source code of ZeroBind from GitHub:
```bash
git clone https://github.com/myprecioushh/ZeroBind.git
```
Install ZeroBind dependencies as following:
```bash
conda install pytorch==2.0.0 torchvision==0.15.0 torchaudio==2.0.0 cu102 -c pytorch
conda install pyg==2.3.0 -c pyg
conda install lightning==2.0.1 -c conda-forge or pip install lightning==2.0.1
conda install -c conda-forge rdkit
pip install graphein
pip install fair-esm
```

## Crawl the data
```bash
nohup python -u ./crawler.py > crawler.log &
```

## Train
Multiple hyperparameters can be selected in meta.py. 
```bash
python metaentry.py  --batch_size=4  --num_workers=16 --num_inner_steps=5 --k_query=50
```
After model training starts, the progress bar will be automatically shown on your command line， and the trained model parameters will be saved in "checkpoints" dictory for every epoch.
## Prediction
```bash
python metaentry.py  --test --num_workers=16 --k_query=50 --checkpoint_path="your_model_path"
```
After predicting with your well trained model, the predicting output will be saved in a "txt" file and the prediction metrics will be shown on your command line.
## Online service and benchmark datasets
Online retrieval service and benchmark datasets are in [(http://www.csbio.sjtu.edu.cn/bioinf/ZeroBind/)](http://www.csbio.sjtu.edu.cn/bioinf/ZeroBind/index.html).

## License
This project is covered under the Apache 2.0 License.
