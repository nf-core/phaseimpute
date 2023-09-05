```bash
conda create -n txt2image
conda activate txt2image
conda install -c "conda-forge/label/cf202003" pyvips
python ./txt2img.py MetroMap.xml test.png
```