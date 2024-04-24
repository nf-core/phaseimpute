# Install desktop app

Got to `https://github.com/jgraph/drawio-desktop/releases/` and download the latest version for your OS.

To install it on wsl

```bash
sudo apt install /mnt/c/Users/llenezet/Dowlnoads/drawio-amd64-21.6.8.deb
```

To use drawio

```bash
drawio --version
drawio docs/images/metro/MetroMap.xml --export --format png --page-index 0 --output docs/images/metro/MetroMap.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --layers 0 --page-index 1 --output docs/images/metro/PostProcessing.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --layers 1 --page-index 1 --output docs/images/metro/Concordance2.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --layers 2 --page-index 1 --output docs/images/metro/Simulate.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --layers 3 --page-index 1 --output docs/images/metro/Phase.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --layers 4 --page-index 1 --output docs/images/metro/PreProcessing.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --layers 5 --page-index 1 --output docs/images/metro/Concordance.png --scale 3
```
