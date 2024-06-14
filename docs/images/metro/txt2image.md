# Install desktop app

Got to `https://github.com/jgraph/drawio-desktop/releases/` and download the latest version for your OS.

To install it on wsl

```bash
sudo apt install /mnt/c/Users/llenezet/Dowlnoads/drawio-amd64-21.6.8.deb
```

To use drawio

```bash
drawio --version
drawio docs/images/metro/MetroMap.xml --export --format png --page-index 2 --layers 1 --output docs/images/metro/MetroMap.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --page-index 3 --layers 1 --output docs/images/metro/Simulate.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --page-index 4 --layers 0 --output docs/images/metro/Validate.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --page-index 5 --layers 0 --output docs/images/metro/PanelPrep.png --scale 3
drawio docs/images/metro/MetroMap.xml --export --format png --page-index 6 --layers 1 --output docs/images/metro/Impute.png --scale 3
```
