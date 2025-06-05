#!/bin/bash

../xcftools_static view --i reference.bcf -o reference_xcf.bcf -O sh -r 20 -T 8 -m 0.03125

../impute5_static --h reference_xcf.bcf --g target.bcf --o imputed.bcf --r 20:1200000-3800000 --m chr20.b37.gmap.gz --buffer-region 20:1000000-4000000

