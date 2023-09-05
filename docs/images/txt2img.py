#!/usr/bin/env python3

import sys
import pyvips
from xml.sax.saxutils import escape

# load first arg as a string
txt = open(sys.argv[1], "r").read()

# pyvips allows pango markup in strings -- you can write stuff like
# text("hello <i>sailor!</i>")
# so we need to escape < > & in the text file
txt = escape(txt)

img = pyvips.Image.text(txt)

# save to second arg
img.write_to_file(sys.argv[2])