#!/bin/sh
# Note: Plotly, for inscrutable reasons, creates an extra blank page
# in the PDF files of most of these plots.
# I delete those extra pages manually in Preview.
./get_stats.py --sankey -H --stacked --bar --font Optima,Arial
# To get a png, run
# magick -density 600 layer_class_sankey.pdf layer_class_sankey.png
# To get eps, run
# pdftops layer_class_sankey.pdf layer_class_sankey.ps
# ps2eps layer_class_sankey.ps
