#!/bin/sh
# Note: Plotly, for inscrutable reasons, creates an extra blank page
# in the PDF files of most of these plots.
# I delete those extra pages manually in Preview.
./get_stats.py --sankey -H --stacked --bar --font Optima,Arial
