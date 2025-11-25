#!/bin/bash

# Generate scStyleBulk data
python ZINB_Transform.py \
    --bulk_data 'data/gexpr/Bulk/Bulk_Etoposide_HVG.csv' \
    --sc_data 'data/gexpr/singlecell/Etoposide_HVG.csv' \
    --save_path './' \
    --drug 'Etoposide'

