
#!/bin/bash

# type: Bulk or scStyleBulk


######################
# Simple Classifier
# Best configurations with gene expression
######################

# GSE117872HN120_Pri
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type Bulk --lr 0.0001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# GSE117872HN120_Met
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type Bulk --lr 0.0001 --hdim 256 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# GSE117872HN137_Pri
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type Bulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# GSE117872HN137_Met
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type Bulk --lr 0.0001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# GSE110894
python Simple_classifier.py --data_name "GSE110894" --type Bulk --lr 0.0001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE110894" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# GSE112274
python Simple_classifier.py --data_name "GSE112274" --type Bulk --lr 0.0001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE112274" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# GSE140440
python Simple_classifier.py --data_name "GSE140440" --type Bulk --lr 0.0001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE140440" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# GSE149383
python Simple_classifier.py --data_name "GSE149383" --type Bulk --lr 0.00001 --hdim 32 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "GSE149383" --type scStyleBulk --lr 0.00001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# Etoposide
python Simple_classifier.py --data_name "Etoposide" --type Bulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "Etoposide" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# PLX4720_A375
python Simple_classifier.py --data_name "PLX4720_A375" --type Bulk --lr 0.00001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "PLX4720_A375" --type scStyleBulk --lr 0.00001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# PLX4720_451Lu
python Simple_classifier.py --data_name "PLX4720_451Lu" --type Bulk --lr 0.0001 --hdim 32 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "PLX4720_451Lu" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# Afatinib
python Simple_classifier.py --data_name "Afatinib" --type Bulk --lr 0.0001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "Afatinib" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# AR-42
python Simple_classifier.py --data_name "AR-42" --type Bulk --lr 0.00001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "AR-42" --type scStyleBulk --lr 0.00001 --hdim 256 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# Cetuximab
python Simple_classifier.py --data_name "Cetuximab" --type Bulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "Cetuximab" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# Gefitinib
python Simple_classifier.py --data_name "Gefitinib" --type Bulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "Gefitinib" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# NVP-TAE684
python Simple_classifier.py --data_name "NVP-TAE684" --type Bulk --lr 0.0001 --hdim 128 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "NVP-TAE684" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# Sorafenib
python Simple_classifier.py --data_name "Sorafenib" --type Bulk --lr 0.00001 --hdim 64 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "Sorafenib" --type scStyleBulk --lr 0.00001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr

# Vorinostat
python Simple_classifier.py --data_name "Vorinostat" --type Bulk --lr 0.00001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr
python Simple_classifier.py --data_name "Vorinostat" --type scStyleBulk --lr 0.00001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model gexpr



######################
# Simple Classifier w/scRobust *
# Best configurations with scRobust foundation model
######################


# GSE117872HN120_Pri
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type Bulk --lr 0.00001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# GSE117872HN120_Met
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type Bulk --lr 0.0001 --hdim 128 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# GSE117872HN137_Pri
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type Bulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# GSE117872HN137_Met
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type Bulk --lr 0.00001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# GSE110894
python Simple_classifier.py --data_name "GSE110894" --type Bulk --lr 0.00001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE110894" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# GSE112274
python Simple_classifier.py --data_name "GSE112274" --type Bulk --lr 0.0001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE112274" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# GSE140440
python Simple_classifier.py --data_name "GSE140440" --type Bulk --lr 0.0001 --hdim 128 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE140440" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# GSE149383
python Simple_classifier.py --data_name "GSE149383" --type Bulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "GSE149383" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# Etoposide
python Simple_classifier.py --data_name "Etoposide" --type Bulk --lr 0.00001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "Etoposide" --type scStyleBulk --lr 0.00001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# PLX4720_A375
python Simple_classifier.py --data_name "PLX4720_A375" --type Bulk --lr 0.0001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "PLX4720_A375" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# PLX4720_451Lu
python Simple_classifier.py --data_name "PLX4720_451Lu" --type Bulk --lr 0.0001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "PLX4720_451Lu" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# Afatinib
python Simple_classifier.py --data_name "Afatinib" --type Bulk --lr 0.00001 --hdim 128 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "Afatinib" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# AR-42
python Simple_classifier.py --data_name "AR-42" --type Bulk --lr 0.00001 --hdim 256 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "AR-42" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# Cetuximab
python Simple_classifier.py --data_name "Cetuximab" --type Bulk --lr 0.0001 --hdim 256 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "Cetuximab" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# Gefitinib
python Simple_classifier.py --data_name "Gefitinib" --type Bulk --lr 0.0001 --hdim 256 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "Gefitinib" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# NVP-TAE684
python Simple_classifier.py --data_name "NVP-TAE684" --type Bulk --lr 0.00001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "NVP-TAE684" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# Sorafenib
python Simple_classifier.py --data_name "Sorafenib" --type Bulk --lr 0.00001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "Sorafenib" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust

# Vorinostat
python Simple_classifier.py --data_name "Vorinostat" --type Bulk --lr 0.00001 --hdim 256 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust
python Simple_classifier.py --data_name "Vorinostat" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scRobust



######################
# Simple Classifier w/scGPT *
# Best configurations with scGPT foundation model
######################

# GSE117872HN120_Pri
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type Bulk --lr 0.0001 --hdim 32 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# GSE117872HN120_Met
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type Bulk --lr 0.00001 --hdim 32 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type scStyleBulk --lr 0.00001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# GSE117872HN137_Pri
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type Bulk --lr 0.00001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# GSE117872HN137_Met
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type Bulk --lr 0.00001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# GSE110894
python Simple_classifier.py --data_name "GSE110894" --type Bulk --lr 0.00001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE110894" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# GSE112274
python Simple_classifier.py --data_name "GSE112274" --type Bulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE112274" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# GSE140440
python Simple_classifier.py --data_name "GSE140440" --type Bulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE140440" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# GSE149383
python Simple_classifier.py --data_name "GSE149383" --type Bulk --lr 0.00001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "GSE149383" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# Etoposide
python Simple_classifier.py --data_name "Etoposide" --type Bulk --lr 0.0001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "Etoposide" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# PLX4720_A375
python Simple_classifier.py --data_name "PLX4720_A375" --type Bulk --lr 0.0001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "PLX4720_A375" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# PLX4720_451Lu
python Simple_classifier.py --data_name "PLX4720_451Lu" --type Bulk --lr 0.0001 --hdim 64 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "PLX4720_451Lu" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# Afatinib
python Simple_classifier.py --data_name "Afatinib" --type Bulk --lr 0.0001 --hdim 64 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "Afatinib" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# AR-42
python Simple_classifier.py --data_name "AR-42" --type Bulk --lr 0.00001 --hdim 32 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "AR-42" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# Cetuximab
python Simple_classifier.py --data_name "Cetuximab" --type Bulk --lr 0.0001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "Cetuximab" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# Gefitinib
python Simple_classifier.py --data_name "Gefitinib" --type Bulk --lr 0.0001 --hdim 64 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "Gefitinib" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# NVP-TAE684
python Simple_classifier.py --data_name "NVP-TAE684" --type Bulk --lr 0.00001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "NVP-TAE684" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# Sorafenib
python Simple_classifier.py --data_name "Sorafenib" --type Bulk --lr 0.00001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "Sorafenib" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT

# Vorinostat
python Simple_classifier.py --data_name "Vorinostat" --type Bulk --lr 0.00001 --hdim 64 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT
python Simple_classifier.py --data_name "Vorinostat" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scGPT



######################
# Simple Classifier w/scFoundation *
# Best configurations with scFoundation foundation model
######################

# GSE117872HN120_Pri
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type Bulk --lr 0.00001 --hdim 32 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE117872HN120_Pri" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# GSE117872HN120_Met
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type Bulk --lr 0.00001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE117872HN120_Met" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# GSE117872HN137_Pri
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type Bulk --lr 0.00001 --hdim 32 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE117872HN137_Pri" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# GSE117872HN137_Met
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type Bulk --lr 0.00001 --hdim 64 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE117872HN137_Met" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# GSE110894
python Simple_classifier.py --data_name "GSE110894" --type Bulk --lr 0.00001 --hdim 32 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE110894" --type scStyleBulk --lr 0.0001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# GSE112274
python Simple_classifier.py --data_name "GSE112274" --type Bulk --lr 0.00001 --hdim 32 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE112274" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# GSE140440
python Simple_classifier.py --data_name "GSE140440" --type Bulk --lr 0.0001 --hdim 64 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE140440" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# GSE149383
python Simple_classifier.py --data_name "GSE149383" --type Bulk --lr 0.0001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "GSE149383" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# Etoposide
python Simple_classifier.py --data_name "Etoposide" --type Bulk --lr 0.00001 --hdim 32 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "Etoposide" --type scStyleBulk --lr 0.0001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# PLX4720_A375
python Simple_classifier.py --data_name "PLX4720_A375" --type Bulk --lr 0.00001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "PLX4720_A375" --type scStyleBulk --lr 0.00001 --hdim 32 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# PLX4720_451Lu
python Simple_classifier.py --data_name "PLX4720_451Lu" --type Bulk --lr 0.0001 --hdim 32 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "PLX4720_451Lu" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# Afatinib
python Simple_classifier.py --data_name "Afatinib" --type Bulk --lr 0.00001 --hdim 32 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "Afatinib" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# AR-42
python Simple_classifier.py --data_name "AR-42" --type Bulk --lr 0.0001 --hdim 32 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "AR-42" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# Cetuximab
python Simple_classifier.py --data_name "Cetuximab" --type Bulk --lr 0.00001 --hdim 256 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "Cetuximab" --type scStyleBulk --lr 0.00001 --hdim 256 --batch 32 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# Gefitinib
python Simple_classifier.py --data_name "Gefitinib" --type Bulk --lr 0.00001 --hdim 32 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "Gefitinib" --type scStyleBulk --lr 0.0001 --hdim 256 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# NVP-TAE684
python Simple_classifier.py --data_name "NVP-TAE684" --type Bulk --lr 0.0001 --hdim 32 --batch 64 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "NVP-TAE684" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# Sorafenib
python Simple_classifier.py --data_name "Sorafenib" --type Bulk --lr 0.00001 --hdim 256 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "Sorafenib" --type scStyleBulk --lr 0.0001 --hdim 32 --batch 16 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation

# Vorinostat
python Simple_classifier.py --data_name "Vorinostat" --type Bulk --lr 0.00001 --hdim 256 --batch 128 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation
python Simple_classifier.py --data_name "Vorinostat" --type scStyleBulk --lr 0.00001 --hdim 128 --batch 8 --wait 5 --dropout 0.5 --bn False --foundation_model scFoundation