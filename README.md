We predict the electrospray ionization tandem mass spectrum of Chemically Derived Molecules (CDMs) using deep learning, named as DeepCDM. Rather than developing a new predicting tool from scratch, DeepCDM focuses on re-training existing algorisms using a small set of experimentally obtained CDM spectra via transfer learning. The performance is evaluated by the weighted cosine similarity (WCS) between experimental spectra and predicted spectra. 

This Code Package is Supplementary Information for the article "Deep learning enables high-quality MS/MS spectrum prediction for chemically derived molecules"

## Required Packages:
Python 3.6, RDKit, Tensorflow=1.13.2 is used for model training. 

## DeepCDM:
Ref codes from "Jennifer N. Wei, David Belanger, Ryan P. Adams, and D. Sculley. Rapid Prediction of Electron–Ionization Mass Spectrometry Using Neural Networks. ACS Central Science 2019 5 (4), 700-708 DOI: 10.1021/acscentsci.9b00085" are modified for model re-training. `molecule_estimator_transfer.py` and `molecule_predictors_transfer.py` are used for fine tuning. 

### Folders:
Folder “Ref Codes” contains the original codes of NEIMS. 

Folder “Tuning” contains codes modified from NEIMS and codes for transfer learning. 

“SI Codes” contains codes for library construction and spectra matching.

“Example” contains examples of molecular structural information for spectrum prediction.

“Datasets” contains example datasets for model pre-training, fine-tuning and validation. 

`merged_MoNA_ESI-MSMS_spectra.sdf` is used for model pre-training, `replicate_from_mona.sdf` is used for validation in pre-training.

`fine_tuning_example.sdf` is a demo training set used for model fine-tuning, while `replicate_from_dansylation_example.sdf` is a demo for validation. Larger training set is recommended for better performance. Moreover, the training set for fine-tuning can be changed for variable CDMs.


### 1.	Training, validation and test data split:
`make_train_test_split.py` is used to randomly split datasets for model re-training and fine-tuning.

#### 1)	Data split for model training:
```
python make_train_test_split.py \
--main_sdf_name=test_data/merged_MoNA_ESI-MSMS_spectra.sdf \
--replicates_sdf_name=test_data/replicate_from_mona.sdf \
--output_master_dir=tmp/pretrain/spectra_tf_records
```

#### 2)	Data split for fine-tuning:
```
python make_train_test_split.py \
--main_sdf_name=test_data/fine_tuning_example.sdf \
--replicates_sdf_name=test_data/replicate_from_dansylation_example.sdf \
--output_master_dir=tmp/finetuning/spectra_tf_records
```

### 2.	Model pre-training:
The MoNA dataset is used to re-train the architecture of NEIMS by `molecule_estimator.py` to establish ESI-MLP.
```
python molecule_estimator.py \
--dataset_config_file=tmp/pretrain/spectra_tf_records/query_replicates_val_predicted_replicates_val.json \
--train_steps=1000 \
--model_dir=models/pre_train \
--alsologtostderr
```

### 3.	Model fine tuning:
ESI-MLP is fine-tuned by `molecule_estimator_transfer.py` using the dansylated training set to make it specialized for dansylated molecule (Dns-MS). 
```
python molecule_estimator_transfer.py \
--dataset_config_file=tmp/finetuning/spectra_tf_records/query_replicates_val_predicted_replicates_val.json \
--train_steps=1000 \
--model_dir=models/fine_tuning \
-- warm_start_dir=models/pre_train \
--alsologtostderr
```

### 4.	Spectra prediction:
The CDM-specific model, such as Dns-MS, predicts spectra for CDMs using `make_spectra_prediction.py`.
```
python make_spectra_prediction.py \
--input_file=examples/emample_molecules.sdf \
--output_file=examples/emample_molecules_pred.sdf \
--weights_dir=models/finetuning
```

## Library Construction:
`1_dansylation_filtering.py` is used to extract molecules contains 14 chemical elements and whose exact mass are ≤ 1,000 Da.

`2_SMART_reaction.py` is used for virtual dansylation reaction.  

## Spectra Matching:
`3_spectra_matching.py` is used for DnsBank spectra matching.

`4_sort_match_res.py` is used to sort molecule annotations according to spectra similarity. 
