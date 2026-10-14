The model is trained using the method implemented in [CBICA/SurrealGAN](https://github.com/CBICA/SurrealGAN).
`model_agreements.csv` contains the results generated using the `Model_M6_Lam0.1_Gam0.8` configuration.

## Training sanity checks

SurrealGAN's patient batch iterator (`SurrealGAN/data_loading.py`, `PTIterator`) runs `int(n_patients / batchsize) - 1` training steps per epoch. With 411 patients and `batchsize=300` that is **zero steps**: the networks keep their random starting weights, every saved epoch is identical, and `Rindices_corr` in `model_agreements.csv` is the same at every epoch. Depending on the installed version, such a run either crashes at epoch 100 or saves untrained models without any message.

- `P1_Model_test.py` now stops before training if the batch size leaves no steps, and takes `--batchsize` (default 300, unchanged).
- Before submitting: `python check_model_training.py --data dt_Model_learning_1014.csv --batchsize 100` prints the steps per epoch and the largest usable batch size (205 for one step per epoch, 102 for at least three with 411 patients).
- After training: `python check_model_training.py /path/to/Multidisease` flags every configuration folder whose agreement or saved weights never change (`FROZEN`).
- `pytest test_check_model_training.py` tests the checks (needs SurrealGAN and PyTorch for the step and weight checks).
- The results in `model_agreements.csv` above (`Model_M6_Lam0.1_Gam0.8`) show this problem: all 33 checkpoints, epochs 4000 to 68000, are identical. Do not use R-indices from such runs.
