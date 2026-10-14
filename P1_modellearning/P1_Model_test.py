import os
import argparse
import itertools
import pandas as pd
from SurrealGAN.Surreal_GAN_representation_learning import repetitive_representation_learning

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task_id", type=int, required=True, help="0-based task index, range: 0-1349")
    args = parser.parse_args()

    # =========================
    # =========================
    base_output_dir = "/cpfs01/projects-HDD/cfff-afe2df89e32e_HDD/public/gyj_data/yanf/STAAR_QSM/Part2_MODEL_learning/Multidisease"
    os.makedirs(base_output_dir, exist_ok=True)
    
    data_input_dir = "/cpfs01/projects-HDD/cfff-afe2df89e32e_HDD/public/gyj_data/yanf/STAAR_QSM/Part2_MODEL_learning"
    train_data_path = os.path.join(data_input_dir, "dt_Model_learning_1014.csv")
    train_data = pd.read_csv(train_data_path)

    # =========================
    # =========================
    batchsize = 300
    lr = 0.0008
    repetition = 50

    tar_size = len(train_data[train_data["diagnosis"] == 1])
    if tar_size == 0:
        raise ValueError("diagnosis == 1 的 TAR 组样本数为 0，请检查输入数据。")
    final_saving_epoch = int(1500000 * (batchsize / tar_size))

    # =========================
    # =========================
    npattern_list = [2, 3, 4, 5, 6]
    lam_list = [0.1, 0.5, 0.8]
    gamma_list = [0.1, 0.5, 0.8]

    param_grid = list(itertools.product(npattern_list, lam_list, gamma_list))
    total_tasks = len(param_grid) * repetition

    if args.task_id < 0 or args.task_id >= total_tasks:
        raise ValueError(f"task_id 超出范围，应在 0 到 {total_tasks-1} 之间。")

    grid_idx = args.task_id // repetition
    rep_idx = (args.task_id % repetition) # 0-based repetition index

    npattern, lam, gamma = param_grid[grid_idx]

    current_out_dir = os.path.join(
        base_output_dir,
        f"Model_M{npattern}_Lam{lam}_Gam{gamma}"
    )
    os.makedirs(current_out_dir, exist_ok=True)

    print(f"当前任务 task_id = {args.task_id}")
    print(f"grid_idx = {grid_idx}, repetition = {rep_idx}")
    print(f"模式数(M): {npattern} | Lambda: {lam} | Gamma: {gamma}")
    print(f"自适应 Epochs: {final_saving_epoch}")
    print(f"输出目录: {current_out_dir}")

    repetitive_representation_learning(
        train_data,
        npattern,
        repetition=repetition,
        start_repetition=rep_idx,
        stop_repetition=rep_idx + 1,
        final_saving_epoch=final_saving_epoch,
        output_dir=current_out_dir,
        fraction=1,
        lr=lr,
        batchsize=batchsize,
        verbose=False,
        lipschitz_k=0.5,
        covariate=None,
        lam=lam,
        gamma=gamma
    )

    print("Single repetition complated")

if __name__ == "__main__":
    main()