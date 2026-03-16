HPhi Post-Processing Toolkit: Physical Observable Calculator

简介:
本工具用于处理 HPhi（专家模式）生成的输出文件，计算零温、有限温及全对角化下的多种物理量。它能够自动识别并处理由 HPhiBuild （Run.py）产生的以 KG_thetaxxxPI_hxxx 命名的输出文件，并将计算结果分类存储到对应的文件夹中，便于后续分析与可视化。

注意:
1. main.cpp 中的 “SAMPLE_NUM” 的参数控制是基态还是激发态，0：基态；1：第一激发态；2：第二激发态...
2. 各个计算结果的文件必须命名为 "KG_xxxPI_xxxh"，必须与 Common.cpp 中 "GetAllFilesWithPrefix" 函数中的设置相对应。

主要功能
1. 零温计算 
基于 HPhi 的 LOBPCG 结果，可计算以下物理量：

保真度（Fidelity）

磁化强度（Magnetization）

单体自旋关联 <S^α_i>

二体自旋关联 <S^α_i S^β_j>，包括矢量手性（vector chirality）

标量手性关联（scalar chirality）：完整关联以及去除了低阶关联贡献的部分

四体关联：例如手性-手性关联 <χ_i χ_j>

Wilson Loop 算符：Wx, Wy, Wz

六体算符：Wp（plaquette 算符）及其关联 <Wp Wp>

纠缠熵（Entanglement Entropy）

拓扑纠缠熵（Topological Entanglement Entropy）

注意：零温计算提供了两种实现方式：

基于波函数：直接读取 HPhi 输出的本征态文件（仅支持自旋 1/2 系统）。

基于格林函数：利用 HPhi 生成的多体格林函数文件（支持任意自旋大小）。

2. 有限温计算
基于 HPhi 的 cTPQ（典型热纯态）方法生成的数据，计算：

比热（Specific Heat）

热熵（Thermal Entropy）

有限温单体自旋关联 <S^α_i>

有限温二体自旋关联 <S^α_i S^β_j>

有限温六体关联（如 Wp 的期望值）

有限温磁化率（Magnetic Susceptibility）

对于有限温结果，除了对多次 cTPQ 样本进行简单平均外，还提供了 Bootstrap 重采样 处理方式（HPhi 官方推荐），以更准确地估计统计误差，消除随机初态带来的涨落。

3. 全对角化计算 (Full Diagonalization)
基于 HPhi 的全对角化输出，计算与有限温相同的物理量，包括比热、热熵、单体/二体/六体关联、磁化率。全对角化结果可用作精确基准，用于验证 cTPQ 方法的收敛性。

文件组织
工具自动按照计算类型将输入文件和输出结果存入以下文件夹：

输入文件夹：存放 HPhi 生成的原始输出文件。文件名格式应为 KG_xxxPI_xxh（由 HPhiBuild 中的 Run.py 脚本定义）。

零温输入：inputCG/

有限温输入：inputTPQ/

全对角化输入：inputFull/

输出文件夹：程序运行后，计算结果将以文本文件形式保存在当前工作目录（或指定输出路径），文件名标识物理量、温度及系统参数。

