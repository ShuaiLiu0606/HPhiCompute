#include <common.h>
#include <TPQ.h>

TPQ::TPQ(const string &ParaFilePath, const string file_theta, const string file_h, int sample_num) : Nsites(0), Sdim(0), Hdim(0), NumAve(0), Lanczos_Max(0), ExpecInterval(0), LargeValue(0.0), Sample_num(sample_num), File_theta(file_theta), File_h(file_h)
{
    const string ParaFile = ParaFilePath + "/../modpara.def";

    // NumAve: TPQ独立运行的次数
    string AA = "NumAve";
    GetPara(ParaFile, AA, NumAve);

    // Lanczos_max: TPQ每次运行的最大迭代次数
    string BB = "Lanczos_max";
    GetPara(ParaFile, BB, Lanczos_Max);

    // LargeValue: 虚时演化的最大值
    string CC = "LargeValue";
    GetPara(ParaFile, CC, LargeValue);

    // ExpecInterval: 输出关联函数的步长
    string FF = "ExpecInterval";
    GetPara(ParaFile, FF, ExpecInterval);

    string DD = "Nsite";           // 关键字
    GetPara(ParaFile, DD, Nsites); // 获取参数

    const string LocspnFile = ParaFilePath + "/../locspn.def";

    string EE = to_string(Nsites - 1); // 关键字
    int S2 = 0.0;
    GetPara(LocspnFile, EE, S2); // 获取参数

    Spin = S2 / 2.0;
    Sdim = S2 + 1;
    Hdim = pow(Sdim, Nsites);

    cout << "Nsite: " << Nsites << endl;
    cout << "Spin: " << Spin << endl;
    cout << "Sdim: " << Sdim << endl;
    cout << "Hdim: " << Hdim << endl;
    cout << "NumAve: " << NumAve << endl;
    cout << "Lanczos_Max: " << Lanczos_Max << endl;
    cout << "ExpecInterval: " << ExpecInterval << endl;
    cout << "LargeValue: " << LargeValue << endl;
    cout << endl;
}

void TPQ::mainPhys(const string &TPQ_FilePath)
{
    my_mat SxCorr, SyCorr, SzCorr;
    my_mat SxSxCorr, SySyCorr, SzSzCorr;

    my_mat Wp_mat;

    SxCorr.resize(NumAve, my_vec(Lanczos_Max / ExpecInterval));
    SyCorr.resize(NumAve, my_vec(Lanczos_Max / ExpecInterval));
    SzCorr.resize(NumAve, my_vec(Lanczos_Max / ExpecInterval));
    SxSxCorr.resize(NumAve, my_vec(Lanczos_Max / ExpecInterval));
    SySyCorr.resize(NumAve, my_vec(Lanczos_Max / ExpecInterval));
    SzSzCorr.resize(NumAve, my_vec(Lanczos_Max / ExpecInterval));
    Wp_mat.resize(NumAve, my_vec(Lanczos_Max / ExpecInterval));

    S_Corr_TPQ(TPQ_FilePath, SxCorr, SyCorr, SzCorr);
    SS_Corr_TPQ(TPQ_FilePath, SxSxCorr, SySyCorr, SzSzCorr);
    WpMoments_TPQ(TPQ_FilePath, Wp_mat);

    my_vec ave_Sx = cla_mean(SxCorr, 0);
    my_vec err_Sx = cal_standard_deviation(SxCorr, 0, 1);
    my_vec ave_Sy = cla_mean(SyCorr, 0);
    my_vec err_Sy = cal_standard_deviation(SyCorr, 0, 1);
    my_vec ave_Sz = cla_mean(SzCorr, 0);
    my_vec err_Sz = cal_standard_deviation(SzCorr, 0, 1);

    //

    my_vec ave_SxSx = cla_mean(SxSxCorr, 0);
    my_vec err_SxSx = cal_standard_deviation(SxSxCorr, 0, 1);
    my_vec ave_SySy = cla_mean(SySyCorr, 0);
    my_vec err_SySy = cal_standard_deviation(SySyCorr, 0, 1);
    my_vec ave_SzSz = cla_mean(SzSzCorr, 0);
    my_vec err_SzSz = cal_standard_deviation(SzSzCorr, 0, 1);

    //

    my_vec ave_Wp = cla_mean(Wp_mat, 0);
    my_vec err_Wp = cal_standard_deviation(Wp_mat, 0, 1);

    // 写入TPQ_SimpleAverage.dat文件

    const string SimpleFile = "cTPQ_Physics_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream my_file1(SimpleFile);

    my_file1 << "#  T    Sx    Sx_err    Sy    Sy_err    Sz    Sz_err    SxSx    SxSx_err    SySy    SySy_err    SzSz    SzSz_err    Chi_x    Chi_y    Chi_z    Wp    Wp_err" << endl;

    for (int k = 1; k < Lanczos_Max / ExpecInterval; k++)
    {
        double temp = 1 / (k * ExpecInterval * (1 / LargeValue));

        // 正确的磁化率计算
        double chiX = (ave_SxSx[k] - ave_Sx[k] * ave_Sx[k]) / temp / Nsites;
        double chiY = (ave_SySy[k] - ave_Sy[k] * ave_Sy[k]) / temp / Nsites;
        double chiZ = (ave_SzSz[k] - ave_Sz[k] * ave_Sz[k]) / temp / Nsites;

        my_file1 << fixed << setprecision(16) << left << setw(22) << temp // 1
                 << setw(22) << ave_Sx[k]                                 // 2
                 << setw(22) << err_Sx[k]                                 // 3
                 << setw(22) << ave_Sy[k]                                 // 4
                 << setw(22) << err_Sy[k]                                 // 5
                 << setw(22) << ave_Sz[k]                                 // 6
                 << setw(22) << err_Sz[k]                                 // 7
                 << setw(22) << ave_SxSx[k]                               // 8
                 << setw(22) << err_SxSx[k]                               // 9
                 << setw(22) << ave_SySy[k]                               // 10
                 << setw(22) << err_SySy[k]                               // 11
                 << setw(22) << ave_SzSz[k]                               // 12
                 << setw(22) << err_SzSz[k]                               // 13
                 << setw(22) << chiX                                      // 14
                 << setw(22) << chiY                                      // 15
                 << setw(22) << chiZ                                      // 16
                 << setw(22) << ave_Wp[k]
                 << setw(22) << err_Wp[k]
                 << endl;
    }
}

void TPQ::mainBasic(const string &TPQ_FilePath)
{
    // 读取TPQ的数据文件
    //  波函数的范数的平方 ，逆温度， 能量， 能量平方， Sz, Sz2, 磁化率; 矩阵的大小为(NumAve * Lanczos_Max)
    my_mat Norm, InvTemp, Ene, Ene2, Spc, Sz, Sz2, chi;
    Read_TPQ_file(TPQ_FilePath, Norm, InvTemp, Ene, Ene2, Spc, Sz, Sz2, chi);

    // 计算并写入Norm.dat文件
    my_vec log_Z(Lanczos_Max, 0.0);
    double tmp_log_Z = -log(Norm[0][0]);
    for (int k = 0; k < Lanczos_Max; k++)
    {
        tmp_log_Z += log(Norm[0][k]); // 第0个cTPQ态的平方范数求ln
        log_Z[k] = tmp_log_Z;         // log_Z[0] = 1;即：Norm[n][0] = 1; 即：<phi_rand|phi_rand> = 1.0
    }

    my_mat Ent(NumAve, my_vec(Lanczos_Max, 0.0));
    my_mat NewNorm = Norm;
    for (int i = 0; i < NumAve; i++)
    {
        double logZ = -log(NewNorm[i][0]);
        for (int j = 0; j < Lanczos_Max; j++)
        {
            logZ += log(NewNorm[i][j]);

            Ent[i][j] = (logZ + InvTemp[i][j] * Ene[i][j] + log(Hdim)) / log(Hdim); // Ent = (log(Z) + <H>/T) / log(Hdim)
        }
    }

    // 其中每一行代表一次TPQ方法计算的物理量，每一列代表一种物理值
    // 计算平均值和标准差
    my_vec ave_InvTemp = cla_mean(InvTemp, 0);

    my_vec ave_Ene = cla_mean(Ene, 0);
    my_vec err_Ene = cal_standard_deviation(Ene, 0, 1); // ddof: 如果为 true，使用无偏估计，即分母使用 N-1（样本大小减 1）
    my_vec ave_Spc = cla_mean(Spc, 0);
    my_vec err_Spc = cal_standard_deviation(Spc, 0, 1);
    my_vec ave_Sz = cla_mean(Sz, 0);
    my_vec err_Sz = cal_standard_deviation(Sz, 0, 1);
    my_vec ave_Sz2 = cla_mean(Sz2, 0);
    my_vec err_Sz2 = cal_standard_deviation(Sz2, 0, 1);
    my_vec ave_chi = cla_mean(chi, 0);
    my_vec err_chi = cal_standard_deviation(chi, 0, 1);
    my_vec ave_En = cla_mean(Ent, 0);
    my_vec err_En = cal_standard_deviation(Ent, 0, 1);

    // 写入TPQ_SimpleAverage.dat文件

    const string SimpleFile = "cTPQ_Simple_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream my_file1(SimpleFile);

    my_file1 << "#  T    E    Cv   Cv_err    S    S_err    Sz    Sz_err    Sz2    Sz2_err    Chi_z    Chi_err" << endl;
    for (int k = 1; k < Lanczos_Max; k++)
    {

        double temp = 1.0 / ave_InvTemp[k];
        my_file1 << fixed << setprecision(16) << left << setw(22) << temp // 1
                 << setw(22) << ave_Ene[k]                                // 2
                 << setw(22) << ave_Spc[k]                                // 3
                 << setw(22) << err_Spc[k]                                // 4
                 << setw(22) << ave_En[k]                                 // 5
                 << setw(22) << err_En[k]                                 // 6
                 << setw(22) << ave_Sz[k]                                 // 7
                 << setw(22) << err_Sz[k]                                 // 8
                 << setw(22) << ave_Sz2[k]                                // 9
                 << setw(22) << err_Sz2[k]                                // 10
                 << setw(22) << ave_chi[k]                                // 11
                 << setw(22) << err_chi[k]                                // 12
                 << endl;
    }

    // 计算规范化物理量
    CalcBasic(Norm, Ene, Ene2, Sz, Sz2, InvTemp, log_Z); // CalcBasic 不需要log_Z
}

void TPQ::CalcBasic(const my_mat &Norm, const my_mat &Ene, const my_mat &Ene2, const my_mat &Sz, const my_mat &Sz2, const my_mat &InvTemp, my_vec &log_Z)
{
    my_mat phys_Z, phys_Ene, phys_Ene2, phys_Sz, phys_Sz2, phys_InvTemp;

    phys_Z.resize(NumAve, my_vec(Lanczos_Max, 0.0));
    phys_Ene.resize(NumAve, my_vec(Lanczos_Max, 0.0));
    phys_Ene2.resize(NumAve, my_vec(Lanczos_Max, 0.0));
    phys_Sz.resize(NumAve, my_vec(Lanczos_Max, 0.0));
    phys_Sz2.resize(NumAve, my_vec(Lanczos_Max, 0.0));
    phys_InvTemp.resize(NumAve, my_vec(Lanczos_Max, 0.0));

    // 第一个样本的处理
    for (int k = 0; k < Lanczos_Max; k++)
    {
        // cnt_samp=0, Z is always 1
        // 第‘0’个样本,Z为1，是因为要将第‘0’个样本的平方范数作为 W_beat = <phi^0_beat|phi^0_beat>
        // 然后让所有样本的平方范数除以 W_beat，达到归一化的目的

        phys_Z[0][k] = 1.0;
        phys_Ene[0][k] = Ene[0][k];
        phys_Ene2[0][k] = Ene2[0][k];
        phys_Sz[0][k] = Sz[0][k];
        phys_Sz2[0][k] = Sz2[0][k];

        phys_InvTemp[0][k] = InvTemp[0][k];

        // 对于第一个文件
    }

    // 处理其他样本
    for (int cnt_samp = 1; cnt_samp < NumAve; cnt_samp++)
    {
        // k=0 for each sample, the 1st norm should be 1； <rand|rand> = 1.0

        // D0=N0=<phi_rand|phi_rand>=1.0
        double tot_Z = 1;
        double tot_Ene = Ene[cnt_samp][0];
        double tot_Ene2 = Ene2[cnt_samp][0];
        double tot_Sz = Sz[cnt_samp][0];
        double tot_Sz2 = Sz2[cnt_samp][0];

        phys_Z[cnt_samp][0] = tot_Z;
        phys_Ene[cnt_samp][0] = tot_Ene;
        phys_Ene2[cnt_samp][0] = tot_Ene2;
        phys_Sz[cnt_samp][0] = tot_Sz;
        phys_Sz2[cnt_samp][0] = tot_Sz2;

        phys_InvTemp[cnt_samp][0] = InvTemp[cnt_samp][0];

        // 处理其他步骤
        const string tpq_type = "cTPQ"; // 固定为cTPQ
        for (int k = 1; k < Lanczos_Max; k++)
        {
            double IPL_Z, IPL_Ene, IPL_Ene2, IPL_Sz, IPL_Sz2, IPL_InvTemp, IPL_Wp;

            if (tpq_type == "mTPQ")
            {
                // this part will be modified, now not used
                int ext_k;
                if (InvTemp[0][k] > InvTemp[cnt_samp][k])
                {
                    ext_k = k + 1;
                    if (InvTemp[0][k] > InvTemp[cnt_samp][k + 1])
                    {
                        ext_k = k + 2;
                    }
                }
                else if (InvTemp[0][k] < InvTemp[cnt_samp][k])
                {
                    ext_k = k - 1;
                    if (InvTemp[0][k] < InvTemp[cnt_samp][k - 1])
                    {
                        ext_k = k - 2;
                    }
                }
                else
                {
                    cerr << "fatal" << endl;
                }
                double ratio_beta = (InvTemp[0][k] - InvTemp[cnt_samp][k]) / (InvTemp[cnt_samp][ext_k] - InvTemp[cnt_samp][k]);
                IPL_Z = (Norm[cnt_samp][ext_k] - Norm[cnt_samp][k]) * ratio_beta + Norm[cnt_samp][k];
                IPL_Ene = (Ene[cnt_samp][ext_k] - Ene[cnt_samp][k]) * ratio_beta + Ene[cnt_samp][k];
                IPL_Ene2 = (Ene2[cnt_samp][ext_k] - Ene2[cnt_samp][k]) * ratio_beta + Ene2[cnt_samp][k];
                IPL_Sz = (Sz[cnt_samp][ext_k] - Sz[cnt_samp][k]) * ratio_beta + Sz[cnt_samp][k];
                IPL_Sz2 = (Sz2[cnt_samp][ext_k] - Sz2[cnt_samp][k]) * ratio_beta + Sz2[cnt_samp][k];
                IPL_InvTemp = InvTemp[0][k];
            }
            if (tpq_type == "cTPQ")
            {

                IPL_Z = Norm[cnt_samp][k]; // cnt_samp 和 k 都大于1
                IPL_Ene = Ene[cnt_samp][k];
                IPL_Ene2 = Ene2[cnt_samp][k];
                IPL_Sz = Sz[cnt_samp][k];
                IPL_Sz2 = Sz2[cnt_samp][k];
                IPL_InvTemp = InvTemp[cnt_samp][k];
            }

            // 归一化和累积计算
            // normalized by Z @ samp=0
            // N_k = <phi_k|phi_k>

            tot_Z = tot_Z * IPL_Z / Norm[0][k]; // 为什么是连乘呢？？？
            phys_Z[cnt_samp][k] = tot_Z;
            phys_Ene[cnt_samp][k] = IPL_Ene * tot_Z;
            phys_Ene2[cnt_samp][k] = IPL_Ene2 * tot_Z;
            phys_Sz[cnt_samp][k] = IPL_Sz * tot_Z;
            phys_Sz2[cnt_samp][k] = IPL_Sz2 * tot_Z;

            phys_InvTemp[cnt_samp][k] = IPL_InvTemp;
        }
    }

    int max_BS = NumAve;

    BS_Basic(max_BS, log_Z, phys_Z, phys_Ene, phys_Ene2, phys_Sz, phys_Sz2, phys_InvTemp);
}

void TPQ::BS_Basic(int &max_BS, const my_vec &log_Z, const my_mat &phys_Z, const my_mat &phys_Ene, const my_mat &phys_Ene2, const my_mat &phys_Sz, const my_mat &phys_Sz2, const my_mat &phys_InvTemp)
{
    int NumSample = max_BS;

    cout << "max_BS: " << max_BS << endl; // 抽样次数
    cout << "NumSample: " << NumSample << endl;

    // 初始化Bootstrap矩阵
    my_mat BS_Z(max_BS, my_vec(Lanczos_Max, 0.0));
    my_mat BS_Ene(max_BS, my_vec(Lanczos_Max, 0.0));
    my_mat BS_Spc(max_BS, my_vec(Lanczos_Max, 0.0));
    my_mat BS_Sz(max_BS, my_vec(Lanczos_Max, 0.0));
    my_mat BS_Chi(max_BS, my_vec(Lanczos_Max, 0.0));
    my_mat BS_Ent(max_BS, my_vec(Lanczos_Max, 0.0));

    std::random_device my_rd;
    std::mt19937 my_gen(my_rd());

    // 从NumAve个样本文件中选择max_BS个样本文件，共进行NumSample次抽样
    //  进行Bootstrap采样
    for (int set_BS = 0; set_BS < max_BS; set_BS++) // P
    {
        vector<int> sample_count(NumAve, 0); // 初始化计数器，每个样本的计数初始为0
        for (int cnt = 0; cnt < Lanczos_Max; cnt++)
        {
            double tmp_z = 0.0;
            double tmp_ene = 0.0;
            double tmp_ene2 = 0.0;
            double tmp_Sz = 0.0;
            double tmp_Sz2 = 0.0;

            for (int cnt_BS = 0; cnt_BS < NumSample; cnt_BS++) // M
            {
                // 每次抽样，将所有样本中相同迭代次数的值累加
                std::uniform_int_distribution<> dis(0, NumAve - 1);

                int rand_set = dis(my_gen);

                sample_count[rand_set]++; // 更新样本被抽中的计数器

                tmp_z += phys_Z[rand_set][cnt]; // 从NumAve样本中任意选择一个样本，选择NumSample次，将相同迭代次数的值加起来
                tmp_ene += phys_Ene[rand_set][cnt];
                tmp_ene2 += phys_Ene2[rand_set][cnt];
                tmp_Sz += phys_Sz[rand_set][cnt];
                tmp_Sz2 += phys_Sz2[rand_set][cnt];
            }

            BS_Z[set_BS][cnt] = tmp_z / NumSample; // 将其他物理量 A/tmp_z = <A> 即得到力学量平均值
            BS_Ene[set_BS][cnt] = tmp_ene / tmp_z;
            BS_Spc[set_BS][cnt] = (tmp_ene2 / tmp_z - pow(tmp_ene / tmp_z, 2));
            BS_Sz[set_BS][cnt] = tmp_Sz / tmp_z;
            BS_Chi[set_BS][cnt] = tmp_Sz2 / tmp_z - pow(tmp_Sz / tmp_z, 2);

            BS_Ent[set_BS][cnt] = log(tmp_z / NumSample) + tmp_ene / tmp_z * phys_InvTemp[0][cnt] + log_Z[cnt] + log(Hdim);
            //  S = (log(z) + <H>/T) /log(dH) dH是矩阵的维度
        }

        // 输出每个样本的被抽中次数
        cout << "Bootstrap set " << set_BS + 1 << ": sample counts -> ";
        for (int i = 0; i < NumAve; i++)
        {
            cout << sample_count[i] << "    ";
        }
        cout << endl;
    }

    // 计算平均值和标准差  力学量平局值除以总的抽样数
    my_vec ave_BS_Z = cla_mean(BS_Z, 0);
    my_vec err_BS_Z = cal_standard_deviation(BS_Z, 0, 1);
    my_vec ave_BS_Ene = cla_mean(BS_Ene, 0);
    my_vec err_BS_Ene = cal_standard_deviation(BS_Ene, 0, 1);
    my_vec ave_BS_Spc = cla_mean(BS_Spc, 0);
    my_vec err_BS_Spc = cal_standard_deviation(BS_Spc, 0, 1);
    my_vec ave_BS_Sz = cla_mean(BS_Sz, 0);
    my_vec err_BS_Sz = cal_standard_deviation(BS_Sz, 0, 1);
    my_vec ave_BS_Chi = cla_mean(BS_Chi, 0);
    my_vec err_BS_Chi = cal_standard_deviation(BS_Chi, 0, 1);
    my_vec ave_BS_Ent = cla_mean(BS_Ent, 0);
    my_vec err_BS_Ent = cal_standard_deviation(BS_Ent, 0, 1);
    // 输出结果到文件
    const string cTPQF = "cTPQ_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream file(cTPQF);
    for (int cnt = 0; cnt < Lanczos_Max; cnt++)
    {
        double beta = phys_InvTemp[0][cnt];
        double ave_Z = ave_BS_Z[cnt];
        double err_Z = err_BS_Z[cnt];
        double ave_Ene = ave_BS_Ene[cnt];
        double err_Ene = err_BS_Ene[cnt];
        double ave_Spc = pow(beta, 2) * ave_BS_Spc[cnt] / Nsites; // 处理与时间的关系
        double err_Spc = pow(beta, 2) * err_BS_Spc[cnt] / Nsites;
        double ave_Sz = ave_BS_Sz[cnt];
        double err_Sz = err_BS_Sz[cnt];
        double ave_Chi = beta * ave_BS_Chi[cnt] / Nsites;
        double err_Chi = beta * err_BS_Chi[cnt] / Nsites;
        double ave_Ent = ave_BS_Ent[cnt] / Nsites / log(Sdim);
        double err_Ent = err_BS_Ent[cnt] / Nsites / log(Sdim);

        if (cnt == 0)
        {
            file << "#  T   En   Cv   Cv_err   S   S_err   Chi   Chi_err   Sz   Sz_err   Z   Z_err   step" << endl;
        }

        if (cnt > 0)
        {
            file << left << fixed << setprecision(16)
                 << setw(22) << 1.0 / phys_InvTemp[0][cnt]
                 << setw(22) << ave_Ene
                 << setw(22) << ave_Spc << setw(22) << err_Spc
                 << setw(22) << ave_Ent << setw(22) << err_Ent
                 << setw(22) << ave_Chi << setw(22) << err_Chi
                 << setw(22) << ave_Sz << setw(22) << err_Sz
                 << setw(22) << ave_Z << setw(22) << err_Z
                 << setw(22) << cnt << endl;
        }
    }
}

void TPQ::Read_TPQ_file(const string &pri_file, my_mat &Norm, my_mat &InvTemp, my_mat &Ene, my_mat &Ene2, my_mat &Spc, my_mat &Sz, my_mat &Sz2, my_mat &chi)
{
    // 读取三个TPQ文件的内容

    Norm.resize(NumAve, my_vec(Lanczos_Max));    // 波函数的范数的平方
    InvTemp.resize(NumAve, my_vec(Lanczos_Max)); // 逆温度
    Ene.resize(NumAve, my_vec(Lanczos_Max));     // <E>
    Ene2.resize(NumAve, my_vec(Lanczos_Max));    // <E^2>
    Spc.resize(NumAve, my_vec(Lanczos_Max));     // 比热
    Sz.resize(NumAve, my_vec(Lanczos_Max));      //<Sz> = \sum_i(<S^z_i>) / N
    Sz2.resize(NumAve, my_vec(Lanczos_Max));     //<Sz^2>= <\sum_i(S^z_i)^2> / N
    chi.resize(NumAve, my_vec(Lanczos_Max));     // 磁化率

    for (int cnt_samp = 0; cnt_samp < NumAve; cnt_samp++)
    {
        // 处理 Norm_rand 文件
        string norm_file = pri_file + "/Norm_rand" + to_string(cnt_samp) + ".dat";

        ifstream norm_stream(norm_file);
        if (!norm_stream.is_open())
        {
            cerr << "Function Read_TPQ_file; Failed to open file: " << norm_file << endl;
            return;
        }

        string line;
        int line_number = 0;
        while (getline(norm_stream, line))
        {
            // 第一行是： # inv_temp, global_norm, global_1st_norm, step_i，舍弃
            if (line_number > 0)
            { // 跳过第一行
                stringstream ss(line);
                vector<string> tokens;
                string token;
                while (ss >> token)
                {
                    tokens.push_back(token);
                }
                if (tokens.size() > 1)
                {
                    Norm[cnt_samp][line_number - 1] = pow(stod(tokens[1]), 2);
                }
            }
            line_number++;
        }
        norm_stream.close();

        // 处理 SS_rand 文件
        string ss_file = pri_file + "/SS_rand" + to_string(cnt_samp) + ".dat";

        ifstream ss_stream(ss_file);
        if (!ss_stream.is_open())
        {
            cerr << "Function Read_TPQ_file; Failed to open file: " << ss_file << endl;
            return;
        }

        line_number = 0;
        while (getline(ss_stream, line))
        {
            if (line_number > 0)
            {
                stringstream ss(line);
                vector<string> tokens;
                string token;
                while (ss >> token)
                {
                    tokens.push_back(token);
                }
                if (tokens.size() > 3)
                {
                    InvTemp[cnt_samp][line_number - 1] = stod(tokens[0]);                                                            // 逆温度
                    Ene[cnt_samp][line_number - 1] = stod(tokens[1]);                                                                //<H>
                    Ene2[cnt_samp][line_number - 1] = stod(tokens[2]);                                                               //<H^2>
                    Spc[cnt_samp][line_number - 1] = pow(stod(tokens[0]), 2) * (stod(tokens[2]) - pow(stod(tokens[1]), 2)) / Nsites; // specific heat
                    // Spc = [<E^2> - (<E>^2)] / T^2 / Nsites
                }
            }
            line_number++;
        }
        ss_stream.close();

        // 处理 Flct_rand 文件
        string flct_file = pri_file + "/Flct_rand" + to_string(cnt_samp) + ".dat";

        ifstream flct_stream(flct_file);
        if (!flct_stream.is_open())
        {
            cerr << "Function Read_TPQ_file; Failed to open file: " << flct_file << endl;
            return;
        }

        line_number = 0;
        while (getline(flct_stream, line))
        {
            if (line_number > 0)
            {
                stringstream ss(line);
                vector<string> tokens;
                string token;
                while (ss >> token)
                {
                    tokens.push_back(token);
                }
                if (tokens.size() > 6)
                {
                    // 这里取-Sz是因为HPhi中将1当成0了
                    Sz[cnt_samp][line_number - 1] = -stod(tokens[5]);
                    Sz2[cnt_samp][line_number - 1] = stod(tokens[6]);
                    chi[cnt_samp][line_number - 1] = stod(tokens[0]) * (stod(tokens[6]) - pow(-stod(tokens[5]), 2)) / Nsites;

                    // chi = (<Sz^2> - <Sz>^2) / T / Nsites
                }
            }
            line_number++;
        }
        flct_stream.close();
    }
}

void TPQ::WpMoments_TPQ(const string &SixBodyPath, my_mat &Wp_mat)
{
    int data_num = 64;

    int num_wpFile = Lanczos_Max / ExpecInterval; // 一个NumAve中Wp文件的数目
    cout << "the num of Wp File: " << num_wpFile << endl;

    // plaquette 列表
    my_mat Plaqs;
    Mat_Read("check_plaq.dat", Plaqs);
    int NumPlaq = Plaqs.size();

    for (int cnt_samp = 0; cnt_samp < NumAve; cnt_samp++) // 遍历所有 set
    {
        for (int nFile = 0; nFile < num_wpFile; nFile++) // 忽略step = 0，即逆温度为0的一项
        {
            int steps = nFile * ExpecInterval; // 当前 step 的文件

            // SixBody 文件
            string SixBodyFile = SixBodyPath + "/zvo_SixBody_set" + to_string(cnt_samp) + "step" + to_string(steps) + ".dat";

            // 检查文件是否存在
            if (!ifstream(SixBodyFile))
            {
                cerr << "Warning: file " << SixBodyFile << " not found. Setting correlation to 0." << endl;
                continue; // 保持 SxCorr[cnt_samp][nFile] 等为 0（已初始化）
            }

            // SixBody 数据
            my_mat SixBody;
            Mat_Read(SixBodyFile, SixBody);

            // 构建基础自旋算符
            Operator Sx = build_Sx(Spin);
            Operator Sy = build_Sy(Spin);
            Operator Sz = build_Sz(Spin);

            // 构建 Wp = Sx * Sy * Sz * Sx * Sy * Sz
            Operator Wp = multiply_terms(multiply_terms(multiply_terms(multiply_terms(multiply_terms(Sx, Sy), Sz), Sx), Sy), Sz);

            double TotalWp = 0.0;

            for (int b = 0; b < NumPlaq; b++) // 遍历每个 Wp
            {
                int j = Plaqs[b][0];
                int k = Plaqs[b][1];
                int l = Plaqs[b][2];
                int m = Plaqs[b][3];
                int o = Plaqs[b][4];
                int p = Plaqs[b][5];

                double Wpreal = 0.0, Wpimag = 0.0;

                // 只传 plaquette index b
                ComSSSSSS_Corr(Nsites, Sdim, Wp, SixBody, b, Wpreal, Wpimag);

                my_cplx Wp = my_cplx(Wpreal, Wpimag);

                TotalWp += Wp.real();
            }

            TotalWp = TotalWp / NumPlaq; // 文件内所有 Wp 的平均值

            Wp_mat[cnt_samp][nFile] = TotalWp; // 填值到矩阵，行表示 set，列表示 step
        }
    }
}

void TPQ::S_Corr_TPQ(const string &OneBodyPath, my_mat &SxCorr, my_mat &SyCorr, my_mat &SzCorr)
{

    int numFile = Lanczos_Max / ExpecInterval; // 一个NumAve中OneBody文件的数目
    cout << "the num of OneBody File: " << numFile << endl;

    for (int cnt_samp = 0; cnt_samp < NumAve; cnt_samp++) // 遍历所有 set
    {
        for (int nFile = 0; nFile < numFile; nFile++) // 忽略step = 0，即逆温度为0的一项
        {
            int steps = nFile * ExpecInterval; // 当前 step 的文件

            string filename_CG = OneBodyPath + "/zvo_cisajs_set" + to_string(cnt_samp) + "step" + to_string(steps) + ".dat";

            // 检查文件是否存在
            if (!ifstream(filename_CG))
            {
                cerr << "Warning: file " << filename_CG << " not found. Setting correlation to 0." << endl;
                continue; // 保持 SxCorr[cnt_samp][nFile] 等为 0（已初始化）
            }

            my_mat OneBody;
            Mat_Read(filename_CG, OneBody);

            Operator SxOp = build_Sx(Spin);
            Operator SyOp = build_Sy(Spin);
            Operator SzOp = build_Sz(Spin);

            double SxNorm = 0.0, SyNorm = 0.0, SzNorm = 0.0;
            for (int i = 0; i < Nsites; i++)
            {
                double SxR = 0, SxI = 0;
                double SyR = 0, SyI = 0;
                double SzR = 0, SzI = 0;

                ComS_Corr(Sdim, SxOp, OneBody, i, SxR, SxI);
                ComS_Corr(Sdim, SyOp, OneBody, i, SyR, SyI);
                ComS_Corr(Sdim, SzOp, OneBody, i, SzR, SzI);

                my_cplx Sx(SxR, SxI);
                my_cplx Sy(SyR, SyI);
                my_cplx Sz(SzR, SzI);

                SxNorm += Sx.real();
                SyNorm += Sy.real();
                SzNorm += Sz.real();
            }

            SxCorr[cnt_samp][nFile] = SxNorm;
            SyCorr[cnt_samp][nFile] = SyNorm;
            SzCorr[cnt_samp][nFile] = SzNorm;
        }
    }
}

void TPQ::SS_Corr_TPQ(const string &TwoBodyPath, my_mat &SxSxCorr, my_mat &SySyCorr, my_mat &SzSzCorr)
{

    int numFile = Lanczos_Max / ExpecInterval; // 一个NumAve中TwoBody文件的数目
    cout << "the num of TwoBody File: " << numFile << endl;

    for (int cnt_samp = 0; cnt_samp < NumAve; cnt_samp++) // 遍历所有 set
    {
        for (int nFile = 0; nFile < numFile; nFile++) //
        {
            int steps = nFile * ExpecInterval; // 当前 step 的文件

            string filename_CG = TwoBodyPath + "/zvo_cisajscktalt_set" + to_string(cnt_samp) + "step" + to_string(steps) + ".dat";

            // 检查文件是否存在
            if (!ifstream(filename_CG))
            {
                cerr << "Warning: file " << filename_CG << " not found. Setting correlation to 0." << endl;
                continue; // 保持 SxCorr[cnt_samp][nFile] 等为 0（已初始化）
            }

            my_mat TwoBody;
            Mat_Read(filename_CG, TwoBody);

            Operator Sx = build_Sx(Spin);
            Operator Sy = build_Sy(Spin);
            Operator Sz = build_Sz(Spin);

            Operator SxSx = multiply_terms(Sx, Sx);
            Operator SxSy = multiply_terms(Sx, Sy);
            Operator SxSz = multiply_terms(Sx, Sz);

            Operator SySx = multiply_terms(Sy, Sx);
            Operator SySy = multiply_terms(Sy, Sy);
            Operator SySz = multiply_terms(Sy, Sz);

            Operator SzSx = multiply_terms(Sz, Sx);
            Operator SzSy = multiply_terms(Sz, Sy);
            Operator SzSz = multiply_terms(Sz, Sz);

            // s1 = sigma1 * Sdim + sigma2;    s2 = sigma3 * Sdim + sigma4
            // column0 = i * (Sdim^2 * N * Sdim^2 )+s1 * (NSites * Sdim^2) + j * Sdim^2 + s2

            my_cplx *SxSx_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SxSy_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SxSz_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SySx_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SySy_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SySz_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SzSx_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SzSy_Corr = new my_cplx[Nsites * Nsites];
            my_cplx *SzSz_Corr = new my_cplx[Nsites * Nsites];

            double TotalSxSx = 0.0;
            double TotalSySy = 0.0;
            double TotalSzSz = 0.0;

            for (int i = 0; i < Nsites; i++)
            {
                for (int j = 0; j < Nsites; j++)
                {

                    double SxSxReal = 0.0, SxSxImag = 0.0;
                    double SxSyReal = 0.0, SxSyImag = 0.0;
                    double SxSzReal = 0.0, SxSzImag = 0.0;

                    double SySxReal = 0.0, SySxImag = 0.0;
                    double SySyReal = 0.0, SySyImag = 0.0;
                    double SySzReal = 0.0, SySzImag = 0.0;

                    double SzSxReal = 0.0, SzSxImag = 0.0;
                    double SzSyReal = 0.0, SzSyImag = 0.0;
                    double SzSzReal = 0.0, SzSzImag = 0.0;

                    ComSS_Corr(Nsites, Sdim, SxSx, TwoBody, i, j, SxSxReal, SxSxImag);
                    ComSS_Corr(Nsites, Sdim, SxSy, TwoBody, i, j, SxSyReal, SxSyImag);
                    ComSS_Corr(Nsites, Sdim, SxSz, TwoBody, i, j, SxSzReal, SxSzImag);

                    ComSS_Corr(Nsites, Sdim, SySx, TwoBody, i, j, SySxReal, SySxImag);
                    ComSS_Corr(Nsites, Sdim, SySy, TwoBody, i, j, SySyReal, SySyImag);
                    ComSS_Corr(Nsites, Sdim, SySz, TwoBody, i, j, SySzReal, SySzImag);

                    ComSS_Corr(Nsites, Sdim, SzSx, TwoBody, i, j, SzSxReal, SzSxImag);
                    ComSS_Corr(Nsites, Sdim, SzSy, TwoBody, i, j, SzSyReal, SzSyImag);
                    ComSS_Corr(Nsites, Sdim, SzSz, TwoBody, i, j, SzSzReal, SzSzImag);

                    SxSx_Corr[i * Nsites + j] = my_cplx(SxSxReal, SxSxImag);
                    SxSy_Corr[i * Nsites + j] = my_cplx(SxSyReal, SxSyImag);
                    SxSz_Corr[i * Nsites + j] = my_cplx(SxSzReal, SxSzImag);

                    SySx_Corr[i * Nsites + j] = my_cplx(SySxReal, SySxImag);
                    SySy_Corr[i * Nsites + j] = my_cplx(SySyReal, SySyImag);
                    SySz_Corr[i * Nsites + j] = my_cplx(SySzReal, SySzImag);

                    SzSx_Corr[i * Nsites + j] = my_cplx(SzSxReal, SzSxImag);
                    SzSy_Corr[i * Nsites + j] = my_cplx(SzSyReal, SzSyImag);
                    SzSz_Corr[i * Nsites + j] = my_cplx(SzSzReal, SzSzImag);

                    TotalSxSx += SxSx_Corr[i * Nsites + j].real();
                    TotalSySy += SySy_Corr[i * Nsites + j].real();
                    TotalSzSz += SzSz_Corr[i * Nsites + j].real();
                }
            }

            SxSxCorr[cnt_samp][nFile] = TotalSxSx;
            SySyCorr[cnt_samp][nFile] = TotalSySy;
            SzSzCorr[cnt_samp][nFile] = TotalSzSz;

            //**************************释放内存**************************

            delete[] SxSy_Corr;
            delete[] SxSz_Corr;
            delete[] SySx_Corr;
            delete[] SySz_Corr;
            delete[] SzSx_Corr;
            delete[] SzSy_Corr;
            delete[] SxSx_Corr;
            delete[] SySy_Corr;
            delete[] SzSz_Corr;
        }
    }
}