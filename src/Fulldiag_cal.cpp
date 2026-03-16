#include <common.h>
#include <LOBPCG.h>
#include <FullDiag.h>

FullDiag::FullDiag(const string &ParaFilePath, const string file_theta, const string file_h) : Nsites(0), Sdim(0), Hdim(0), File_theta(file_theta), File_h(file_h)
{
    const string ParaFile = ParaFilePath + "/../modpara.def";

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
    cout << endl;
}

void FullDiag::Fulldiag_cal(const string &in_file, double T)
{

    double *Wps = new double[Hdim](); // 加括号进行值初始化，所有元素为 0
    double *SxCorr = new double[Hdim]();
    double *SyCorr = new double[Hdim]();
    double *SzCorr = new double[Hdim]();
    double *SxSxCorr = new double[Hdim]();
    double *SySyCorr = new double[Hdim]();
    double *SzSzCorr = new double[Hdim]();

    S_CorrFullDiag(in_file, SxCorr, SyCorr, SzCorr);
    cout << "SxCorr[0]: " << setprecision(6) << SxCorr[0] << endl;
    cout << "SyCorr[0]: " << setprecision(6) << SyCorr[0] << endl;
    cout << "SzCorr[0]: " << setprecision(6) << SzCorr[0] << endl;

    SS_CorrFullDiag(in_file, SxSxCorr, SySyCorr, SzSzCorr);
    cout << "SxSxCorr[0]: " << setprecision(6) << SxSxCorr[0] << endl;
    cout << "SySyCorr[0]: " << setprecision(6) << SySyCorr[0] << endl;
    cout << "SzSzCorr[0]: " << setprecision(6) << SzSzCorr[0] << endl;

    WpMoments_Fulldiag(in_file, Wps);

    string Energy_file = in_file + "/Eigenvalue.dat";
    my_mat matrix;
    Mat_Read(Energy_file, matrix);

    // 读取能量
    my_vec Energy(Hdim, 0.0);
    for (int j = 0; j < Hdim; j++)
    {
        Energy[j] = matrix[j][1]; // 能量在矩阵的第二列
    }

    my_vec Szs(Hdim, 0.0);
    my_vec S2s(Hdim, 0.0);
    Read_FullDiag_Phys(in_file, Szs, S2s);

    // 输出文件

    const string out_file = "Fulldiag_" + File_theta + "_" + File_h + ".dat";
    ofstream ofc(out_file);

    // 温度点生成的参数
    // double e_temp = -5; // 最低温：0.00001

    double e_temp = -5; // 最低温：0.01
    int int_temp = 0;
    double temp = -1;

    int p_N = 50;
    while (temp < T)
    {
        if (int_temp == p_N) // 控制在10^e_temp范围内生成p_N个温度点
        {
            int_temp = 0;
            e_temp += 1;
        }

        temp = (1 + (9.0 / (1.0 * p_N) * int_temp)) * pow(10, e_temp);

        double beta = 1.0 / temp;
        int_temp++;

        double Z = 0.0;
        double all_E = 0.0;
        double all_E2 = 0.0;
        double all_Wp = 0.0;
        double all_SzH = 0.0;
        double all_S2 = 0.0;

        //
        double all_Sx = 0.0, all_Sy = 0.0, all_Sz = 0.0;
        double all_Sx2 = 0.0, all_Sy2 = 0.0, all_Sz2 = 0.0;

        for (int en = 0; en < Hdim; en++) // 行号
        {

            double norm_ene = (Energy[en] - Energy[0]);
            double exp_factor = std::exp(-beta * norm_ene);

            Z += exp_factor;
            all_E += Energy[en] * exp_factor;
            all_E2 += Energy[en] * Energy[en] * exp_factor;

            all_SzH += Szs[en] * exp_factor;
            all_S2 += S2s[en] * exp_factor;

            all_Sx += SxCorr[en] * exp_factor;
            all_Sy += SyCorr[en] * exp_factor;
            all_Sz += SzCorr[en] * exp_factor;

            all_Sx2 += SxSxCorr[en] * exp_factor;
            all_Sy2 += SySyCorr[en] * exp_factor;
            all_Sz2 += SzSzCorr[en] * exp_factor;

            all_Wp += Wps[en] * exp_factor;
        }

        // 计算物理量
        double ob_H = all_E / Z;
        double ob_H2 = all_E2 / Z;
        double C_v = (ob_H2 - ob_H * ob_H) * (beta * beta) / Nsites;
        // (<E^2> - <E>^2) / T^2 / Nsites

        double SzH = all_SzH / Z;
        double S2 = all_S2 / Z;

        double Ent = (log(Z) + ob_H * beta - Energy[0] * beta) / Nsites / log(Sdim);
        // Ent = ln(Z) + <E>/T - E0/T

        double ob_Wp = all_Wp / Z / (Nsites / 2);

        double obSx = all_Sx / Z;
        double obSy = all_Sy / Z;
        double obSz = all_Sz / Z;

        double obSx2 = all_Sx2 / Z;
        double obSy2 = all_Sy2 / Z;
        double obSz2 = all_Sz2 / Z;

        double SxCHi = (obSx2 - obSx * obSx) * beta;
        double SyCHi = (obSy2 - obSy * obSy) * beta;
        double SzCHi = (obSz2 - obSz * obSz) * beta;

        //  tmp_Chi = beta*all_S2/Z
        double tmp_Chi = beta * all_S2 / Z;

        // 写入文件
        ofc << setw(22) << fixed << setprecision(16) << right << 1 / beta;              // 1
        ofc << setw(26) << fixed << setprecision(16) << right << ob_H;                  // 2
        ofc << setw(26) << fixed << setprecision(16) << right << C_v;                   // 3
        ofc << setw(26) << fixed << setprecision(16) << right << Ent;                   // 4
        ofc << setw(26) << fixed << setprecision(16) << right << SzH;                   // 5
        ofc << setw(26) << fixed << setprecision(16) << right << S2;                    // 6
        ofc << setw(26) << fixed << setprecision(16) << right << ob_Wp;                 // 7
        ofc << setw(26) << fixed << setprecision(16) << right << obSx;                  // 8
        ofc << setw(26) << fixed << setprecision(16) << right << obSy;                  // 9
        ofc << setw(26) << fixed << setprecision(16) << right << obSz;                  // 10
        ofc << setw(26) << fixed << setprecision(16) << right << obSx2;                 // 11
        ofc << setw(26) << fixed << setprecision(16) << right << obSy2;                 // 12
        ofc << setw(26) << fixed << setprecision(16) << right << obSz2;                 // 13
        ofc << setw(26) << fixed << setprecision(16) << right << obSx2 + obSy2 + obSz2; // 14
        ofc << setw(26) << fixed << setprecision(16) << right << SxCHi / Nsites;        // 15
        ofc << setw(26) << fixed << setprecision(16) << right << SyCHi / Nsites;        // 16
        ofc << setw(26) << fixed << setprecision(16) << right << SzCHi / Nsites;        // 17
        ofc << setw(26) << fixed << setprecision(16) << right << tmp_Chi / Nsites / 3;  // 18

        ofc << endl;
    }

    delete[] Wps;

    delete[] SxCorr;
    delete[] SyCorr;
    delete[] SzCorr;

    delete[] SxSxCorr;
    delete[] SySyCorr;
    delete[] SzSzCorr;

    ofc.close();
}

void FullDiag::S_CorrFullDiag(const string OneBodyPath, double *SxCorr, double *SyCorr, double *SzCorr)
{

    cout << "Calculating FullDiag <S^a_j>" << endl;

    for (int cnt_samp = 0; cnt_samp < Hdim; cnt_samp++) // 行号
    {

        string filename_CG = OneBodyPath + "/zvo_cisajs_eigen" + to_string(cnt_samp) + ".dat";
        if (!ifstream(filename_CG))
        { // 检查文件是否存在
            cerr << "Warning: file " << filename_CG << " not found. Setting correlation to 0." << endl;
            continue; // 保持数组对应位置为 0（已初始化）
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

        SxCorr[cnt_samp] = SxNorm;
        SyCorr[cnt_samp] = SyNorm;
        SzCorr[cnt_samp] = SzNorm;
    }
}

void FullDiag::SS_CorrFullDiag(const string TwoBodyPath, double *SxSxCorr, double *SySyCorr, double *SzSzCorr)
{

    cout << "Calculating FullDiag <S^a_i S^b_j>" << endl;

    for (int cnt_samp = 0; cnt_samp < Hdim; cnt_samp++) // 行号
    {

        string filename_CG = TwoBodyPath + "/zvo_cisajscktalt_eigen" + to_string(cnt_samp) + ".dat";
        if (!ifstream(filename_CG))
        { // 检查文件是否存在
            cerr << "Warning: file " << filename_CG << " not found. Setting correlation to 0." << endl;
            continue; // 保持数组对应位置为 0（已初始化）
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

        SxSxCorr[cnt_samp] = TotalSxSx;
        SySyCorr[cnt_samp] = TotalSySy;
        SzSzCorr[cnt_samp] = TotalSzSz;

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

void FullDiag::WpMoments_Fulldiag(const string &input_file, double *Wps)
{
    // plaquette 列表
    my_mat Plaqs;
    Mat_Read("check_plaq.dat", Plaqs);
    int NumPlaq = Plaqs.size();

    for (int cnt_samp = 0; cnt_samp < Hdim; cnt_samp++) // 行号
    {
        // SixBody 文件
        string SixBodyFile = input_file + "/zvo_SixBody_eigen" + to_string(cnt_samp) + ".dat";
        if (!ifstream(SixBodyFile))
        { // 检查文件是否存在
            cerr << "Warning: file " << SixBodyFile << " not found. Setting correlation to 0." << endl;
            continue; // 保持数组对应位置为 0（已初始化）
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
        for (int b = 0; b < NumPlaq; b++)
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
        Wps[cnt_samp] = TotalWp;
    }
}

void FullDiag::Read_FullDiag_Phys(const string &file_path, vector<double> &Sz, vector<double> &S2)
{
    const string inputFile = file_path + "/zvo_phys.dat";

    ifstream Infile(inputFile);
    if (!Infile.is_open())
    {
        cerr << "ERROR!!! Failed to open file: " << inputFile << endl;
        return;
    }

    string line;
    int line_number = 0;
    while (getline(Infile, line))
    {
        // 第一行是：  <H> <N> <Sz> <S2> <D>，舍弃
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
                // 取负号
                Sz[line_number - 1] = -stod(tokens[2]);
                S2[line_number - 1] = stod(tokens[3]);
            }
        }
        line_number++;
    }
    Infile.close();
}
