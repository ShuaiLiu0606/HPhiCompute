#include <common.h>
#include <LOBPCG.h>

//******************************************************************************************************

void LOBPCG::S_Corr(const string OneBodyPath)
{
    cout << "Calculating <S^a_j>" << endl;

    string filename_CG = OneBodyPath + "/zvo_cisajs_eigen" + to_string((int)Sample_num) + ".dat";
    string filename_Lanczos = OneBodyPath + "/zvo_cisajs.dat";

    ifstream file_CG(filename_CG);
    string OneBodyFile = file_CG.good() ? filename_CG : filename_Lanczos;

    ofstream SxCorrFile("SxCorr_All_" + to_string((int)Sample_num) + ".dat", ios::app);
    ofstream SyCorrFile("SyCorr_All_" + to_string((int)Sample_num) + ".dat", ios::app);
    ofstream SzCorrFile("SzCorr_All_" + to_string((int)Sample_num) + ".dat", ios::app);

    my_mat OneBody;
    Mat_Read(OneBodyFile, OneBody);

    Operator SxOp = build_Sx(Spin);
    Operator SyOp = build_Sy(Spin);
    Operator SzOp = build_Sz(Spin);

    double totSx = 0.0, totSy = 0.0, totSz = 0.0;
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

        totSx += Sx.real();
        totSy += Sy.real();
        totSz += Sz.real();

        SxCorrFile << fixed << setprecision(16) << Sx << " ";
        SyCorrFile << fixed << setprecision(16) << Sy << " ";
        SzCorrFile << fixed << setprecision(16) << Sz << " ";
    }

    cout << "Total <Sx>: " << fixed << setprecision(6) << totSx << endl;
    cout << "Total <Sy>: " << fixed << setprecision(6) << totSy << endl;
    cout << "Total <Sz>: " << fixed << setprecision(6) << totSz << endl;
    cout << endl;
}

void LOBPCG::SS_PCorr(const string TwoBodyPath)
{
    // 计算SxSx, SySy, SzSz的关联函数，通过双体格林函数zvo_cisajscktalt.dat文件，由索引规则得到每个算符所对应的项
    // TwoBodyPath: zvo_cisajscktalt.dat文件的路径
    // SxSx_Corr, SySy_Corr, SzSz_Corr: 保存关联函数的值

    const string P_corrPath = "SS_PCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream P_corrFile(P_corrPath, ios::app);

    P_corrFile << left << fixed << setprecision(12) << setfill(' ');
    P_corrFile << setw(4) << "j" << setw(4) << "k" << setw(4) << "a" << setw(32) << "<P^a_{j,k}>" << endl;

    cout << "Calculating <P^a_{j,k}> -> " << P_corrPath << " ... " << flush;

    cout << "Calculating <P^a_{j,k}>" << endl;

    string filename_CG = TwoBodyPath + "/zvo_cisajscktalt_eigen" + to_string((int)Sample_num) + ".dat";
    string filename_Lanczos = TwoBodyPath + "/zvo_cisajscktalt.dat";

    // 检查文件是否存在
    ifstream file_CG(filename_CG);

    string TwoBodyFile;

    if (file_CG.good())
    {
        TwoBodyFile = filename_CG;
    }
    else
    {
        TwoBodyFile = filename_Lanczos;
    }

    my_mat TwoBody;
    Mat_Read(TwoBodyFile, TwoBody);

    Operator Sx = build_Sx(Spin);
    Operator Sy = build_Sy(Spin);
    Operator Sz = build_Sz(Spin);

    Operator SySz = multiply_terms(Sy, Sz);
    Operator SzSy = multiply_terms(Sz, Sy);

    Operator SzSx = multiply_terms(Sz, Sx);
    Operator SxSz = multiply_terms(Sx, Sz);

    Operator SxSy = multiply_terms(Sx, Sy);
    Operator SySx = multiply_terms(Sy, Sx);

    // s1 = sigma1 * Sdim + sigma2;    s2 = sigma3 * Sdim + sigma4
    // column0 = i * (Sdim^2 * N * Sdim^2 )+s1 * (NSites * Sdim^2) + j * Sdim^2 + s2

    // 计算矢量手性的关联

    for (int i = 0; i < Nsites; i++)
    {
        for (int j = 0; j < Nsites; j++)
        {
            if (j <= i)
            {
                continue;
            }

            double SySzReal = 0.0, SySzImag = 0.0;
            double SzSyReal = 0.0, SzSyImag = 0.0;

            double SzSxReal = 0.0, SzSxImag = 0.0;
            double SxSzReal = 0.0, SxSzImag = 0.0;

            double SxSyReal = 0.0, SxSyImag = 0.0;
            double SySxReal = 0.0, SySxImag = 0.0;

            // Px = (S^y_j S^z_k - S^z_j S^y_k)
            // Py = (S^z_j S^x_k - S^x_j S^z_k)
            // Pz = (S^x_j S^y_k - S^y_j S^x_k)

            ComSS_Corr(Nsites, Sdim, SySz, TwoBody, i, j, SySzReal, SySzImag);
            ComSS_Corr(Nsites, Sdim, SzSy, TwoBody, i, j, SzSyReal, SzSyImag);

            ComSS_Corr(Nsites, Sdim, SzSx, TwoBody, i, j, SzSxReal, SzSxImag);
            ComSS_Corr(Nsites, Sdim, SxSz, TwoBody, i, j, SxSzReal, SxSzImag);

            ComSS_Corr(Nsites, Sdim, SxSy, TwoBody, i, j, SxSyReal, SxSyImag);
            ComSS_Corr(Nsites, Sdim, SySx, TwoBody, i, j, SySxReal, SySxImag);

            my_cplx PxCorr = my_cplx(SySzReal - SzSyReal, SySzImag - SzSyImag);
            my_cplx PyCorr = my_cplx(SzSxReal - SxSzReal, SzSxImag - SxSzImag);
            my_cplx PzCorr = my_cplx(SxSyReal - SySxReal, SxSyImag - SySxImag);

            P_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 0;
            P_corrFile << right << fixed << setprecision(16) << setw(24) << PxCorr.real() << right << fixed << setprecision(16) << setw(24) << PxCorr.imag() << endl; // <Phi|P^a_{j,k}|Phi>

            P_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 1;
            P_corrFile << right << fixed << setprecision(16) << setw(24) << PyCorr.real() << right << fixed << setprecision(16) << setw(24) << PyCorr.imag() << endl; // <Phi|P^a_{j,k}|Phi>

            P_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 2;
            P_corrFile << right << fixed << setprecision(16) << setw(24) << PzCorr.real() << right << fixed << setprecision(16) << setw(24) << PzCorr.imag() << endl; // <Phi|P^a_{j,k}|Phi>
        }
    }
}

void LOBPCG::SS_QCorr(const string TwoBodyPath)
{
    // 计算SxSx, SySy, SzSz的关联函数，通过双体格林函数zvo_cisajscktalt.dat文件，由索引规则得到每个算符所对应的项
    // TwoBodyPath: zvo_cisajscktalt.dat文件的路径
    // SxSx_Corr, SySy_Corr, SzSz_Corr: 保存关联函数的值

    const string Q_corrPath = "SS_QCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream Q_corrFile(Q_corrPath, ios::app);

    Q_corrFile << left << fixed << setprecision(12) << setfill(' ');
    Q_corrFile << setw(4) << "j" << setw(4) << "k" << setw(4) << "a" << setw(4) << "b" << setw(32) << "<Q^{a,b}_{j,k}>" << endl;
    cout << "Calculating <Q^{a,b}_{j,k}> -> " << Q_corrPath << " ... " << flush;

    string filename_CG = TwoBodyPath + "/zvo_cisajscktalt_eigen" + to_string((int)Sample_num) + ".dat";
    string filename_Lanczos = TwoBodyPath + "/zvo_cisajscktalt.dat";

    // 检查文件是否存在
    ifstream file_CG(filename_CG);

    string TwoBodyFile;

    if (file_CG.good())
    {
        TwoBodyFile = filename_CG;
    }
    else
    {
        TwoBodyFile = filename_Lanczos;
    }

    my_mat TwoBody;
    Mat_Read(TwoBodyFile, TwoBody);
    Operator Sx = build_Sx(Spin);
    Operator Sy = build_Sy(Spin);
    Operator Sz = build_Sz(Spin);

    Operator SxSx = multiply_terms(Sx, Sx);
    Operator SySy = multiply_terms(Sy, Sy);
    Operator SzSz = multiply_terms(Sz, Sz);

    Operator SxSy = multiply_terms(Sx, Sy);
    Operator SySx = multiply_terms(Sy, Sx);

    Operator SxSz = multiply_terms(Sx, Sz);
    Operator SzSx = multiply_terms(Sz, Sx);

    Operator SySz = multiply_terms(Sy, Sz);
    Operator SzSy = multiply_terms(Sz, Sy);

    // s1 = sigma1 * Sdim + sigma2;    s2 = sigma3 * Sdim + sigma4
    // column0 = i * (Sdim^2 * N * Sdim^2 )+s1 * (NSites * Sdim^2) + j * Sdim^2 + s2

    for (int i = 0; i < Nsites; i++)
    {
        for (int j = 0; j < Nsites; j++)
        {
            if (j <= i)
            {
                continue;
            }

            double SxSxReal = 0.0, SxSxImag = 0.0;
            double SySyReal = 0.0, SySyImag = 0.0;
            double SzSzReal = 0.0, SzSzImag = 0.0;

            double SxSyReal = 0.0, SxSyImag = 0.0;
            double SxSzReal = 0.0, SxSzImag = 0.0;
            double SySzReal = 0.0, SySzImag = 0.0;

            double SySxReal = 0.0, SySxImag = 0.0;
            double SzSxReal = 0.0, SzSxImag = 0.0;
            double SzSyReal = 0.0, SzSyImag = 0.0;

            ComSS_Corr(Nsites, Sdim, SxSx, TwoBody, i, j, SxSxReal, SxSxImag);
            ComSS_Corr(Nsites, Sdim, SySy, TwoBody, i, j, SySyReal, SySyImag);
            ComSS_Corr(Nsites, Sdim, SzSz, TwoBody, i, j, SzSzReal, SzSzImag);

            ComSS_Corr(Nsites, Sdim, SxSy, TwoBody, i, j, SxSyReal, SxSyImag);
            ComSS_Corr(Nsites, Sdim, SxSz, TwoBody, i, j, SxSzReal, SxSzImag);
            ComSS_Corr(Nsites, Sdim, SySz, TwoBody, i, j, SySzReal, SySzImag);

            ComSS_Corr(Nsites, Sdim, SySx, TwoBody, i, j, SySxReal, SySxImag);
            ComSS_Corr(Nsites, Sdim, SzSx, TwoBody, i, j, SzSxReal, SzSxImag);
            ComSS_Corr(Nsites, Sdim, SzSy, TwoBody, i, j, SzSyReal, SzSyImag);

            my_cplx SxSx = my_cplx(SxSxReal, SxSxImag);
            my_cplx SySy = my_cplx(SySyReal, SySyImag);
            my_cplx SzSz = my_cplx(SzSzReal, SzSzImag);

            my_cplx SxSy = my_cplx(SxSyReal, SxSyImag);
            my_cplx SxSz = my_cplx(SxSzReal, SxSzImag);
            my_cplx SySz = my_cplx(SySzReal, SySzImag);

            my_cplx SySx = my_cplx(SySxReal, SySxImag);
            my_cplx SzSx = my_cplx(SzSxReal, SzSxImag);
            my_cplx SzSy = my_cplx(SzSyReal, SzSyImag);

            // xx -> (2 xx - yy - zz)/3

            my_cplx QxxCorr = (2.0 * SxSx - SySy - SzSz) / 3.0;

            // yy -> (- xx + 2 yy - zz)/3
            my_cplx QyyCorr = (-SxSx + 2.0 * SySy - SzSz) / 3.0;

            // zz -> (-xx - yy + 2 zz)/3
            my_cplx QzzCorr = (-SxSx - SySy + 2.0 * SzSz) / 3.0;

            // xy -> (xy + yx)/2
            my_cplx QxyCorr = (SxSy + SySx) / 2.0;

            // xz -> (xz + zx)/2
            my_cplx QxzCorr = (SxSz + SzSx) / 2.0;

            // yz -> (yz + zy)/2
            my_cplx Qyzorr = (SySz + SzSy) / 2.0;

            // <Phi|Q^{a,b}_{j,k}|Phi>

            Q_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 0 << setw(4) << 0;
            Q_corrFile << right << fixed << setprecision(16) << setw(24) << QxxCorr.real() << right << fixed << setprecision(16) << setw(24) << QxxCorr.imag() << endl; // <Phi|P^a_{j,k}|Phi>

            Q_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 0 << setw(4) << 1;
            Q_corrFile << right << fixed << setprecision(16) << setw(24) << QxyCorr.real() << right << fixed << setprecision(16) << setw(24) << QxyCorr.imag() << endl;

            Q_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 0 << setw(4) << 2;
            Q_corrFile << right << fixed << setprecision(16) << setw(24) << QxzCorr.real() << right << fixed << setprecision(16) << setw(24) << QxzCorr.imag() << endl;

            Q_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 1 << setw(4) << 1;
            Q_corrFile << right << fixed << setprecision(16) << setw(24) << QyyCorr.real() << right << fixed << setprecision(16) << setw(24) << QyyCorr.imag() << endl; // <Phi|P^a_{j,k}|Phi>

            Q_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 1 << setw(4) << 2;
            Q_corrFile << right << fixed << setprecision(16) << setw(24) << Qyzorr.real() << right << fixed << setprecision(16) << setw(24) << Qyzorr.imag() << endl;

            Q_corrFile << left << setw(4) << i << setw(4) << j << setw(4) << 2 << setw(4) << 2;
            Q_corrFile << right << fixed << setprecision(16) << setw(24) << QzzCorr.real() << right << fixed << setprecision(16) << setw(24) << QzzCorr.imag() << endl; // <Phi|P^a_{j,k}|Phi>
        }
    }
}

void LOBPCG::SS_Corr(const string TwoBodyPath)
{

    // 计算SxSx, SySy, SzSz的关联函数，通过双体格林函数zvo_cisajscktalt.dat文件，由索引规则得到每个算符所对应的项
    // TwoBodyPath: zvo_cisajscktalt.dat文件的路径
    // SxSx_Corr, SySy_Corr, SzSz_Corr: 保存关联函数的值

    cout << "Calculating <S^a_j S^b_k>" << endl;

    string filename_CG = TwoBodyPath + "/zvo_cisajscktalt_eigen" + to_string((int)Sample_num) + ".dat";
    string filename_Lanczos = TwoBodyPath + "/zvo_cisajscktalt.dat";

    // 检查文件是否存在
    ifstream file_CG(filename_CG);

    string TwoBodyFile;

    if (file_CG.good())
    {
        TwoBodyFile = filename_CG;
    }
    else
    {
        TwoBodyFile = filename_Lanczos;
    }

    my_mat TwoBody;
    Mat_Read(TwoBodyFile, TwoBody);

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
        }
    }

    //******************保存关联函数*******************
    const string SxSxCorrPath = "SxSxCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SySyCorrPath = "SySyCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SzSzCorrPath = "SzSzCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";

    const string SxSyCorrPath = "SxSyCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SxSzCorrPath = "SxSzCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SySxCorrPath = "SySxCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SySzCorrPath = "SySzCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SzSxCorrPath = "SzSxCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SzSyCorrPath = "SzSyCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";

    ofstream SxSxCorrFile(SxSxCorrPath);
    ofstream SySyCorrFile(SySyCorrPath);
    ofstream SzSzCorrFile(SzSzCorrPath);

    ofstream SxSyCorrFile(SxSyCorrPath);
    ofstream SySxCorrFile(SySxCorrPath);
    ofstream SySzCorrFile(SySzCorrPath);
    ofstream SzSyCorrFile(SzSyCorrPath);
    ofstream SxSzCorrFile(SxSzCorrPath);
    ofstream SzSxCorrFile(SzSxCorrPath);

    my_cplx TotalSxSx = 0.0, TotalSySy = 0.0, TotalSzSz = 0.0;
    for (int i = 0; i < Nsites; i++)
    {
        for (int j = 0; j < Nsites; j++)
        {
            int inde = i * Nsites + j;

            SxSxCorrFile << fixed << setprecision(16) << SxSx_Corr[inde];
            SySyCorrFile << fixed << setprecision(16) << SySy_Corr[inde];
            SzSzCorrFile << fixed << setprecision(16) << SzSz_Corr[inde];

            SxSyCorrFile << fixed << setprecision(16) << SxSy_Corr[inde];
            SySxCorrFile << fixed << setprecision(16) << SySx_Corr[inde];
            SySzCorrFile << fixed << setprecision(16) << SySz_Corr[inde];
            SzSyCorrFile << fixed << setprecision(16) << SzSy_Corr[inde];
            SxSzCorrFile << fixed << setprecision(16) << SxSz_Corr[inde];
            SzSxCorrFile << fixed << setprecision(16) << SzSx_Corr[inde];

            if (j != Nsites - 1)
            {
                SxSxCorrFile << " ";
                SySyCorrFile << " ";
                SzSzCorrFile << " ";
                SxSyCorrFile << " ";
                SySxCorrFile << " ";
                SySzCorrFile << " ";
                SzSyCorrFile << " ";
                SxSzCorrFile << " ";
                SzSxCorrFile << " ";
            }

            TotalSxSx += SxSx_Corr[inde];
            TotalSySy += SySy_Corr[inde];
            TotalSzSz += SzSz_Corr[inde];
        }
        SxSxCorrFile << endl;
        SySyCorrFile << endl;
        SzSzCorrFile << endl;

        SxSyCorrFile << endl;
        SySxCorrFile << endl;
        SySzCorrFile << endl;
        SzSyCorrFile << endl;
        SxSzCorrFile << endl;
        SzSxCorrFile << endl;
    }

    //******************保存总关联函数*******************
    const string TotalSxSxCorrPath = "TotalSxSxCorr_All_" + to_string((int)Sample_num) + ".dat";
    const string TotalSySyCorrPath = "TotalSySyCorr_All_" + to_string((int)Sample_num) + ".dat";
    const string TotalSzSzCorrPath = "TotalSzSzCorr_All_" + to_string((int)Sample_num) + ".dat";

    ofstream TotalSxSxCorrFile(TotalSxSxCorrPath, ios::app);
    ofstream TotalSySyCorrFile(TotalSySyCorrPath, ios::app);
    ofstream TotalSzSzCorrFile(TotalSzSzCorrPath, ios::app);

    TotalSxSxCorrFile << File_theta << "    " << File_h << "    ";
    TotalSxSxCorrFile << fixed << setprecision(16) << right << TotalSxSx.real() << endl;
    TotalSySyCorrFile << File_theta << "    " << File_h << "    ";
    TotalSySyCorrFile << fixed << setprecision(16) << right << TotalSySy.real() << endl;
    TotalSzSzCorrFile << File_theta << "    " << File_h << "    ";
    TotalSzSzCorrFile << fixed << setprecision(16) << right << TotalSzSz.real() << endl;

    //******************输出总关联函数*******************
    cout << fixed << setprecision(8);
    cout << "Total <SxSx>: " << "(" << setw(11) << right << TotalSxSx.real() << "," << setw(11) << right << TotalSxSx.imag() << ")" << endl;
    cout << "Total <SySy>: " << "(" << setw(11) << right << TotalSySy.real() << "," << setw(11) << right << TotalSySy.imag() << ")" << endl;
    cout << "Total <SzSz>: " << "(" << setw(11) << right << TotalSzSz.real() << "," << setw(11) << right << TotalSzSz.imag() << ")" << endl;
    cout << endl;

    //**************************************************************************** */

    // 转换矩阵

    vector<double> a = {1 / sqrt(6), 1 / sqrt(6), -2 / sqrt(6)}; // A = 1 / sqrt(6) * (1, 1, -2)
    vector<double> b = {-1 / sqrt(2), 1 / sqrt(2), 0};           // B = 1 / sqrt(2) * (-1, 1, 0)
    vector<double> c = {1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)};  // C = 1 / sqrt(3) * (1, 1, 1)

    my_cplx *SaSa_Corr = new my_cplx[Nsites * Nsites];
    my_cplx *SbSb_Corr = new my_cplx[Nsites * Nsites];
    my_cplx *ScSc_Corr = new my_cplx[Nsites * Nsites];

    // 保存转换后的关联函数
    const string SaSaPath = "SaSaCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SbSbPath = "SbSbCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string ScScPath = "ScScCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";

    ofstream SaSaFile(SaSaPath);
    ofstream SbSbFile(SbSbPath);
    ofstream ScScFile(ScScPath);

    // 计算转换后的关联函数
    for (size_t i = 0; i < Nsites; i++)
    {
        for (size_t j = 0; j < Nsites; j++)
        {
            int inde = i * Nsites + j;

            // 计算规则

            SaSa_Corr[inde] = a[0] * a[0] * SxSx_Corr[inde] + a[1] * a[1] * SySy_Corr[inde] + a[2] * a[2] * SzSz_Corr[inde] +
                              a[0] * a[1] * SxSy_Corr[inde] + a[1] * a[0] * SySx_Corr[inde] + a[0] * a[2] * SxSz_Corr[inde] +
                              a[2] * a[0] * SzSx_Corr[inde] + a[1] * a[2] * SySz_Corr[inde] + a[2] * a[1] * SzSy_Corr[inde];
            SbSb_Corr[inde] = b[0] * b[0] * SxSx_Corr[inde] + b[1] * b[1] * SySy_Corr[inde] + b[2] * b[2] * SzSz_Corr[inde] +
                              b[0] * b[1] * SxSy_Corr[inde] + b[1] * b[0] * SySx_Corr[inde] + b[0] * b[2] * SxSz_Corr[inde] +
                              b[2] * b[0] * SzSx_Corr[inde] + b[1] * b[2] * SySz_Corr[inde] + b[2] * b[1] * SzSy_Corr[inde];
            ScSc_Corr[inde] = c[0] * c[0] * SxSx_Corr[inde] + c[1] * c[1] * SySy_Corr[inde] + c[2] * c[2] * SzSz_Corr[inde] +
                              c[0] * c[1] * SxSy_Corr[inde] + c[1] * c[0] * SySx_Corr[inde] + c[0] * c[2] * SxSz_Corr[inde] +
                              c[2] * c[0] * SzSx_Corr[inde] + c[1] * c[2] * SySz_Corr[inde] + c[2] * c[1] * SzSy_Corr[inde];

            SaSaFile << fixed << setprecision(16) << SaSa_Corr[inde];
            SbSbFile << fixed << setprecision(16) << SbSb_Corr[inde];
            ScScFile << fixed << setprecision(16) << ScSc_Corr[inde];

            if (j != Nsites - 1)
            {
                SaSaFile << " ";
                SbSbFile << " ";
                ScScFile << " ";
            }
        }
        SaSaFile << endl;
        SbSbFile << endl;
        ScScFile << endl;
    }

    SaSaFile.close();
    SbSbFile.close();
    ScScFile.close();

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

    delete[] SaSa_Corr;
    delete[] SbSb_Corr;
    delete[] ScSc_Corr;
}

void LOBPCG::SSS_Corr(const string ThreeBodyPath, int numh)
{
    cout << "Calculating <S^a_i S^b_j S^c_k>" << endl;

    ofstream out_full("SSS_Corr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat");

    ofstream out_conn("SSS_Corr_NoMag_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat");

    out_full << left << fixed << setprecision(12) << setfill(' ');
    out_full << setw(4) << "j" << setw(4) << "a" << setw(4) << "k" << setw(4) << "b" << setw(4) << "l" << setw(4) << "c" << setw(32) << "<S^a_j S^b_k S^c_l>" << endl;

    out_conn << left << fixed << setprecision(12) << setfill(' ');
    out_conn << setw(4) << "j" << setw(4) << "a" << setw(4) << "k" << setw(4) << "b" << setw(4) << "l" << setw(4) << "c" << setw(32) << "<S^a_j S^b_k S^c_l>_c" << endl;

    // ===============================
    // 读取三体文件
    // ===============================
    string ThreeBodyFile = ThreeBodyPath + "/zvo_ThreeBody_eigen" + to_string((int)Sample_num) + ".dat";

    if (!ifstream(ThreeBodyFile).good())
        ThreeBodyFile = ThreeBodyPath + "/zvo_ThreeBody.dat";

    my_mat ThreeBody;
    Mat_Read(ThreeBodyFile, ThreeBody);

    // ===============================
    // 读取单体关联
    // ===============================
    my_cplx_mat SxCoor, SyCoor, SzCoor;

    Mat_Read("SxCorr_All_" + to_string((int)Sample_num) + ".dat", SxCoor);
    Mat_Read("SyCorr_All_" + to_string((int)Sample_num) + ".dat", SyCoor);
    Mat_Read("SzCorr_All_" + to_string((int)Sample_num) + ".dat", SzCoor);

    // ===============================
    // 读取二体关联
    // ===============================
    my_cplx_mat SxSyCoor, SxSzCoor, SySxCoor, SySzCoor, SzSxCoor, SzSyCoor;

    Mat_Read("SxSyCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat", SxSyCoor);

    Mat_Read("SxSzCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat", SxSzCoor);

    Mat_Read("SySxCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat", SySxCoor);

    Mat_Read("SySzCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat", SySzCoor);

    Mat_Read("SzSxCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat", SzSxCoor);

    Mat_Read("SzSyCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat", SzSyCoor);

    // ===============================
    // 构造算符
    // ===============================
    Operator Sx = build_Sx(Spin);
    Operator Sy = build_Sy(Spin);
    Operator Sz = build_Sz(Spin);

    vector<Operator> ops = {
        multiply_terms(multiply_terms(Sx, Sy), Sz),
        multiply_terms(multiply_terms(Sx, Sz), Sy),
        multiply_terms(multiply_terms(Sy, Sx), Sz),
        multiply_terms(multiply_terms(Sy, Sz), Sx),
        multiply_terms(multiply_terms(Sz, Sx), Sy),
        multiply_terms(multiply_terms(Sz, Sy), Sx)};

    vector<tuple<int, int, int>> labels = {
        {0, 1, 2},
        {0, 2, 1},
        {1, 0, 2},
        {1, 2, 0},
        {2, 0, 1},
        {2, 1, 0}};

    my_cplx_mat *singleMats[3] = {&SxCoor, &SyCoor, &SzCoor};

    my_cplx_mat *twoBodyMats[3][3] = {
        {nullptr, &SxSyCoor, &SxSzCoor},
        {&SySxCoor, nullptr, &SySzCoor},
        {&SzSxCoor, &SzSyCoor, nullptr}};

    // ===============================
    // 主循环
    // ===============================
    for (int i = 0; i < Nsites; i++)
    {
        for (int j = i + 1; j < Nsites; j++)
        {
            for (int k = j + 1; k < Nsites; k++)
            {
                for (int m = 0; m < 6; m++)
                {
                    double r = 0.0, im = 0.0;

                    ComSSS_Corr(Nsites, Sdim, ops[m], ThreeBody, i, j, k, r, im);

                    my_cplx corr(r, im);

                    auto [a, b, c] = labels[m];

                    my_cplx sa = (*singleMats[a])[numh][i];
                    my_cplx sb = (*singleMats[b])[numh][j];
                    my_cplx sc = (*singleMats[c])[numh][k];

                    my_cplx sab = (*twoBodyMats[a][b])[i][j];
                    my_cplx sac = (*twoBodyMats[a][c])[i][k];
                    my_cplx sbc = (*twoBodyMats[b][c])[j][k];

                    my_cplx connected = corr - sa * sbc - sb * sac - sc * sab + 2.0 * sa * sb * sc;

                    out_full << left << setw(4) << i << setw(4) << a << setw(4) << j << setw(4) << b << setw(4) << k << setw(4) << c;
                    out_full << right << setprecision(16) << fixed << setw(24) << corr.real() << setw(24) << corr.imag() << endl;

                    out_conn << left << setw(4) << i << setw(4) << a << setw(4) << j << setw(4) << b << setw(4) << k << setw(4) << c;
                    out_conn << right << setprecision(16) << fixed << setw(24) << connected.real() << setw(24) << connected.imag() << endl;
                }
            }
        }
    }
}

void LOBPCG::SSSS_PPCorr(const string FourBodyPath)
{

    // 输出文件准备
    string PP_CoorPath = "SSSS_PP_Corr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream PPCorrFile(PP_CoorPath, ios::app);
    PPCorrFile << left << fixed << setprecision(12) << setfill(' ');
    PPCorrFile << setw(4) << "j" << setw(4) << "k" << setw(4) << "a" << setw(4) << "l" << setw(4) << "m" << setw(4) << "b" << setw(32) << "<P^a_{j,k} P^b_{l,m}>" << endl;
    cout << "Calculating <P^a_{j,k} P^b_{l,m}> -> " << PP_CoorPath << " ... " << endl;
    // 文件处理
    string FourBodyFile = FourBodyPath + "/zvo_FourBody_eigen" + to_string((int)Sample_num) + ".dat";
    if (!ifstream(FourBodyFile).good())
    {
        FourBodyFile = FourBodyPath + "/zvo_FourBody.dat";
    }

    const string bondF = "check_bond.dat";
    my_mat bonds;
    Mat_Read(bondF, bonds);
    int NumBond = bonds.size();

    my_mat FourBody;
    Mat_Read(FourBodyFile, FourBody);

    // 构建基础自旋算符（已包含正确归一化系数）
    Operator Sx = build_Sx(Spin);
    Operator Sy = build_Sy(Spin);
    Operator Sz = build_Sz(Spin);

    // 辅助函数：两个算符的乘积（直接返回 multiply_terms 的结果，不合并）
    auto mul = [&](const Operator &a, const Operator &b) -> Operator
    {
        return multiply_terms(a, b);
    };

    // 辅助函数：算符的线性组合 (a - b) —— 不合并重复项
    auto sub = [&](const Operator &a, const Operator &b) -> Operator
    {
        Operator result = a;
        for (const auto &term : b)
        {
            result.emplace_back(term.first, -term.second);
        }
        return result; // 允许重复项存在，后续累加时自动处理
    };

    // 构建 bond 算符 P^a
    Operator Px = sub(multiply_terms(Sy, Sz), multiply_terms(Sz, Sy)); // P^x = S^y*S^z - S^z*S^y
    Operator Py = sub(multiply_terms(Sz, Sx), multiply_terms(Sx, Sz)); // P^y = S^z*S^x - S^x*S^z
    Operator Pz = sub(multiply_terms(Sx, Sy), multiply_terms(Sy, Sx)); // P^z = S^x*S^y - S^y*S^x

    // 构建 9 个关联函数算符：P^a * P^b (a,b = 0,1,2)
    vector<vector<Operator>> PP(3, vector<Operator>(3));
    PP[0][0] = multiply_terms(Px, Px);
    PP[0][1] = multiply_terms(Px, Py);
    PP[0][2] = multiply_terms(Px, Pz);
    PP[1][0] = multiply_terms(Py, Px);
    PP[1][1] = multiply_terms(Py, Py);
    PP[1][2] = multiply_terms(Py, Pz);
    PP[2][0] = multiply_terms(Pz, Px);
    PP[2][1] = multiply_terms(Pz, Py);
    PP[2][2] = multiply_terms(Pz, Pz);

    TimeVar ti = timeNow();

    // 遍历所有键对
    for (uint32_t b = 0; b < NumBond; b++)
    {
        uint32_t j = bonds[b][1];
        uint32_t k = bonds[b][2];

        for (uint32_t bp = 0; bp < NumBond; bp++)
        {
            uint32_t l = bonds[bp][1];
            uint32_t m = bonds[bp][2];

            // 计算所有 9 个关联函数
            my_cplx results[3][3];
            for (int a = 0; a < 3; a++)
            {
                for (int b = 0; b < 3; b++)
                {
                    double real = 0.0, imag = 0.0;
                    ComSSSS_Corr(Nsites, Sdim, PP[a][b], FourBody, j, k, l, m, real, imag);
                    results[a][b] = my_cplx(real, imag);
                }
            }

            // 按原格式输出
            for (int a = 0; a < 3; a++)
            {
                for (int b = 0; b < 3; b++)
                {
                    PPCorrFile << left << setw(4) << j << setw(4) << k << setw(4) << a << setw(4) << l << setw(4) << m << setw(4) << b;
                    PPCorrFile << right << setprecision(16) << fixed << setw(24) << results[a][b].real() << setw(24) << results[a][b].imag() << endl;
                }
            }
        }
    }

    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;
    PPCorrFile.close();
}

void LOBPCG::SSSS_Corr(const string FourBodyPath)
{
    // 准备输出文件
    string SSSS_CoorPath = "SSSS_Corr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream SSSSCorrFile(SSSS_CoorPath, ios::app);
    SSSSCorrFile << left << fixed << setprecision(12) << setfill(' ');
    SSSSCorrFile << setw(4) << "j" << setw(4) << "a" << setw(4) << "k" << setw(4) << "b" << setw(4) << "l" << setw(4) << "c" << setw(4) << "m" << setw(4) << "d";
    SSSSCorrFile << setw(32) << "<S^a_j S^b_k S^c_l S^d_m >" << endl;
    cout << "Calculating <S^a_j S^b_k S^c_l S^d_m> -> " << SSSS_CoorPath << " ... " << endl;

    // 文件路径处理
    string FourBodyFile = FourBodyPath + "/zvo_FourBody_eigen" + to_string((int)Sample_num) + ".dat";
    if (!ifstream(FourBodyFile).good())
    {
        FourBodyFile = FourBodyPath + "/zvo_FourBody.dat";
    }

    const string bondF = "check_bond.dat";
    my_mat bonds;
    Mat_Read(bondF, bonds);
    int NumBond = bonds.size();

    my_mat FourBody;
    Mat_Read(FourBodyFile, FourBody);

    // 构建单个自旋算符的 Operator
    Operator Sx = build_Sx(Spin);
    Operator Sy = build_Sy(Spin);
    Operator Sz = build_Sz(Spin);

    // 定义所有四算符组合的向量
    vector<vector<Operator>> SSSS_Pairs = {
        {Sx, Sx, Sx, Sx},
        {Sx, Sx, Sx, Sy},
        {Sx, Sx, Sx, Sz},
        {Sx, Sx, Sy, Sx},
        {Sx, Sx, Sy, Sy},
        {Sx, Sx, Sy, Sz},
        {Sx, Sx, Sz, Sx},
        {Sx, Sx, Sz, Sy},
        {Sx, Sx, Sz, Sz},

        {Sx, Sy, Sx, Sx},
        {Sx, Sy, Sx, Sy},
        {Sx, Sy, Sx, Sz},
        {Sx, Sy, Sy, Sx},
        {Sx, Sy, Sy, Sy},
        {Sx, Sy, Sy, Sz},
        {Sx, Sy, Sz, Sx},
        {Sx, Sy, Sz, Sy},
        {Sx, Sy, Sz, Sz},

        {Sx, Sz, Sx, Sx},
        {Sx, Sz, Sx, Sy},
        {Sx, Sz, Sx, Sz},
        {Sx, Sz, Sy, Sx},
        {Sx, Sz, Sy, Sy},
        {Sx, Sz, Sy, Sz},
        {Sx, Sz, Sz, Sx},
        {Sx, Sz, Sz, Sy},
        {Sx, Sz, Sz, Sz},

        {Sy, Sx, Sx, Sx},
        {Sy, Sx, Sx, Sy},
        {Sy, Sx, Sx, Sz},
        {Sy, Sx, Sy, Sx},
        {Sy, Sx, Sy, Sy},
        {Sy, Sx, Sy, Sz},
        {Sy, Sx, Sz, Sx},
        {Sy, Sx, Sz, Sy},
        {Sy, Sx, Sz, Sz},

        {Sy, Sy, Sx, Sx},
        {Sy, Sy, Sx, Sy},
        {Sy, Sy, Sx, Sz},
        {Sy, Sy, Sy, Sx},
        {Sy, Sy, Sy, Sy},
        {Sy, Sy, Sy, Sz},
        {Sy, Sy, Sz, Sx},
        {Sy, Sy, Sz, Sy},
        {Sy, Sy, Sz, Sz},

        {Sy, Sz, Sx, Sx},
        {Sy, Sz, Sx, Sy},
        {Sy, Sz, Sx, Sz},
        {Sy, Sz, Sy, Sx},
        {Sy, Sz, Sy, Sy},
        {Sy, Sz, Sy, Sz},
        {Sy, Sz, Sz, Sx},
        {Sy, Sz, Sz, Sy},
        {Sy, Sz, Sz, Sz},

        {Sz, Sx, Sx, Sx},
        {Sz, Sx, Sx, Sy},
        {Sz, Sx, Sx, Sz},
        {Sz, Sx, Sy, Sx},
        {Sz, Sx, Sy, Sy},
        {Sz, Sx, Sy, Sz},
        {Sz, Sx, Sz, Sx},
        {Sz, Sx, Sz, Sy},
        {Sz, Sx, Sz, Sz},

        {Sz, Sy, Sx, Sx},
        {Sz, Sy, Sx, Sy},
        {Sz, Sy, Sx, Sz},
        {Sz, Sy, Sy, Sx},
        {Sz, Sy, Sy, Sy},
        {Sz, Sy, Sy, Sz},
        {Sz, Sy, Sz, Sx},
        {Sz, Sy, Sz, Sy},
        {Sz, Sy, Sz, Sz},

        {Sz, Sz, Sx, Sx},
        {Sz, Sz, Sx, Sy},
        {Sz, Sz, Sx, Sz},
        {Sz, Sz, Sy, Sx},
        {Sz, Sz, Sy, Sy},
        {Sz, Sz, Sy, Sz},
        {Sz, Sz, Sz, Sx},
        {Sz, Sz, Sz, Sy},
        {Sz, Sz, Sz, Sz}};

    // 标签可以重用，但不再用于系数乘法，仅用于输出标识
    const vector<tuple<int, int, int, int>> SSSSLabels = {
        // 前两个运算符为 Sx, Sx
        {0, 0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 0, 2},
        {0, 0, 1, 0},
        {0, 0, 1, 1},
        {0, 0, 1, 2},
        {0, 0, 2, 0},
        {0, 0, 2, 1},
        {0, 0, 2, 2},

        // 前两个运算符为 Sx, Sy
        {0, 1, 0, 0},
        {0, 1, 0, 1},
        {0, 1, 0, 2},
        {0, 1, 1, 0},
        {0, 1, 1, 1},
        {0, 1, 1, 2},
        {0, 1, 2, 0},
        {0, 1, 2, 1},
        {0, 1, 2, 2},

        // 前两个运算符为 Sx, Sz
        {0, 2, 0, 0},
        {0, 2, 0, 1},
        {0, 2, 0, 2},
        {0, 2, 1, 0},
        {0, 2, 1, 1},
        {0, 2, 1, 2},
        {0, 2, 2, 0},
        {0, 2, 2, 1},
        {0, 2, 2, 2},

        // 前两个运算符为 Sy, Sx
        {1, 0, 0, 0},
        {1, 0, 0, 1},
        {1, 0, 0, 2},
        {1, 0, 1, 0},
        {1, 0, 1, 1},
        {1, 0, 1, 2},
        {1, 0, 2, 0},
        {1, 0, 2, 1},
        {1, 0, 2, 2},

        // 前两个运算符为 Sy, Sy
        {1, 1, 0, 0},
        {1, 1, 0, 1},
        {1, 1, 0, 2},
        {1, 1, 1, 0},
        {1, 1, 1, 1},
        {1, 1, 1, 2},
        {1, 1, 2, 0},
        {1, 1, 2, 1},
        {1, 1, 2, 2},

        // 前两个运算符为 Sy, Sz
        {1, 2, 0, 0},
        {1, 2, 0, 1},
        {1, 2, 0, 2},
        {1, 2, 1, 0},
        {1, 2, 1, 1},
        {1, 2, 1, 2},
        {1, 2, 2, 0},
        {1, 2, 2, 1},
        {1, 2, 2, 2},

        // 前两个运算符为 Sz, Sx
        {2, 0, 0, 0},
        {2, 0, 0, 1},
        {2, 0, 0, 2},
        {2, 0, 1, 0},
        {2, 0, 1, 1},
        {2, 0, 1, 2},
        {2, 0, 2, 0},
        {2, 0, 2, 1},
        {2, 0, 2, 2},

        // 前两个运算符为 Sz, Sy
        {2, 1, 0, 0},
        {2, 1, 0, 1},
        {2, 1, 0, 2},
        {2, 1, 1, 0},
        {2, 1, 1, 1},
        {2, 1, 1, 2},
        {2, 1, 2, 0},
        {2, 1, 2, 1},
        {2, 1, 2, 2},

        // 前两个运算符为 Sz, Sz
        {2, 2, 0, 0},
        {2, 2, 0, 1},
        {2, 2, 0, 2},
        {2, 2, 1, 0},
        {2, 2, 1, 1},
        {2, 2, 1, 2},
        {2, 2, 2, 0},
        {2, 2, 2, 1},
        {2, 2, 2, 2}};

    // 预计算所有四算符乘积序列
    vector<Operator> SSSS_Results(SSSS_Pairs.size());
    for (size_t i = 0; i < SSSS_Pairs.size(); ++i)
    {
        const auto &ops = SSSS_Pairs[i];
        // 从四个 Operator 逐个相乘得到乘积序列
        Operator prod = ops[0];
        for (size_t j = 1; j < ops.size(); ++j)
        {
            prod = multiply_terms(prod, ops[j]);
        }
        SSSS_Results[i] = prod;
    }

    TimeVar ti = timeNow();
    for (uint32_t b = 0; b < NumBond; b++)
    {
        uint32_t j = bonds[b][1];
        uint32_t k = bonds[b][2];

        for (uint32_t bp = 0; bp < NumBond; bp++)
        {
            uint32_t l = bonds[bp][1];
            uint32_t m = bonds[bp][2];

            for (size_t num = 0; num < SSSS_Results.size(); num++)
            {
                // 获取标签中的分量值用于输出
                auto &ssssLabel = SSSSLabels[num];
                int a = get<0>(ssssLabel), b = get<1>(ssssLabel);
                int c = get<2>(ssssLabel), d = get<3>(ssssLabel);

                double real = 0.0, imag = 0.0;
                ComSSSS_Corr(Nsites, Sdim, SSSS_Results[num], FourBody, j, k, l, m, real, imag);
                my_cplx SSSS_Corr(real, imag); // 直接使用累加结果，无需额外系数

                // 输出结果
                SSSSCorrFile << left << setw(4) << j << setw(4) << a << setw(4) << k << setw(4) << b;
                SSSSCorrFile << left << setw(4) << l << setw(4) << c << setw(4) << m << setw(4) << d;
                SSSSCorrFile << right << setprecision(16) << fixed << setw(24) << SSSS_Corr.real() << setprecision(16) << fixed << setw(24) << SSSS_Corr.imag() << endl;
            }
        }
    }
    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;
    // 可选输出耗时
    SSSSCorrFile.close();
}

void LOBPCG::SSSSSS_WpCorr(const string SixBodyPath)
{
    // 输出文件
    string Wp_CorrPath = "Wp_Corr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream Wp_CorrFile(Wp_CorrPath);

    cout << "Calculating <WpMoments> -> " << Wp_CorrPath << " ... " << endl;

    // SixBody 文件
    string SixBodyFile = SixBodyPath + "/zvo_SixBody_eigen" + to_string((int)Sample_num) + ".dat";
    if (!ifstream(SixBodyFile).good())
    {
        SixBodyFile = SixBodyPath + "/zvo_SixBody.dat";
    }

    // plaquette 列表
    my_mat Plaqs;
    Mat_Read("check_plaq.dat", Plaqs);
    int NumPlaq = Plaqs.size();

    // SixBody 数据
    my_mat SixBody;
    Mat_Read(SixBodyFile, SixBody);

    // 构建基础自旋算符
    Operator Sx = build_Sx(Spin);
    Operator Sy = build_Sy(Spin);
    Operator Sz = build_Sz(Spin);

    // 构建 Wp = Sx * Sy * Sz * Sx * Sy * Sz
    Operator Wp = multiply_terms(multiply_terms(multiply_terms(multiply_terms(multiply_terms(Sx, Sy), Sz), Sx), Sy), Sz);

    if (!Wp_CorrFile.is_open())
    {
        cerr << "Failed to open " << Wp_CorrPath << endl;
        return;
    }

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

        Wp_CorrFile << b << setw(18) << fixed << setprecision(12) << Wp.real() << endl;
    }

    Wp_CorrFile.close();
}

void LOBPCG::S(int a, int j, my_cplx *input, my_cplx *output)
{
    // 这是使用波函数的方式来计Sx|phi>, Sy|phi>, Sz|phi>
    // input是波函数，output = Sx|input>, Sy|input>, Sz|input>
    // 在函数外部分配空间!!!!!
    // 将output数组的所有元素初始化为0
    std::fill(output, output + Hdim, my_cplx(0, 0));

    //
    assert(a == 0 || a == 1 || a == 2);

    if (a == 0) // a = 0 = Sx
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            // 1 << ind：找到S^x_j中j所对应的态的位置  state ^ (1 << ind):将state中j位置的自旋进行翻转
            uint32_t state = c;
            output[c] += 0.5 * input[state ^ (1 << j)];
        }
    }
    if (a == 1) // Sy
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            // 注意：在HPhi中，自旋向上是0，自旋向下是1
            uint32_t state = c;

            if ((state >> j) & 1) // |down> came from |up>
            {

                // 如果 （state >> j) == 1，则表明’当前‘的态,即算符作用以后的态（state = c）的 j 位置是下，那么原来的 j 位置是上
                output[c] += 0.5 * I * input[state ^ (1UL << j)];
            }
            else //  |up> came from |down>
            {
                output[c] += -0.5 * I * input[state ^ (1UL << j)];
            }
        }
    }
    if (a == 2) // Sz
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            uint32_t state = c;
            if ((state >> j) & 1) // |down>   //HPhi
            {
                output[c] += -0.5 * input[c];
            }
            else
            {
                output[c] += 0.5 * input[c];
            }
        }
    }
}

/* void LOBPCG::S(int a, int j, my_cplx *input, my_cplx *output)
{
    // 这是使用波函数的方式来计Sx|phi>, Sy|phi>, Sz|phi>
    // input是波函数，output = Sx|input>, Sy|input>, Sz|input>
    // 在函数外部分配空间!!!!!
    // 将output数组的所有元素初始化为0
    std::fill(output, output + Hdim, my_cplx(0, 0));

    //
    assert(a == 0 || a == 1 || a == 2);

    if (a == 0) // a = 0 = Sx
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            // 1 << ind：找到S^x_j中j所对应的态的位置  state ^ (1 << ind):将state中j位置的自旋进行翻转
            uint32_t state = c;
            output[state ^ (1 << j)] += 0.5 * input[c];
        }
    }
    if (a == 1) // Sy
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            // 注意：在HPhi中，自旋向上是0，自旋向下是1
            uint32_t state = c;

            if ((state >> j) & 1) // |down> came from |up>
            {

                // 如果 （state >> j) == 1，则表明’当前‘的态（state = c）的 j 位置是下，那么原来的 j 位置是上
                output[state ^ (1UL << j)] += -0.5 * I * input[c];
            }
            else //  |up> came from |down>
            {
                output[state ^ (1UL << j)] += +0.5 * I * input[c];
            }
        }
    }
    if (a == 2) // Sz
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            uint32_t state = c;
            if ((state >> j) & 1) // |down>   //HPhi
            {
                output[c] += -0.5 * input[c];
            }
            else
            {
                output[c] += 0.5 * input[c];
            }
        }
    }
} */

void LOBPCG::SpSm(int a, int j, my_cplx *input, my_cplx *output)
{
    // 将 output 数组的所有元素初始化为 0
    std::fill(output, output + Hdim, my_cplx(0, 0));

    // 检查 a 的值是否合法
    assert(a == 0 || a == 1);

    if (a == 0)
    { // S⁺ (升算符)
#pragma omp parallel for
        for (int c = 0; c < Hdim; c++)
        {
            int state = c;
            if ((state >> j) & 1)
            {                                           // 如果是 1（下自旋），才能执行 S⁺
                output[c] += input[state ^ (1UL << j)]; // 翻转 j 号自旋
            }
        }
    }
    else if (a == 1)
    { // S⁻ (降算符)
#pragma omp parallel for
        for (int c = 0; c < Hdim; c++)
        {
            int state = c;
            if (!((state >> j) & 1))
            {                                           // 如果是 0（上自旋），才能执行 S⁻
                output[c] += input[state ^ (1UL << j)]; // 翻转 j 号自旋
            }
        }
    }
    else
    {
        std::cerr << "Error: SpSm; a must be 0 or 1" << std::endl;
        std::exit(EXIT_FAILURE); // 终止程序
    }
}

void LOBPCG::Sq(int a, int j, my_cplx *input, my_cplx *output, vector<vector<double>> R, pair<double, double> q)
{
    // 这是使用波函数的方式来计算Sq^x|phi>, Sq^y|phi>, Sq^z|phi>  Sq^x = e^(iq*r)Sx Sq^y = e^(iq*r)Sy Sq^z = e^(iq*r)Sz
    // input是波函数，output = Sq^x|input>, Sq^y|input>, Sq^z|input>
    // 在函数外部分配空间!!!!!

    // 将output数组的所有元素初始化为0
    std::fill(output, output + Hdim, my_cplx(0, 0));

    assert(a == 0 || a == 1 || a == 2);

    double rx = R[j][1];
    double ry = R[j][2];
    double qx_point = q.first;
    double qy_point = q.second;
    double phase = qx_point * rx + qy_point * ry;             // q * r_i
    complex<double> e_iqr = exp(complex<double>(0.0, phase)); // e^[iq * r_i - r_j]

    if (a == 0) // Sx
    {
#pragma omp parallel for
        for (int c = 0; c < Hdim; c++)
        {

            int state = c;
            output[c] += 0.5 * input[state ^ (1U << j)] * e_iqr;
        }
    }
    if (a == 1) // Sy
    {
#pragma omp parallel for
        for (int c = 0; c < Hdim; c++)
        {
            int state = c;
            if ((state >> j) & 1) //|down>
            {
                output[c] += -0.5 * I * input[state ^ (1UL << j)] * e_iqr;
            }
            else
            {
                output[c] += 0.5 * I * input[state ^ (1UL << j)] * e_iqr;
            }
        }
    }
    if (a == 2) // Sz
    {
#pragma omp parallel for
        for (int c = 0; c < Hdim; c++)
        {
            int state = c;
            if ((state >> j) & 1) // |down>
            {
                output[c] += -0.5 * input[c] * e_iqr;
            }
            else
            {
                output[c] += +0.5 * input[c] * e_iqr;
            }
        }
    }
}

void LOBPCG::S_Moments(my_cplx *Phi)
{
    // 通过波函数计算 <S^a_j> = <Phi|S^a_j|Phi>
    // a = 0, 1, 2 分别对应 Sx, Sy, Sz
    // Phi是波函数，Sx_Corr, Sy_Corr, Sz_Corr是输出的结果
    //<Phi|S^a_j|Phi>

    const string SxCorrPath = "SxMoments_All_" + to_string((int)Sample_num) + ".dat";
    const string SyCorrPath = "SyMoments_All_" + to_string((int)Sample_num) + ".dat";
    const string SzCorrPath = "SzMoments_All_" + to_string((int)Sample_num) + ".dat";

    cout << "Calculating <S^a_j> -> " << SxCorrPath << "; " << SyCorrPath << "; " << SzCorrPath << " ... " << flush;

    ofstream SxCorrFile(SxCorrPath, ios::app);
    ofstream SyCorrFile(SyCorrPath, ios::app);
    ofstream SzCorrFile(SzCorrPath, ios::app);

    TimeVar ti = timeNow();

    my_cplx *V = new my_cplx[Hdim];
    for (int j = 0; j < Nsites; j++)
    {
        for (int a = 0; a < 3; a++)
        {
            // 将V数组的所有元素初始化为0
            fill(V, V + Hdim, my_cplx(0.0, 0.0));

            S(a, j, Phi, V);                         // V = S^a_j|Phi>
            my_cplx z = Vec_Dot(Hdim, Phi, 1, V, 1); // <Phi|S^a_j|Phi>

#pragma omp critical
            {
                if (a == 0) // Sx
                {
                    SxCorrFile << fixed << setprecision(16) << right << z << "    ";
                }
                else if (a == 1) // Sy
                {

                    SyCorrFile << fixed << setprecision(16) << right << z << "    ";
                }
                else if (a == 2) // Sz
                {
                    SzCorrFile << fixed << setprecision(16) << right << z << "    ";
                }
            }
        }
    }

    SxCorrFile << endl;
    SyCorrFile << endl;
    SzCorrFile << endl;

    delete[] V;

    SxCorrFile.close();
    SyCorrFile.close();
    SzCorrFile.close();

    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;
    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::SS_Correlators(my_cplx *Phi)
{
    // 计算 <S^a_j S^b_k> = <Phi|S^a_j S^b_k|Phi>
    // a, b = 0, 1, 2 分别对应 Sx, Sy, Sz
    // Phi是波函数，SxSx_Corr, SySy_Corr, SzSz_Corr是输出的结果

    const string SxSxCorrPath = "SxSx_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SySyCorrPath = "SySy_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SzSzCorrPath = "SzSz_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";

    const string SxSyCorrPath = "SxSy_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SySzCorrPath = "SySz_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SxSzCorrPath = "SxSz_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SySxCorrPath = "SySx_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SzSyCorrPath = "SzSy_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    const string SzSxCorrPath = "SzSx_Moments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";

    cout << "Calculating <S^a_j S^b_k> -> " << SxSxCorrPath << "; " << SxSyCorrPath << "; " << SxSzCorrPath << "; " << endl;
    cout << "                             " << SySxCorrPath << "; " << SySyCorrPath << "; " << SySzCorrPath << "; " << endl;
    cout << "                             " << SzSxCorrPath << "; " << SzSyCorrPath << "; " << SzSzCorrPath << " ... " << flush;

    ofstream SxSxCorrFile(SxSxCorrPath);
    ofstream SySyCorrFile(SySyCorrPath);
    ofstream SzSzCorrFile(SzSzCorrPath);

    ofstream SxSyCorrFile(SxSyCorrPath);
    ofstream SySzCorrFile(SySzCorrPath);
    ofstream SxSzCorrFile(SxSzCorrPath);
    ofstream SySxCorrFile(SySxCorrPath);
    ofstream SzSyCorrFile(SzSyCorrPath);
    ofstream SzSxCorrFile(SzSxCorrPath);

    my_cplx *V1 = new my_cplx[Hdim];
    my_cplx *V2 = new my_cplx[Hdim];
    my_cplx *V3 = new my_cplx[Hdim];

    TimeVar ti = timeNow();

    for (int j = 0; j < Nsites; j++)
    {
        for (int k = 0; k < Nsites; k++)
        {
            for (int a = 0; a < 3; a++)
            {
                for (int b = 0; b < 3; b++)
                {

                    fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
                    fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
                    fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));

                    S(b, k, Phi, V2); // V2 = S^b_k|Phi>
                    S(a, j, V2, V3);  // V3 = S^a_j S^b_k|Phi>

                    my_cplx z1 = Vec_Dot(Hdim, Phi, 1, V3, 1); // <Phi|S^a_j S^b_k|Phi>

#pragma omp critical
                    {
                        if (a == 0 && b == 0) // Sx*Sx
                        {

                            SxSxCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SxSxCorrFile << " ";
                            }
                        }
                        else if (a == 1 && b == 1) // Sy*Sy
                        {

                            SySyCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SySyCorrFile << " ";
                            }
                        }
                        else if (a == 2 && b == 2) // Sz*Sz
                        {

                            SzSzCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SzSzCorrFile << " ";
                            }
                        }
                        else if (a == 0 && b == 1) // Sx*Sy
                        {

                            SxSyCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SxSyCorrFile << " ";
                            }
                        }
                        else if (a == 1 && b == 0) // Sy*Sx
                        {

                            SySxCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SySxCorrFile << " ";
                            }
                        }
                        else if (a == 0 && b == 2) // Sx*Sz
                        {

                            SxSzCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SxSzCorrFile << " ";
                            }
                        }
                        else if (a == 2 && b == 0) // Sz*Sx
                        {

                            SzSxCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SzSxCorrFile << " ";
                            }
                        }
                        else if (a == 1 && b == 2) // Sy*Sz
                        {

                            SySzCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SySzCorrFile << " ";
                            }
                        }
                        else if (a == 2 && b == 1) // Sz*Sy
                        {

                            SzSyCorrFile << fixed << setprecision(16) << z1;
                            if (k != Nsites - 1)
                            {
                                SzSyCorrFile << " ";
                            }
                        }
                    }
                }
            }
        }

#pragma omp critical
        {
            SxSxCorrFile << endl;
            SySyCorrFile << endl;
            SzSzCorrFile << endl;

            SxSyCorrFile << endl;
            SxSzCorrFile << endl;
            SySzCorrFile << endl;

            SySxCorrFile << endl;
            SzSxCorrFile << endl;
            SzSyCorrFile << endl;
        }
    }
    delete[] V1;
    delete[] V2;
    delete[] V3;

    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::ChiralityMoments(my_cplx *Phi)
{

    const string SSS_FilePath = "SSS_Chirality_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream SSS_File(SSS_FilePath, ios::app);

    const string SSS_NoMagFilePath = "SSS_Chirality_NoMag_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream SSS_NoMagFile(SSS_NoMagFilePath, ios::app);

    SSS_File << left << fixed << setprecision(12) << setfill(' ');
    SSS_File << setw(4) << "j" << setw(4) << "a" << setw(4) << "k" << setw(4) << "b" << setw(4) << "l" << setw(4) << "c" << setw(32) << "<S^a_j S^b_k S^c_l>" << endl;
    cout << "Calculating <S^a_j S^b_k S^c_l> -> " << SSS_FilePath << " ... " << flush;

    SSS_NoMagFile << left << fixed << setprecision(12) << setfill(' ');
    SSS_NoMagFile << setw(4) << "j" << setw(4) << "a" << setw(4) << "k" << setw(4) << "b" << setw(4) << "l" << setw(4) << "c" << setw(32) << "<S^a_j S^b_k S^c_l>_c" << endl;
    cout << "Calculating <S^a_j S^b_k S^c_l>_c -> " << SSS_NoMagFilePath << " ... " << flush;

    my_cplx *V1 = new my_cplx[Hdim];
    my_cplx *V2 = new my_cplx[Hdim];
    my_cplx *V3 = new my_cplx[Hdim];
    my_cplx *V4 = new my_cplx[Hdim];
    my_cplx *V5 = new my_cplx[Hdim];
    my_cplx *V6 = new my_cplx[Hdim];
    my_cplx *V7 = new my_cplx[Hdim];

    TimeVar ti = timeNow();
    for (uint32_t j = 0; j < Nsites; j++)
    {
        for (uint32_t k = j + 1; k < Nsites; k++)
        {
            for (uint32_t l = k + 1; l < Nsites; l++) // Nc3 = N(N-1)(N-2)/3! terms, 2024 for 24 sites
            {
                uint32_t ind[] = {0, 1, 2}; // {a,b,c} permutations, 3! = 6 terms, N(N-1)(N-2) total, 12,144 for 24 sites
                do
                {
                    fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
                    fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
                    fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));
                    fill(V4, V4 + Hdim, my_cplx(0.0, 0.0));
                    fill(V5, V5 + Hdim, my_cplx(0.0, 0.0));
                    fill(V6, V6 + Hdim, my_cplx(0.0, 0.0));
                    fill(V7, V7 + Hdim, my_cplx(0.0, 0.0));

                    S(ind[0], j, Phi, V1); // S^a_j|Phi>
                    S(ind[1], k, Phi, V2); // S^b_k|Phi>

                    S(ind[2], l, Phi, V3); // S^c_l|Phi>

                    S(ind[0], j, V2, V4); // S^a_j S^b_k|Phi>

                    S(ind[0], j, V3, V5);                              // S^a_j S^c_l|Phi>
                    S(ind[1], k, V3, V6);                              // S^b_k S^c_l|Phi>
                    S(ind[0], j, V6, V7);                              // S^a_j S^b_k S^c_l|Phi>
                    complex<double> z1 = Vec_Dot(Hdim, Phi, 1, V7, 1); // <Phi|S^a_j S^b_k S^c_l|Phi>
                                                                       // <S^a_j S^b_k S^c_l> - <S^a_j><S^b_k S^c_l> - <S^b_k><S^a_j S^c_l> - <S^c_l><S^a_j S^b_k> + 2<S^a_j><S^b_k><S^c_l>
                    complex<double> z2 = z1 - (Vec_Dot(Hdim, Phi, 1, V1, 1)) * (Vec_Dot(Hdim, Phi, 1, V6, 1)) -
                                         (Vec_Dot(Hdim, Phi, 1, V2, 1)) * (Vec_Dot(Hdim, Phi, 1, V5, 1)) -
                                         (Vec_Dot(Hdim, Phi, 1, V3, 1)) * (Vec_Dot(Hdim, Phi, 1, V4, 1)) +
                                         2.0 * (Vec_Dot(Hdim, Phi, 1, V1, 1)) * (Vec_Dot(Hdim, Phi, 1, V2, 1)) * (Vec_Dot(Hdim, Phi, 1, V3, 1));
#pragma omp critical
                    {
                        SSS_File << left << setw(4) << j << setw(4) << ind[0] << setw(4) << k << setw(4) << ind[1] << setw(4) << l << setw(4) << ind[2];
                        SSS_File << right << setprecision(16) << fixed << setw(24) << z1.real() << setprecision(16) << fixed << setw(24) << z1.imag() << endl; // <Phi|S^a_j S^b_k S^c_l|Phi>

                        SSS_NoMagFile << left << setw(4) << j << setw(4) << ind[0] << setw(4) << k << setw(4) << ind[1] << setw(4) << l << setw(4) << ind[2];
                        SSS_NoMagFile << right << setprecision(16) << fixed << setw(24) << z2.real() << setprecision(16) << fixed << setw(24) << z2.imag() << endl; // <Phi|S^a_j S^b_k S^c_l|Phi>_c
                    }
                } while (std::next_permutation(ind, ind + 3));
            }
        }
    }
    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    delete[] V1;
    delete[] V2;
    delete[] V3;
    delete[] V4;
    delete[] V5;
    delete[] V6;
    delete[] V7;

    SSS_File.close();
    SSS_NoMagFile.close();

    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::Wp(vector<int> Plaquette, my_cplx *input, my_cplx *output)
{
    //  Wp = S^x_j S^y_k S^z_l S^ ,通过波函数计算 output = Wp|input>
    // Plaquette是一个长度为6的数组，包含了 plaquette 的 6 个格点的索引
    // input是波函数，output是输出的结果

    if (Plaquette.size() != 6)
    {
        std::cerr << "Function Wp, Plaquette.size() != 6" << endl;
        std::exit(EXIT_FAILURE);
    }

    // 将output数组的所有元素初始化为0
    fill(output, output + Hdim, my_cplx(0.0, 0.0));

    // Kitaev 模型惯例 x, y, z, x, y, z
    int j, k, l, m, n, o;
    j = Plaquette[0]; // S_j^x
    k = Plaquette[1]; // S_k^y
    l = Plaquette[2]; // S_l^z
    m = Plaquette[3]; // S_m^x
    n = Plaquette[4]; // S_m^y
    o = Plaquette[5]; // S_o^z

#pragma omp parallel for
    for (uint32_t c = 0; c < Hdim; c++)
    {
        uint32_t state = c;
        // ^ 异或：相同为0，不同为1（翻转操作）
        // & 位与：同1为1， 否则为0 （检测整数 c 的第 k 位是否为 1）
        //__builtin_popcount 是一个GCC内置函数，用于计算整数中二进制 1 的个数。
        // 1UL << k 将k表示为2进制(将1左移k位)
        // 如果他们作用的 4 个态中4个自旋向上，那么会产生负号
        // 如果他们作用的 4 个态中3个自旋向上，那么会产生正号
        // 如果他们作用的 4 个态中2个自旋向上，那么会产生负号
        // 如果他们作用的 4 个态中1个自旋向上，那么会产生正号
        // 这个表达式的结果是一个只在位置 k、l、n、o 上为 1，其余位为 0 的整数。
        // c & ((1UL << k) ^ (1UL << l) ^ (1UL << n) ^ (1UL << o))。判断c上的k,l，n,o位的自旋向上还是向下。奇数，产生-号
        // 只有在 c 的相应位位置（即 k、l、n、o）上为 1 时，结果的这些位才会是 1，其余位为 0。(即找到对应的态)
        uint32_t count = (__builtin_popcount(c & ((1UL << k) ^ (1UL << l) ^ (1UL << n) ^ (1UL << o)))) % 2; // S^y_k S^z_l S^y_n S^z_o
        // 偶数 ：count=0(负号)，奇数：count=1(正号)
        // 计算态中与 S^y 和 S^z 操作相关的位置上自旋为 1 的个数，然后判断奇偶。
        // 更新 output 数组。符号 (-1) 由 count 决定。
        // state ^ ((1UL << j) 翻转state中的第 j 位
        output[c] += -1.0 * (1.0 - 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k) ^ (1UL << m) ^ (1UL << n))]; // (i)^2(-1)^2,  S^x_j S^y_k S^x_m S^y_n

        // 只有Sx和Sy会产生翻转。这里似乎考虑的是泡利算符
        // output[c] += pow(0.5,6)*(-1.0 * (1.0 - 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k) ^ (1UL << m) ^ (1UL << n))]);
    }
}

void LOBPCG::WpMoments(my_cplx *Phi)
{
    // 计算 <Wp_j> = <Phi|Wp_j|Phi>

    vector<vector<int>> Plaqs;
    const string PlaqsF = "check_plaq.dat";
    Mat_Read(PlaqsF, Plaqs);

    const string WpPath = "WpMoments_All_" + to_string((int)Sample_num) + ".dat";
    ofstream WpFile(WpPath, ios::app);

    cout << "Calculating <WpMoments> -> " << WpPath << " ... " << flush;
    TimeVar ti = timeNow();

    cout << endl;

    if (!WpFile.is_open())
    {
        cerr << "Function WpMoments; Failed to open:" << WpPath << endl;
        return;
    }

    int Num_Plaq = Plaqs.size();

    WpFile << File_theta << "    " << File_h << "    ";

    my_cplx *V = new my_cplx[Hdim];
    for (int j = 0; j < Num_Plaq; j++)
    {
        fill(V, V + Hdim, my_cplx(0.0, 0.0));

        Wp(Plaqs[j], Phi, V);                    // W_j|Phi>
        my_cplx z = Vec_Dot(Hdim, Phi, 1, V, 1); // <Phi|W_j|Phi>

#pragma omp critical
        {
            cout << "Wp: " << "(" << setprecision(8) << fixed << right << setw(11) << z.real() << "," << setw(11) << z.imag() << ")" << endl;

            WpFile << fixed << setprecision(16) << setw(24) << right << z.real(); // <Phi|W_j|Phi>
        }
    }

    WpFile << endl; // 为了下一个文件，每个文件占据一行
    delete[] V;

    WpFile.close();

    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::WpWpCorrelators(my_cplx *input)
{
    // 通过波函数计算 <Wp_j Wp_k> = <Phi|Wp_j Wp_k|Phi>
    // input是波函数，PlaqsList是 plaquette 的索引，WpWp是输出的结果

    vector<vector<int>> PlaqsList;
    const string PlaqsF = "check_plaq.dat";
    Mat_Read(PlaqsF, PlaqsList);

    const string WpWpPath = "WpWpMoments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream WpWpFile(WpWpPath, ios::app);

    const string WpWpPathNoMag = "WpWpMomentsNomag_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream WpWpFileNoMag(WpWpPathNoMag);

    cout << "Calculating <Wp_j Wp_k> -> " << WpWpPath << "; " << WpWpPathNoMag << " ... " << flush;

    if (!WpWpFile.is_open() || !WpWpFileNoMag.is_open())
    {
        cerr << "Function WpWpCorrelators, Failed to open: " << WpWpPath << "; " << WpWpPathNoMag << endl;
        return;
    }

    int NumPlaqsList = PlaqsList.size();
    my_cplx *V1 = new my_cplx[Hdim];
    my_cplx *V2 = new my_cplx[Hdim];
    my_cplx *V3 = new my_cplx[Hdim];

    TimeVar ti2 = timeNow();

    for (int j = 0; j < NumPlaqsList; j++)
    {
        for (int k = 0; k < NumPlaqsList; k++)
        {

            fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
            fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
            fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));

            Wp(PlaqsList[j], input, V1); // W_j|input>

            Wp(PlaqsList[k], input, V2); // W_k|input>
            Wp(PlaqsList[j], V2, V3);    // W_j W_k|input>

            my_cplx WpWpVal = Vec_Dot(Hdim, input, 1, V3, 1); // <input|W_j W_k|input>

            my_cplx WpWpVal_Nomag = WpWpVal - (Vec_Dot(Hdim, input, 1, V1, 1) * Vec_Dot(Hdim, input, 1, V2, 1));
            // <input|W_j W_k|input> - <input|W_j|input><input|W_k|input>

#pragma omp critical
            {

                WpWpFile << fixed << right << setprecision(16) << WpWpVal;            // <input|W_j W_k|input>
                WpWpFileNoMag << fixed << right << setprecision(16) << WpWpVal_Nomag; // <input|W_j W_k|input> - <input|W_j|input><input|W_k|input>
            }
        }
#pragma omp critical
        {
            WpWpFile << endl;
            WpWpFileNoMag << endl;
        }
    }

    delete[] V1;
    delete[] V2;
    delete[] V3;

    WpWpFile.close();
    WpWpFileNoMag.close();

    TimeVar tf2 = timeNow();
    auto diff2 = durationms(tf2 - ti2) / 1000.0;

    cout << "done. Took " << left << fixed << setprecision(2) << diff2 << " s." << endl;
    cout << endl;
}
void LOBPCG::Wx(int num, vector<int> Wx_circuit, my_cplx *input, my_cplx *output)
{
    // num是Wx_circuit的长度，Wx_circuit是一个数组，包含了Wx的索引
    //  Wx = S^x_j S^x_k S^x_l S^x_m ..., 通过波函数计算 output = Wx|input>
    //  input是波函数，output是输出的结果

    fill(output, output + Hdim, my_cplx(0.0, 0.0));

    int xor_mask = 0;
    for (int i = 0; i < num; i++)
    {
        // ^ 异或：相同为0，不同为1（翻转操作）
        // 1UL << Wx_circuit[i] 将1左移Wx_circuit[i]位
        // xor_mask = xor_mask ^ (1UL << Wx_circuit[i])，将xor_mask和(1UL << Wx_circuit[i])异或

        xor_mask ^= (1UL << Wx_circuit[i]);
    }

#pragma omp parallel for
    for (int c = 0; c < Hdim; c++)
    {
        int state = c;

        output[c] += input[state ^ xor_mask];
    }
}

void LOBPCG::WxMoments(my_cplx *Phi, vector<vector<int>> Wxs)
{
    // 通过波函数计算 <Wx_j> = <Phi|Wx_j|Phi>
    // Phi是波函数，Wxs是 Wx 的索引

    for (vector<int> wx : Wxs)
    {
        for (int site : wx)
        {
            if (site < 0 || site >= Nsites)
            {
                cerr << "Function WxMoments, site index out of range: " << site << endl;
                return;
            }
        }
    }

    const string WxPath = "Wx_All_" + to_string((int)Sample_num) + ".dat";
    ofstream WxFile(WxPath, ios::app);

    if (!WxFile.is_open())
    {
        cerr << "Function WxMoments, Failed to open: " << WxPath << endl;
        return;
    }

    //
    int Num_Wx = Wxs.size();
    cout << "the num of Wx: " << Num_Wx << endl;

    WxFile << File_theta << "    " << File_h << "    " << flush;

    my_cplx *V = new my_cplx[Hdim];
    for (int j = 0; j < Num_Wx; j++)
    {

        fill(V, V + Hdim, my_cplx(0.0, 0.0));

        int num = Wxs[j].size();

        Wx(num, Wxs[j], Phi, V); // W_j|Phi>

        my_cplx z = Vec_Dot(Hdim, Phi, 1, V, 1); // <Phi|W_j|Phi>

#pragma omp critical
        {
            cout << "Wx: " << "(" << setprecision(8) << fixed << right << setw(11) << -z.real() << "," << setw(11) << -z.imag() << ")" << endl;

            WxFile << fixed << setw(24) << right << setprecision(16) << -z.real(); // <Phi|W_j|Phi>
        }
    }

    delete[] V;
    WxFile << endl;

    WxFile.close();

    cout << endl;
}

void LOBPCG::Wy(int num, vector<int> Wy_circuit, my_cplx *input, my_cplx *output)
{
    // num是Wy_circuit的长度，Wy_circuit是一个数组，包含了Wy的索引
    //  Wy = S^y_j S^y_k S^y_l S^y_m ..., 通过波函数计算 output = Wy|input>
    //  input是波函数，output是输出的结果

    // num是Wy_circuit的长度，必须是偶数
    assert(num % 2 == 0);

    fill(output, output + Hdim, my_cplx(0.0, 0.0));

    // 计算 Wy_circuit 中涉及的位操作的结果
    int xor_mask = 0;
    for (int i = 0; i < num; i++)
    {
        // ^ 异或：相同为0，不同为1（翻转操作）
        // 1UL << Wy_circuit[i] 将1左移Wy_circuit[i]位
        // xor_mask = xor_mask ^ (1UL << Wy_circuit[i])，将xor_mask和(1UL << Wy_circuit[i])异或

        xor_mask ^= (1UL << Wy_circuit[i]);
    }

#pragma omp parallel for
    for (int c = 0; c < Hdim; c++)
    {
        int state = c;
        int count = (__builtin_popcount(c & xor_mask)) % 2;

        if (num % 4 == 2) // Wy_circuit数目为 2,6,10... 取偶数，产生负号
        {
            output[c] += -1.0 * (1.0 - 2.0 * count) * input[state ^ xor_mask];
        }
        else if (num % 4 == 0) // Wy_circuit数目为 4,8,12 取奇数，产生正号
        {
            output[c] += (1.0 - 2.0 * count) * input[state ^ xor_mask];
        }
    }
}

void LOBPCG::WyMoments(my_cplx *Phi, vector<vector<int>> Wys)
{
    // 计算 <Wy_j> = <Phi|Wy_j|Phi>
    // Phi是波函数，Wys是 Wy 的索引
    for (vector<int> wy : Wys)
    {
        for (int site : wy)
        {
            if (site < 0 || site >= Nsites)
            {
                cerr << "Function WyMoments, site index out of range: " << site << endl;
                return;
            }
        }
    }

    const string WyPath = "Wy_All_" + to_string((int)Sample_num) + ".dat";
    ofstream WyFile(WyPath, ios::app);

    if (!WyFile.is_open())
    {
        cerr << "Function WyMoments, Failed to open: " << WyPath << endl;
        return;
    }

    int Num_Wy = Wys.size();
    cout << "the num of Wy: " << Num_Wy << endl;
    my_cplx *V = new my_cplx[Hdim];

    WyFile << File_theta << "    " << File_h << "    ";
    for (int j = 0; j < Num_Wy; j++)
    {
        fill(V, V + Hdim, my_cplx(0.0, 0.0));

        int num = Wys[j].size();

        Wy(num, Wys[j], Phi, V);                 // W_j|Phi>
        my_cplx z = Vec_Dot(Hdim, Phi, 1, V, 1); // <Phi|W_j|Phi>

#pragma omp critical
        {
            cout << "Wy: " << "(" << setprecision(8) << fixed << right << setw(11) << -z.real() << "," << setw(11) << -z.imag() << ")" << endl;

            WyFile << fixed << setw(24) << right << setprecision(16) << -z.real(); // <Phi|W_j|Phi>
        }
    }

    WyFile << endl;
    delete[] V;

    cout << endl;
}

void LOBPCG::Wz(int num, vector<int> Wz_circuit, my_cplx *input, my_cplx *output)
{
    // num是Wz_circuit的长度，Wz_circuit是一个数组，包含了Wz的索引
    //  Wz = S^z_j S^z_k S^z_l S^z_m ..., 通过波函数计算 output = Wz|input>
    //  input是波函数，output是输出的结果

    fill(output, output + Hdim, my_cplx(0.0, 0.0));

    int xor_mask = 0;
    for (int i = 0; i < num; i++)
    {
        xor_mask ^= (1UL << Wz_circuit[i]);
    }

#pragma omp parallel for
    for (int c = 0; c < Hdim; c++)
    {
        int state = c;
        int count = (__builtin_popcount(c & xor_mask)) % 2;
        output[c] += (1 - 2.0 * count) * input[state];
    }
}

void LOBPCG::WzMoments(my_cplx *Phi, vector<vector<int>> Wzs)
{

    for (vector<int> wz : Wzs)
    {
        for (int site : wz)
        {
            if (site < 0 || site >= Nsites)
            {
                cerr << "Function WzMoments, site index out of range: " << site << endl;
                return;
            }
        }
    }

    // 通过波函数计算 <Wz_j> = <Phi|Wz_j|Phi>
    // Phi是波函数，Wzs是 Wz 的索引

    const string WzPath = "Wz_All_" + to_string((int)Sample_num) + ".dat";
    ofstream WzFile(WzPath, ios::app);

    if (!WzFile.is_open())
    {
        cerr << "Function WzMoments, Failed to open: " << WzPath << endl;
        return;
    }

    int Num_Wz = Wzs.size();
    cout << "the num of Wz: " << Num_Wz << endl;
    my_cplx *V = new my_cplx[Hdim];

    WzFile << File_theta << "    " << File_h << "    ";
    for (int j = 0; j < Num_Wz; j++)
    {
        fill(V, V + Hdim, my_cplx(0.0, 0.0));

        int num = Wzs[j].size();
        Wz(num, Wzs[j], Phi, V);                 // W_j|Phi>
        my_cplx z = Vec_Dot(Hdim, Phi, 1, V, 1); // <Phi|W_j|Phi>

#pragma omp critical
        {
            cout << "Wz: " << "(" << setprecision(8) << fixed << right << setw(11) << -z.real() << "," << setw(11) << -z.imag() << ")" << endl;

            WzFile << fixed << setw(24) << right << setprecision(16) << -z.real(); // <Phi|W_j|Phi>
        }
    }
    WzFile << endl;
    delete[] V;
    WzFile.close();

    cout << endl;
}

void LOBPCG::P(int a, int j, int k, my_cplx *input, my_cplx *output)
{
    // 计算 P^a_{j,k} = i/2 * (S^y_j S^z_k - S^z_j S^y_k) 或者 -1/2 * (S^z_j S^x_k -  S^x_j S^z_k) 或者 i/2 * (S^x_j S^y_k - S^y_j S^x_k)
    // a = 0,1,2 分别对应 x,y,z
    // j,k是格点索引，input是波函数，output是输出的结果

    assert(a == 0 || a == 1 || a == 2);

    // 由于P^a_{j,j} 的本质是 描述格点 j和k 之间的自旋交换作用，如果让 j=k，它就变成了单点自旋算符，而不是一个两体相互作用
    if (j == k)
    {
        cout << "Error, P^a_{j,j} = i S^a_j is simply the on-site moment." << endl;
        return;
    }

    if (a == 0) // x  (S^y_j S^z_k - S^z_j S^y_k)
    {
#pragma omp parallel for

        for (int c = 0; c < Hdim; c++)
        {
            // ^:异或：相同为0，不同为1，
            // &：同1为1，否则为0

            // 计算自旋向下(1)的数目
            // down came from up(注意：Sy的符号)
            //**************************************************/
            // 这这里产生一个 -1 是因为这里的 k 没有翻转， 即：down came from down

            int state = c;
            int count = (__builtin_popcount(state & ((1UL << j) ^ (1UL << k)))) % 2; //(S^y_j S^z_k)
            output[c] += 0.25 * (-I) * (1.0 - 2.0 * count) * input[state ^ (1UL << j)];

            count = (__builtin_popcount(state & ((1UL << j) ^ (1UL << k)))) % 2;        // -(S^z_j S^y_k)
            output[c] -= 0.25 * (-I) * (1.0 - 2.0 * count) * input[state ^ (1UL << k)]; // (+i)^1(-1)^1

            /*
            uint32_t state = c;
            uint32_t count = (__builtin_popcount(state & ((1UL << j) ^ (1UL << k)))) % 2;							   // (S^y_j S^z_k - S^z_j S^y_k)/2
            output[c] += 0.125 * (-I) * (1.0 - 2.0 * count) * (input[state ^ (1UL << j)] - input[state ^ (1UL << k)]); // (+i)^1(-1)^1
            */
        }
    }

    if (a == 1) // y  (S^z_j S^x_k - S^x_j S^z_k)
    {
#pragma omp parallel for
        for (int c = 0; c < Hdim; c++)
        {
            int state = c;
            int count = (__builtin_popcount(state & (1UL << j))) % 2;                      // Sz不翻转，只影响符号                  // +S^z_j S^x_k/2
            output[c] += 0.25 * (-1.0) * (-1.0 + 2.0 * count) * input[state ^ (1UL << k)]; // +(+i)^0(-1)^1

            count = (__builtin_popcount(state & (1UL << k))) % 2;                          // -S^x_j S^z_k/2
            output[c] -= 0.25 * (-1.0) * (-1.0 + 2.0 * count) * input[state ^ (1UL << j)]; // -(+i)^0(-1)^1

            /*

            uint32_t state = c;
            uint32_t count = (__builtin_popcount(state & (1UL << j))) % 2;				   // +S^z_j S^x_k/2
            output[c] += 0.125 * (-1.0) * (1.0 - 2.0 * count) * input[state ^ (1UL << k)]; // +(+i)^0(-1)^1
            count = (__builtin_popcount(state & (1UL << k))) % 2;						   // -S^x_j S^z_k/2
            output[c] += 0.125 * (1.0 - 2.0 * count) * input[state ^ (1UL << j)];		   // -(+i)^0(-1)^1

            */
        }
    }

    if (a == 2) // z   (S^x_j S^y_k - S^y_j S^x_k)
    {
#pragma omp parallel for
        for (int c = 0; c < Hdim; c++)
        {
            int state = c;                                                                           // Sx不影响符号，只翻转
            int count = (__builtin_popcount(state & (1UL << k))) % 2;                                // +S^x_j S^y_k/2
            output[c] += 0.25 * I * (-1.0 + 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k))]; // +(+i)^1(-1)^0
            count = (__builtin_popcount(state & (1UL << j))) % 2;                                    // -S^y_j S^x_k/2
            output[c] -= 0.25 * I * (-1.0 + 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k))]; // -(+i)^1(-1)^0

            /*

            uint32_t state = c;
            uint32_t count = (__builtin_popcount(state & (1UL << k))) % 2;								// +S^x_j S^y_k/2
            output[c] += 0.125 * I * (1.0 - 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k))];	// +(+i)^1(-1)^0
            count = (__builtin_popcount(state & (1UL << j))) % 2;										// -S^y_j S^x_k/2
            output[c] += 0.125 * (-I) * (1.0 - 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k))]; // -(+i)^1(-1)^0

            */
        }
    }
}

void LOBPCG::PMoments(my_cplx *input)
{

    const string PMomentPath = "PMoments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";

    ofstream PMomentFile(PMomentPath);

    PMomentFile << left << fixed << setprecision(12) << setfill(' ');
    PMomentFile << setw(4) << "j" << setw(4) << "k" << setw(4) << "a" << setw(32) << "<P^a_{j,k}>" << endl;

    cout << "Calculating <P^a_{j,k}> -> " << PMomentPath << " ... " << flush;

    my_cplx *V = new my_cplx[Hdim];

    TimeVar ti = timeNow();
    for (uint32_t j = 0; j < Nsites; j++)
    {
        for (uint32_t k = j + 1; k < Nsites; k++)
        {
            for (uint32_t a = 0; a < 3; a++)
            {

                fill(V, V + Hdim, complex<double>(0.0, 0.0));
                P(a, j, k, input, V);                              // P^a_{j,k}|Phi>
                complex<double> z = Vec_Dot(Hdim, input, 1, V, 1); // <Phi|P^a_{j,k}|Phi>

#pragma omp critical
                {
                    PMomentFile << left << setw(4) << j << setw(4) << k << setw(4) << a;
                    PMomentFile << right << setw(15) << z.real() << setw(16) << z.imag() << endl; // <Phi|P^a_{j,k}|Phi>
                }
            }
        }
    }
    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    delete[] V;

    PMomentFile.close();

    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::PPCorrelators(my_cplx *input)
{
    // 计算 <Phi|P^a_{j,k} P^b_{l,m}|Phi>
    // Phi是波函数，input是波函数

    // check_bond.dat 第一列是键的类型(0-x, 1-y, 2-z)，第二列是键的两个格点的索引
    const string bondF = "check_bond.dat";
    my_mat bonds;
    Mat_Read(bondF, bonds);
    int NumBond = bonds.size();

    const string PPCorrPath = "PPCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream PPCorrFile(PPCorrPath);

    PPCorrFile << left << fixed << setprecision(12) << setfill(' ');
    PPCorrFile << setw(4) << "j" << setw(4) << "k" << setw(4) << "a" << setw(4) << "l" << setw(4) << "m" << setw(4) << "b" << setw(32) << "<P^a_{j,k} P^b_{l,m}>" << endl;
    cout << "Calculating <P^a_{j,k} P^b_{l,m}> -> " << PPCorrPath << " ... " << flush;

    my_cplx *V1 = new my_cplx[Hdim];
    my_cplx *V2 = new my_cplx[Hdim];

    TimeVar ti = timeNow();
    for (uint32_t b = 0; b < NumBond; b++)
    {
        for (uint32_t bp = 0; bp < NumBond; bp++)
        {
            uint32_t j = bonds[b][1];
            uint32_t k = bonds[b][2];

            uint32_t l = bonds[bp][1];
            uint32_t m = bonds[bp][2];

            for (uint32_t a = 0; a < 3; a++)
            {
                for (uint32_t b = 0; b < 3; b++)
                {

                    fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
                    P(b, l, m, input, V1); // P^b_{l,m}|Phi>

                    fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
                    P(a, j, k, V1, V2);                                 // P^a_{j,k} P^b_{l,m}|Phi>
                    complex<double> z = Vec_Dot(Hdim, input, 1, V2, 1); // <Phi|P^a_{j,k} P^b_{l,m}|Phi>

#pragma omp critical
                    {
                        PPCorrFile << left << setw(4) << j << setw(4) << k << setw(4) << a << setw(4) << l << setw(4) << m << setw(4) << b;
                        PPCorrFile << right << setw(15) << z.real() << setw(16) << z.imag() << endl; // <Phi|P^a_{j,k} P^b_{l,m}|Phi>
                    }
                }
            }
        }
    }

    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    delete[] V1;
    delete[] V2;

    PPCorrFile.close();
    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::Q(int a, int b, int j, int k, my_cplx *input, my_cplx *output)
{
    // 四极矩算符
    // 计算 Q^{a,b}_{j,k} = S^a_j S^b_k
    // a,b = 0,1,2 分别对应 x,y,z
    // j,k是格点索引，input是波函数，output是输出的结果

    if (j == k)
    {
        cout << "Error, Q^{a,b}_{j,j} = 0 for S = 1/2. Offence at j = " << j << "." << endl;
        return;
    }

    if ((a == 0) & (b == 0)) // xx -> (2 xx - yy - zz)/3
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            uint32_t state = c;
            output[c] += (2.0 / 3.0) * 0.25 * input[state ^ ((1UL << j) ^ (1UL << k))]; // 2 S^x_j S^x_k/3

            // down came from up
            uint32_t count = (__builtin_popcount(state & ((1UL << j) ^ (1UL << k)))) % 2;
            output[c] += (-1.0 / 3.0) * 0.25 * (-1.0) * (1.0 - 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k))]; // (+i)^2(-1)^0, -S^y_j S^y_k/3

            output[c] += (-1.0 / 3.0) * 0.25 * (1.0 - 2.0 * count) * input[c]; // (+i)^0(-1)^2, -S^z_j S^z_k/3
        }
    }
    if ((a == 1) & (b == 1)) // yy -> (- xx + 2 yy - zz)/3
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            uint32_t state = c;
            output[c] += (-1.0 / 3.0) * 0.25 * input[state ^ ((1UL << j) ^ (1UL << k))]; // - S^x_j S^x_k/3
            uint32_t count = (__builtin_popcount(state & ((1UL << j) ^ (1UL << k)))) % 2;
            output[c] += (2.0 / 3.0) * 0.25 * (-1.0) * (1.0 - 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k))]; // (+i)^2(-1)^0, 2 S^y_j S^y_k/3
            output[c] += (-1.0 / 3.0) * 0.25 * (1.0 - 2.0 * count) * input[c];                                         // (+i)^0(-1)^2, -S^z_j S^z_k/3
        }
    }
    if ((a == 2) & (b == 2)) // zz -> (-xx - yy + 2 zz)/3
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            uint32_t state = c;
            output[c] += (-1.0 / 3.0) * 0.25 * input[state ^ ((1UL << j) ^ (1UL << k))]; // - S^x_j S^x_k/3
            uint32_t count = (__builtin_popcount(state & ((1UL << j) ^ (1UL << k)))) % 2;
            output[c] += (-1.0 / 3.0) * 0.25 * (-1.0) * (1.0 - 2.0 * count) * input[state ^ ((1UL << j) ^ (1UL << k))]; // (+i)^2(-1)^0, -S^y_j S^y_k/3
            output[c] += (2.0 / 3.0) * 0.25 * (1.0 - 2.0 * count) * input[c];                                           // (+i)^0(-1)^2, 2S^z_j S^z_k/3
        }
    }
    if ((a == 0) & (b == 1)) // xy -> (xy + yx)/2
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            uint32_t state = c;
            // 只确定这一个位置 j 的自旋方向
            // down came from up
            uint32_t count1 = (__builtin_popcount(state & (1UL << k))) % 2;                                   // S^x_j S^y_k/2
            output[c] += 0.50 * 0.250 * I * (-1.0 + 2.0 * count1) * input[state ^ ((1UL << k) ^ (1UL << j))]; // (+i)^1(-1)^0

            uint32_t count2 = (__builtin_popcount(state & (1UL << j))) % 2;                                   // S^y_j S^x_k/2
            output[c] += 0.50 * 0.250 * I * (-1.0 + 2.0 * count2) * input[state ^ ((1UL << k) ^ (1UL << j))]; // (+i)^1(-1)^0
        }
    }
    if ((a == 0) & (b == 2)) // xz -> (xz + zx)/2
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {

            //**************************************************/
            // 这这里产生一个 -1 是因为这里的 k 没有翻转， 即：down came from down
            uint32_t state = c;
            uint32_t count1 = (__builtin_popcount(state & (1UL << k))) % 2;                // S^x_j S^z_k/2
            output[c] += 0.125 * (-1) * (-1.0 + 2.0 * count1) * input[state ^ (1UL << j)]; // (+i)^0(-1)^1

            uint32_t count2 = (__builtin_popcount(state & (1UL << j))) % 2;                // S^z_j S^x_k/2
            output[c] += 0.125 * (-1) * (-1.0 + 2.0 * count2) * input[state ^ (1UL << k)]; // (+i)^0(-1)^1
        }
    }
    if ((a == 1) & (b == 2)) // yz -> (yz + zy)/2
    {
#pragma omp parallel for
        for (uint32_t c = 0; c < Hdim; c++)
        {
            uint32_t state = c;
            uint32_t count = (__builtin_popcount(state & ((1UL << j) ^ (1UL << k)))) % 2;  // (S^y_j S^z_k + S^z_j S^y_k)/2
            output[c] += 0.125 * (-I) * (1.0 - 2.0 * count) * (input[state ^ (1UL << j)]); // (+i)^1(-1)^1
            output[c] += 0.125 * (-I) * (1.0 - 2.0 * count) * (input[state ^ (1UL << k)]);
        }
    }

    if (a > b)
    {
        cout << "switching" << endl;
        Q(b, a, j, k, input, output);
    }
}

void LOBPCG::QMoments(my_cplx *input)
{
    // 计算 <Q^{a,b}_{j,k}> = <Phi|Q^{a,b}_{j,k}|Phi>
    // Phi是波函数，input是波函数

    const string QMomentsPath = "QMoments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream QMomentsFile(QMomentsPath);

    QMomentsFile << left << fixed << setprecision(12) << setfill(' ');
    QMomentsFile << setw(4) << "j" << setw(4) << "k" << setw(4) << "a" << setw(4) << "b" << setw(32) << "<Q^{a,b}_{j,k}>" << endl;
    cout << "Calculating <Q^{a,b}_{j,k}> -> " << QMomentsPath << " ... " << flush;

    my_cplx *V = new my_cplx[Hdim];

    TimeVar ti = timeNow();
    for (uint32_t j = 0; j < Nsites; j++)
    {
        for (uint32_t k = j + 1; k < Nsites; k++)
        {
            for (uint32_t a = 0; a < 3; a++)
            {
                for (uint32_t b = a; b < 3; b++)
                {

                    fill(V, V + Hdim, complex<double>(0.0, 0.0));
                    Q(a, b, j, k, input, V);                           // Q^{a,b}_{j,k}|Phi>
                    complex<double> z = Vec_Dot(Hdim, input, 1, V, 1); // <Phi|Q^{a,b}_{j,k}|Phi>

#pragma omp critical
                    {
                        QMomentsFile << left << setw(4) << j << setw(4) << k << setw(4) << a << setw(4) << b;
                        QMomentsFile << right << setw(15) << z.real() << setw(16) << z.imag() << endl; // <Phi|Q^{a,b}_{j,k}|Phi>
                    }
                }
            }
        }
    }
    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    delete[] V;

    QMomentsFile.close();
    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::QQCorrelators(my_cplx *input)
{
    // check_bond.dat 第一列是键的类型(0-x, 1-y, 2-z)，第二列是键的两个格点的索引
    const string bondF = "check_bond.dat";
    my_mat bonds;
    Mat_Read(bondF, bonds);
    int NumBond = bonds.size();

    const string QQMomentsPath = "QQCorr_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream QQMomentsFile(QQMomentsPath);

    QQMomentsFile << left << fixed << setprecision(12) << setfill(' ');
    QQMomentsFile << setw(4) << "j" << setw(4) << "k" << setw(4) << "a" << setw(4) << "b" << setw(4) << "l" << setw(4) << "m" << setw(4) << "c" << setw(4) << "d";
    QQMomentsFile << setw(32) << "<Q^{a,b}_{j,k} Q^{c,d}_{l,m}>" << endl;
    cout << "Calculating <Q^{a,b}_{j,k} Q^{c,d}_{l,m}> -> " << QQMomentsPath << " ... " << flush;

    my_cplx *V1 = new my_cplx[Hdim];
    my_cplx *V2 = new my_cplx[Hdim];

    TimeVar ti = timeNow();
    for (uint32_t b = 0; b < NumBond; b++)
    {
        uint32_t j = bonds[b][1];
        uint32_t k = bonds[b][2];
        for (uint32_t bp = 0; bp < NumBond; bp++)
        {
            uint32_t l = bonds[bp][1];
            uint32_t m = bonds[bp][2];
            for (uint32_t a = 0; a < 3; a++)
            {
                for (uint32_t b = a; b < 3; b++)
                {
                    for (uint32_t c = 0; c < 3; c++)
                    {
                        for (uint32_t d = c; d < 3; d++)
                        {

                            fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
                            Q(c, d, l, m, input, V1); // Q^{c,d}_{l,m}|Phi>

                            fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
                            Q(a, b, j, k, V1, V2);                              // Q^{a,b}_{j,k} Q^{c,d}_{l,m}|Phi>
                            complex<double> z = Vec_Dot(Hdim, input, 1, V2, 1); // <Phi|Q^{a,b}_{j,k} Q^{c,d}_{l,m}|Phi>

#pragma omp critical
                            {
                                QQMomentsFile << left << setw(4) << j << setw(4) << k << setw(4) << a << setw(4) << b;
                                QQMomentsFile << left << setw(4) << l << setw(4) << m << setw(4) << c << setw(4) << d;
                                QQMomentsFile << right << setw(15) << z.real() << setw(16) << z.imag() << endl; // <Phi|Q^{a,b}_{j,k} Q^{c,d}_{l,m}|Phi>
                            }
                        }
                    }
                }
            }
        }
    }
    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    delete[] V1;
    delete[] V2;

    QQMomentsFile.close();
    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}

void LOBPCG::SSSS_Correlators(my_cplx *input)
{
    // check_bond.dat 第一列是键的类型(0-x, 1-y, 2-z)，第二列是键的两个格点的索引
    const string bondF = "check_bond.dat";
    my_mat bonds;
    Mat_Read(bondF, bonds);
    int NumBond = bonds.size();

    const string SSSSMomentsPath = "SSSSMoments_" + File_theta + "_" + File_h + "_" + to_string((int)Sample_num) + ".dat";
    ofstream SSSSMomentsFile(SSSSMomentsPath);

    SSSSMomentsFile << left << fixed << setprecision(12) << setfill(' ');
    SSSSMomentsFile << setw(4) << "j" << setw(4) << "a" << setw(4) << "k" << setw(4) << "b" << setw(4) << "l" << setw(4) << "c" << setw(4) << "m" << setw(4) << "d";
    SSSSMomentsFile << setw(32) << "<S^a_j S^b_k S^c_l S^d_m>" << endl;
    cout << "Calculating <S^a_j S^b_k S^c_l S^d_m> -> " << SSSSMomentsPath << " ... " << flush;

    my_cplx *V1 = new my_cplx[Hdim];
    my_cplx *V2 = new my_cplx[Hdim];
    my_cplx *V3 = new my_cplx[Hdim];
    my_cplx *V4 = new my_cplx[Hdim];

    TimeVar ti = timeNow();
    for (uint32_t b = 0; b < NumBond; b++)
    {
        uint32_t j = bonds[b][1];
        uint32_t k = bonds[b][2];
        for (uint32_t bp = 0; bp < NumBond; bp++)
        {
            uint32_t l = bonds[bp][1];
            uint32_t m = bonds[bp][2];
            for (uint32_t a = 0; a < 3; a++)
            {
                for (uint32_t b = 0; b < 3; b++)
                {
                    for (uint32_t c = 0; c < 3; c++)
                    {
                        for (uint32_t d = 0; d < 3; d++)
                        {

                            fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
                            fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
                            fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));
                            fill(V4, V4 + Hdim, my_cplx(0.0, 0.0));

                            S(d, m, input, V1); // S^d_m|Phi>
                            S(c, l, V1, V2);    // S^c_l S^d_m|Phi>
                            S(b, k, V2, V3);    // S^b_k S^c_l S^d_m|Phi>
                            S(a, j, V3, V4);    // S^a_j S^b_k S^c_l S^d_m|Phi>

                            complex<double> z = Vec_Dot(Hdim, input, 1, V4, 1); // <Phi|S^a_j S^b_k S^c_l S^d_m|Phi>

#pragma omp critical
                            {
                                SSSSMomentsFile << left << setw(4) << j << setw(4) << a << setw(4) << k << setw(4) << b;
                                SSSSMomentsFile << left << setw(4) << l << setw(4) << c << setw(4) << m << setw(4) << d;
                                SSSSMomentsFile << right << setw(15) << z.real() << setw(16) << z.imag() << endl; // <Phi|Q^{a,b}_{j,k} Q^{c,d}_{l,m}|Phi>
                            }
                        }
                    }
                }
            }
        }
    }
    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    delete[] V1;
    delete[] V2;
    delete[] V3;
    delete[] V4;

    SSSSMomentsFile.close();
    cout << "done. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
    cout << endl;
}
