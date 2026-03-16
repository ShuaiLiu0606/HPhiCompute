#include <common.h>
#include <LOBPCG.h>

void LOBPCG::Magnetization(const string OneBodyPath, int num_h)
{

    my_cplx_mat SxCorr, SyCorr, SzCorr;

    const string SxCorrPath = "SxCorr_All_" + to_string((int)Sample_num) + ".dat";
    const string SyCorrPath = "SyCorr_All_" + to_string((int)Sample_num) + ".dat";
    const string SzCorrPath = "SzCorr_All_" + to_string((int)Sample_num) + ".dat";

    Mat_Read(SxCorrPath, SxCorr);
    Mat_Read(SyCorrPath, SyCorr);
    Mat_Read(SzCorrPath, SzCorr);

    // 铁磁序参量
    complex<double> TotalSx(0.0, 0.0);
    complex<double> TotalSy(0.0, 0.0);
    complex<double> TotalSz(0.0, 0.0);

    for (int i = 0; i < Nsites; i++)
    {
        TotalSx += SxCorr[num_h][i];
        TotalSy += SyCorr[num_h][i];
        TotalSz += SzCorr[num_h][i];
    }

    // 计算磁化强度
    const string MgaPath = "Magnetization_All_" + to_string((int)Sample_num) + ".dat";
    complex<double> Magnetization = (TotalSx + TotalSy + TotalSz) / my_cplx(Nsites, 0.0) / my_cplx(sqrt(3.0), 0.0);

    ofstream MagFile(MgaPath, ios::app);

    MagFile << File_theta << "    " << File_h << "    ";
    MagFile << fixed << setw(21) << right << setprecision(16) << TotalSx.real();
    MagFile << fixed << setw(21) << right << setprecision(16) << TotalSy.real();
    MagFile << fixed << setw(21) << right << setprecision(16) << TotalSz.real();
    MagFile << fixed << setw(21) << right << setprecision(16) << Magnetization.real() << endl;

    MagFile.close();

    //**************************************************************transform */

    // a方向[1 1 -2]; b方向[-1 1 0]; c方向[1 1 1]
    vector<double> a = {1 / sqrt(6), 1 / sqrt(6), -2 / sqrt(6)}; // A = 1 / sqrt(6) * (1, 1, -2)
    vector<double> b = {-1 / sqrt(2), 1 / sqrt(2), 0};           // B = 1 / sqrt(2) * (-1, 1, 0)
    vector<double> c = {1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)};  // C = 1 / sqrt(3) * (1, 1, 1)

    // Write the transformed matrices to files
    const string SaCorrPath = "SaCorr_All_" + to_string((int)Sample_num) + ".dat";
    const string SbCorrPath = "SbCorr_All_" + to_string((int)Sample_num) + ".dat";
    const string ScCorrPath = "ScCorr_All_" + to_string((int)Sample_num) + ".dat";

    ofstream SaCorrFile(SaCorrPath, ios::app);
    ofstream SbCorrFile(SbCorrPath, ios::app);
    ofstream ScCorrFile(ScCorrPath, ios::app);

    // Perform the transformation
    for (size_t i = 0; i < Nsites; i++)
    {

        my_cplx SaCorr = a[0] * SxCorr[num_h][i] + a[1] * SyCorr[num_h][i] + a[2] * SzCorr[num_h][i];
        my_cplx SbCorr = b[0] * SxCorr[num_h][i] + b[1] * SyCorr[num_h][i] + b[2] * SzCorr[num_h][i];
        my_cplx ScCorr = c[0] * SxCorr[num_h][i] + c[1] * SyCorr[num_h][i] + c[2] * SzCorr[num_h][i];

        SaCorrFile << fixed << setprecision(16) << right << SaCorr << "    ";
        SbCorrFile << fixed << setprecision(16) << right << SbCorr << "    ";
        ScCorrFile << fixed << setprecision(16) << right << ScCorr << "    ";
    }
    SaCorrFile << endl;
    SbCorrFile << endl;
    ScCorrFile << endl;

    SaCorrFile.close();
    SbCorrFile.close();
    ScCorrFile.close();
}


void LOBPCG::Binder(my_cplx *input)
{
    // 计算BinderCumulant
    // 公式: U = 1 - <M^4> / 3<M^2>^2

    double Mz2 = 0.0;
    double Mz4 = 0.0;

    //*****************************************
    // 情况1：四个点不同 C_{4}^{N}
    //*****************************************
    double result4 = 0.0;
    TimeVar ti0_1 = timeNow();
    cout << "Four different points:" << endl;
    for (int i = 0; i < Nsites - 3; i++)
    {
        for (int j = i + 1; j < Nsites - 2; j++)
        {
            for (int k = j + 1; k < Nsites - 1; k++)
            {
                for (int l = k + 1; l < Nsites; l++)
                {
                    my_cplx *V1 = new my_cplx[Hdim];
                    my_cplx *V2 = new my_cplx[Hdim];
                    my_cplx *V3 = new my_cplx[Hdim];
                    my_cplx *V4 = new my_cplx[Hdim];

                    fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
                    fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
                    fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));
                    fill(V4, V4 + Hdim, my_cplx(0.0, 0.0));

                    S(2, i, input, V1);
                    S(2, j, V1, V2);
                    S(2, k, V2, V3);
                    S(2, l, V3, V4);

                    my_cplx result = Vec_Dot(Hdim, input, 1, V4, 1);
                    result4 += result.real();

                    delete[] V1;
                    delete[] V2;
                    delete[] V3;
                    delete[] V4;
                }
            }
        }
    }
    result4 *= 24;

    TimeVar ti0_2 = timeNow();
    auto diff1 = durationms(ti0_2 - ti0_1) / 1000.0;
    cout << "Four different points done. Took " << left << fixed << setprecision(2) << diff1 << " s." << endl
         << endl;

    //*****************************************
    // 情况2：两个点相同，两个点不同 C_{2}^{N} * C_{2}^{N-2}
    //*****************************************

    double result3_2 = 0.0;
    TimeVar ti1_1 = timeNow();
    cout << "Two points the same, two different:" << endl;
    for (int i = 0; i < Nsites; i++)
    {
        for (int j = 0; j < Nsites - 1; j++)
        {
            if (j == i)
            {
                continue; // 保证j!=i
            }

            for (int k = j + 1; k < Nsites; k++)
            {
                if (k == i)
                {
                    continue; // 保证k!=i
                }

                my_cplx *V1 = new my_cplx[Hdim];
                my_cplx *V2 = new my_cplx[Hdim];
                my_cplx *V3 = new my_cplx[Hdim];
                my_cplx *V4 = new my_cplx[Hdim];

                fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
                fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
                fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));
                fill(V4, V4 + Hdim, my_cplx(0.0, 0.0));

                S(2, i, input, V1);
                S(2, i, V1, V2);
                S(2, j, V2, V3);
                S(2, k, V3, V4);

                my_cplx result = Vec_Dot(Hdim, input, 1, V4, 1);
                result3_2 += result.real();

                delete[] V1;
                delete[] V2;
                delete[] V3;
                delete[] V4;
            }
        }
    }
    result3_2 *= 12;

    TimeVar ti1_2 = timeNow();
    auto diff2 = durationms(ti1_2 - ti1_1) / 1000.0;
    cout << "Two points the same, two different done. Took " << left << fixed << setprecision(2) << diff2 << " s." << endl
         << endl;

    //*****************************************
    // 情况3：两对点相同 C_{2}^{N} * C_{2}^{N-2}
    //*****************************************

    double result2_2 = 0.0;
    TimeVar ti2_1 = timeNow();
    cout << "Two pairs of the same point:" << endl;
    for (int i = 0; i < Nsites - 1; i++)
    {
        for (int j = i + 1; j < Nsites; j++)
        {
            my_cplx *V1 = new my_cplx[Hdim];
            my_cplx *V2 = new my_cplx[Hdim];
            my_cplx *V3 = new my_cplx[Hdim];
            my_cplx *V4 = new my_cplx[Hdim];

            fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
            fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
            fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));
            fill(V4, V4 + Hdim, my_cplx(0.0, 0.0));

            S(2, i, input, V1);
            S(2, i, V1, V2);
            S(2, j, V2, V3);
            S(2, j, V3, V4);

            my_cplx result = Vec_Dot(Hdim, input, 1, V4, 1);
            result2_2 += result.real();

            delete[] V1;
            delete[] V2;
            delete[] V3;
            delete[] V4;
        }
    }
    result2_2 *= 6;

    TimeVar ti2_2 = timeNow();
    auto diff3 = durationms(ti2_2 - ti2_1) / 1000.0;
    cout << "Two pairs of the same point done. Took " << left << fixed << setprecision(2) << diff3 << " s." << endl
         << endl;

    //*****************************************
    // 情况4：三个点相同，一个点不同 C_{3}^{N} * C_{1}^{N-3}
    //*****************************************

    double result3_1 = 0.0;
    TimeVar ti3_1 = timeNow();
    cout << "Three points the same, one different:" << endl;
    for (int i = 0; i < Nsites; i++)
    {
        for (int j = 0; j < Nsites; j++)
        {
            if (j == i)
            {
                continue;
            }

            my_cplx *V1 = new my_cplx[Hdim];
            my_cplx *V2 = new my_cplx[Hdim];
            my_cplx *V3 = new my_cplx[Hdim];
            my_cplx *V4 = new my_cplx[Hdim];

            fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
            fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
            fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));
            fill(V4, V4 + Hdim, my_cplx(0.0, 0.0));

            S(2, i, input, V1);
            S(2, i, V1, V2);
            S(2, i, V2, V3);
            S(2, j, V3, V4);

            my_cplx result = Vec_Dot(Hdim, input, 1, V4, 1);
            result3_1 += result.real();

            delete[] V1;
            delete[] V2;
            delete[] V3;
            delete[] V4;
        }
    }
    result3_1 *= 4;

    TimeVar ti3_2 = timeNow();
    auto diff4 = durationms(ti3_2 - ti3_1) / 1000.0;
    cout << "Three points the same, one different done. Took " << left << fixed << setprecision(2) << diff4 << " s." << endl
         << endl;

    //***********************************************
    // 情况5：四个点相同 C_{4}^{N}
    //***********************************************
    double result4_0 = 0.0;
    TimeVar ti4_1 = timeNow();
    cout << "Four points the same:" << endl;
    for (int j = 0; j < Nsites; j++)
    {
        my_cplx *V1 = new my_cplx[Hdim];
        my_cplx *V2 = new my_cplx[Hdim];
        my_cplx *V3 = new my_cplx[Hdim];
        my_cplx *V4 = new my_cplx[Hdim];

        fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
        fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));
        fill(V3, V3 + Hdim, my_cplx(0.0, 0.0));
        fill(V4, V4 + Hdim, my_cplx(0.0, 0.0));

        S(2, j, input, V1);
        S(2, j, V1, V2);
        S(2, j, V2, V3);
        S(2, j, V3, V4);

        my_cplx result = Vec_Dot(Hdim, input, 1, V4, 1);

        result4_0 += result.real();

        delete[] V1;
        delete[] V2;
        delete[] V3;
        delete[] V4;
    }

    TimeVar ti4_2 = timeNow();
    auto diff5 = durationms(ti4_2 - ti4_1) / 1000.0;
    cout << "Four points the same done. Took " << left << fixed << setprecision(2) << diff5 << " s." << endl
         << endl;

    Mz4 = result4 + result3_2 + result2_2 + result3_1 + result4_0;
    Mz4 /= pow(Nsites, 4);

    for (int i1 = 0; i1 < Nsites; i1++)
    {
        for (int i2 = 0; i2 < Nsites; i2++)
        {
            my_cplx *V1 = new my_cplx[Hdim];
            my_cplx *V2 = new my_cplx[Hdim];
            fill(V1, V1 + Hdim, my_cplx(0.0, 0.0));
            fill(V2, V2 + Hdim, my_cplx(0.0, 0.0));

            S(2, i1, input, V1);
            S(2, i2, V1, V2);

            my_cplx result = Vec_Dot(Hdim, input, 1, V2, 1);

            Mz2 += result.real();

            delete[] V1;
            delete[] V2;
        }
    }

    Mz2 /= pow(Nsites, 2);

    const string binder = "BinderCumulant_All_" + to_string((int)Sample_num) + ".dat";

    ofstream Binde_F(binder, ios::app);

    double binder_cumulant = 1 - (Mz4 / (3.0 * pow(Mz2, 2)));

    cout << "Binder Cumulant is: " << setw(11) << fixed << right << setprecision(8) << binder_cumulant << endl;

    Binde_F << File_theta << "    " << File_h << "    ";
    Binde_F << fixed << setprecision(16) << setw(22) << right << binder_cumulant << endl;

    Binde_F.close();

    TimeVar ti5 = timeNow();
    auto diff = durationms(ti5 - ti0_1) / 1000.0;
    cout << "Binder cumulant done. Took " << left << fixed << setprecision(2) << diff << " s." << endl
         << endl;
}


