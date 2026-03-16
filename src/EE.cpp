#include <common.h>
#include <LOBPCG.h>

void LOBPCG::Halfchain_EE(my_cplx *Phi)
{
    // 计算半链纠缠熵

    uint32_t size_A = Nsites / 2;
    uint32_t size_B = Nsites - size_A;

    uint32_t DimA = 1 << size_A;
    uint32_t DimB = 1 << size_B;
    uint32_t dim_2d = DimA * DimB;
    my_cplx *wf_2d = new my_cplx[dim_2d];

    uint32_t *State = new uint32_t[Hdim];
    for (uint32_t s = 0; s < Hdim; s++)
    {
        State[s] = s;
    }

    for (uint32_t s = 0; s < dim_2d; s++)
    {
        wf_2d[s] = 0;
    }
    for (uint32_t s = 0; s < Hdim; s++)
    {
        wf_2d[State[s]] = Phi[s];
    }
    delete[] State;

    double *lambda = new double[min(DimA, DimB)];
    Mat_Svd(LAPACK_ROW_MAJOR, 'N', 'N', DimB, DimA, wf_2d, DimA, lambda, NULL, DimA, NULL, DimB);

    // singular value is the sqrt of eigenvalues of reduced density martrix
    vdSqr(min(DimA, DimB), lambda, lambda);
    double ee = VN_Entropy(min(DimA, DimB), lambda);
    delete[] wf_2d;
    delete[] lambda;

    cout << "EE_Half: " << setw(11) << fixed << right << setprecision(8) << ee << endl;

    const string EE_FileName = "EE_Half_All_" + to_string((int)Sample_num) + ".dat";
    ofstream EE_File(EE_FileName, ios::app);

    EE_File << File_theta << "    " << File_h << "    ";
    EE_File << setw(22) << fixed << right << setprecision(16) << ee << endl;

    EE_File.close();
}

double LOBPCG::EE(my_cplx *Phi, const uint32_t &size_A, uint32_t *sites_A)
{
    // 计算任意格点区域的纠缠熵
    // size_A: 区域A的大小
    // sites_A: 区域A的格点

    uint32_t size_B = Nsites - size_A; // 区域B的大小
    uint32_t DimA = pow(Sdim, size_A); // 列的维度
    uint32_t DimB = pow(Sdim, size_B); // 行的维度
    uint32_t dim_2d = DimA * DimB;
    my_cplx *wf_2d = new my_cplx[dim_2d];

#pragma omp parallel for
    for (uint32_t k = 0; k < Hdim; k++)
    {
        uint32_t state = k;
        uint32_t ind_a, ind_b;
        Bits_Decomposition(state, Nsites, size_A, sites_A, ind_a, ind_b); // bit分解将state分解为两部分
        wf_2d[ind_b * DimA + ind_a] = Phi[k];
    }

    // lambda: 存储奇异值，大小为 min(m, n)。
    double *lambda = new double[min(DimA, DimB)];
    Mat_Svd(LAPACK_ROW_MAJOR, 'N', 'N', DimB, DimA, wf_2d, DimA, lambda, NULL, DimA, NULL, DimB);

    // singular value is the sqrt of eigenvalues of reduced density martrix
    vdSqr(min(DimA, DimB), lambda, lambda);
    double ee = VN_Entropy(min(DimA, DimB), lambda);
    delete[] wf_2d;
    delete[] lambda;

    cout << " EE: " << setw(11) << fixed << right << setprecision(8) << ee << endl;

    return ee;
}

void LOBPCG::TEE(my_cplx *Phi, const uint32_t &size_A, uint32_t *sites_A, const uint32_t &size_B, uint32_t *sites_B, const uint32_t &size_C, uint32_t *sites_C)
{
    cout << "Calculating <TEE> ... " << flush << endl;
    TimeVar ti = timeNow();

    const uint32_t size_AB = size_A + size_B;
    uint32_t *sites_AB = new uint32_t[size_AB];
    for (uint32_t i = 0; i < size_A; i++)
    {
        sites_AB[i] = sites_A[i];
    }
    for (uint32_t j = 0; j < size_B; j++)
    {
        sites_AB[size_A + j] = sites_B[j];
    }

    const uint32_t size_AC = size_A + size_C;
    uint32_t *sites_AC = new uint32_t[size_AC];
    for (uint32_t i = 0; i < size_A; i++)
    {
        sites_AC[i] = sites_A[i];
    }
    for (uint32_t j = 0; j < size_C; j++)
    {
        sites_AC[size_A + j] = sites_C[j];
    }

    const uint32_t size_BC = size_B + size_C;
    uint32_t *sites_BC = new uint32_t[size_BC];
    for (uint32_t i = 0; i < size_B; i++)
    {
        sites_BC[i] = sites_B[i];
    }
    for (uint32_t j = 0; j < size_C; j++)
    {
        sites_BC[size_B + j] = sites_C[j];
    }

    const uint32_t size_ABC = size_A + size_B + size_C;
    uint32_t *sites_ABC = new uint32_t[size_ABC];
    for (uint32_t i = 0; i < size_A; i++)
    {
        sites_ABC[i] = sites_A[i];
    }
    for (uint32_t j = 0; j < size_B; j++)
    {
        sites_ABC[size_A + j] = sites_B[j];
    }
    for (uint32_t k = 0; k < size_C; k++)
    {
        sites_ABC[size_A + size_B + k] = sites_C[k];
    }

    double EE_A = EE(Phi, size_A, sites_A);
    double EE_B = EE(Phi, size_B, sites_B);
    double EE_C = EE(Phi, size_C, sites_C);
    double EE_AB = EE(Phi, size_AB, sites_AB);
    double EE_AC = EE(Phi, size_AC, sites_AC);
    double EE_BC = EE(Phi, size_BC, sites_BC);
    double EE_ABC = EE(Phi, size_ABC, sites_ABC);

    delete[] sites_AB;
    delete[] sites_AC;
    delete[] sites_BC;
    delete[] sites_ABC;

    double TEE = EE_A + EE_B + EE_C - EE_AB - EE_AC - EE_BC + EE_ABC;

    cout << "TEE:" << fixed << setprecision(8) << setw(12) << TEE << endl;

    // 保存TEE到文件
    const string TEE_FileName = "TEE_All_" + to_string((int)Sample_num) + ".dat";
    ofstream TEE_file(TEE_FileName, ios::app);

    TEE_file << File_theta << "    " << File_h << "    ";
    TEE_file << fixed << setw(24) << right << setprecision(16) << TEE << endl;
    TEE_file.close();

    TimeVar tf = timeNow();
    auto diff = durationms(tf - ti) / 1000.0;

    cout << "Calculate TEE down. Took " << left << fixed << setprecision(2) << diff << " s." << endl;
}
