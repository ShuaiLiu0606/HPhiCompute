#ifndef LOBPCG_H
#define LOBPCG_H

#include <common.h>
#include <LOBPCG.h>

class LOBPCG
{
private:
    int Nsites; // 格点数

    double Spin;
    int Sdim; // 自旋自由度(2*SpinVal + 1)

    uint64_t Hdim; // Hamiltonian的维度

    int Sample_num; // 控制基态还是激发态,0为基态，1为第一激发态

    const string File_theta; // theta的大小，生成文件时的后缀

    const string File_h; // 磁场的大小，生成文件时的后缀

public:
    // 定义 Sx, Sy, Sz
    vector<pair<string, int>> SxBose = {{"01", 1}, {"10", 1}};
    vector<pair<string, int>> SyBose = {{"01", 1}, {"10", -1}};
    vector<pair<string, int>> SzBose = {{"00", 1}, {"11", -1}};

    LOBPCG(const string &ParaFilePath, const string file_theta, const string file_h, int sample_num);

    void Fidelity(my_cplx_mat Phi_Mat, int NumTheta, int Numfile); // 计算保真度

    void EnergySpectrum(const string &EnPath); // 计算能谱

    void Phi_Read(const string &FilePath, my_cplx_vec &EigenVec); // 读取波函数

    void Halfchain_EE(my_cplx *Phi); // 计算半链纠缠熵

    double EE(my_cplx *Phi, const uint32_t &size_A, uint32_t *sites_A); // 计算纠缠熵

    void TEE(my_cplx *Phi, const uint32_t &size_A, uint32_t *sites_A, const uint32_t &size_B, uint32_t *sites_B, const uint32_t &size_C, uint32_t *sites_C); // 计算拓扑纠缠熵

    void S_Corr(const string OneBodyPath); // 通过单体格林函数文件计算单体算符的关联函数<S^x><S^y><S^z>

    void SS_PCorr(const string TwoBodyPath);
    void SS_QCorr(const string TwoBodyPath);
    void SS_Corr(const string TwoBodyPath); // 通过双体格林函数文件计算双体算符的关联函数<S^xS^x><S^yS^y><S^zS^z>....

    void SSS_Corr(const string ThreeBodyPath, int numh); // 通过三体格林函数文件计算三体算符的关联函数<S_i^a S_j^b S_k^c>....

    void SSSS_Corr(const string FourBodyPath);
    void SSSS_PPCorr(const string FourBodyPath);
    void SSSS_QQCorr(const string FourBodyPath);

    void SSSSSS_WpCorr(const string SixBodyPath);

    //*************************************************

    void S(int a, int j, my_cplx *input, my_cplx *output); // 通过波函数计算单体算符S^a_j|phi>

    void SpSm(int a, int j, my_cplx *input, my_cplx *output); // 通过波函数计算单体算符S^+_j|phi>和S^-_j|phi>

    void Sq(int a, int j, my_cplx *input, my_cplx *output, vector<vector<double>> R, pair<double, double> q); // 通过波函数计算单体算符S^a_q|phi> S^a_q = \sum_{j} e^{i*q*R_j} S^a_j

    void S_Moments(my_cplx *Phi); // 通过波函数计算单体算符的关联函数<S^x><S^y><S^z>

    void SS_Correlators(my_cplx *Phi); // 通过波函数计算双体算符的关联函数<S^xS^x><S^yS^y><S^zS^z>

    void P(int a, int j, int k, my_cplx *input, my_cplx *output); // 通过波函数计算P算符 output = P * |input> P^a_{j,k} = i/2 * (S^y_j S^z_k - S^z_j S^y_k) 或者 -1/2 * (S^z_j S^x_k -  S^x_j S^z_k) 或者 i/2 * (S^x_j S^y_k - S^y_j S^x_k)
    void PMoments(my_cplx *input);                                // 通过波函数计算P算符的关联函数<P>
    void PPCorrelators(my_cplx *input);                           // 通过波函数计算PP算符的关联函数<PP>

    void Q(int a, int b, int j, int k, my_cplx *input, my_cplx *output); // 通过波函数计算Q算符 output = Q * |input> Q^a_b_{j,k} = i/2 * (S^a_j S^b_k - S^b_j S^a_k) 或者 -1/2 * (S^b_j S^a_k -  S^a_j S^b_k) 或者 i/2 * (S^b_j S^a_k - S^a_j S^b_k)
    void QMoments(my_cplx *input);                                       // 通过波函数计算Q算符的关联函数<Q>
    void QQCorrelators(my_cplx *input);                                  // 通过波函数计算QQ算符的关联函数<QQ>

    void SSSS_Correlators(my_cplx *input); // 通过波函数计算四体算符的关联函数<S_i^a S_j^b S_k^c S_l^d>....

    void ChiralityMoments(my_cplx *Phi); // 通过波函数计算三体算符的关联函数<S_i^a S_j^b S_k^c>....

    void Wp(vector<int> Plaquette, my_cplx *input, my_cplx *output); // 通过波函数计算环算符Wp output = Wp * |input>
    void WpMoments(my_cplx *Phi); // 通过波函数计算环算符的关联函数<Wp>
    void WpWpCorrelators(my_cplx *input); // 通过波函数计算环算符的关联函数<WpWp>

    void Wx(int num, vector<int> Wx_circuit, my_cplx *input, my_cplx *output); // 通过波函数计算Wx算符 output = Wx * |input>
    void WxMoments(my_cplx *Phi, vector<vector<int>> Wxs); // 通过波函数计算Wx算符的关联函数<Wx>

    void Wy(int num, vector<int> Wy_circuit, my_cplx *input, my_cplx *output); // 通过波函数计算Wy算符 output = Wy * |input>
    void WyMoments(my_cplx *Phi, vector<vector<int>> Wys); // 通过波函数计算Wy算符的关联函数<Wy>

    void Wz(int num, vector<int> Wz_circuit, my_cplx *input, my_cplx *output); // 通过波函数计算Wz算符 output = Wz * |input>
    void WzMoments(my_cplx *Phi, vector<vector<int>> Wzs); // 通过波函数计算Wz算符的关联函数<Wz>

 
    void Magnetization(const string OneBodyPath, int num_h); // 计算磁化强度

    void Binder(my_cplx *input); // 计算BinderCumulant

};

#endif // LOBPCG_HI
