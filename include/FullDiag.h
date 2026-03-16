#ifndef FullDiag_H
#define FullDiag_H

#include <common.h>
#include <FullDiag.h>

class FullDiag
{
private:
    int Nsites; // 格点数
    double Spin;
    int Sdim; // 自旋自由度(2*SpinVal + 1)
    int Hdim; // Hamiltonian的维度

    const string File_theta; // theta的大小，生成文件时的后缀
    const string File_h;     // 磁场的大小，生成文件时的后缀

public:
    FullDiag(const string &ParaFilePath, const string file_theta, const string file_h); // 构造函数，从文件中读取参数并赋值给成员变量

    void Fulldiag_cal(const string &in_file, double T); // 计算物理量

    void WpMoments_Fulldiag(const string &input_file, double *Wps); // 计算Wp矩

    void S_CorrFullDiag(const string OneBodyPath, double *SxCorr, double *SyCorr, double *SzCorr); // 计算 <S> 自旋相关函数

    void SS_CorrFullDiag(const string TwoBodyPath, double *SxSxCorr, double *SySyCorr, double *SzSzCorr); // 计算 <SS> 自旋相关函数

    void Read_FullDiag_Phys(const string &file_path, vector<double> &Sz, vector<double> &S2); // 读取完全对角化的物理量
};

#endif // FullDiag_H