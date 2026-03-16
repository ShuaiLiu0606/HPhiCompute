#ifndef TPQ_H
#define TPQ_H

#include <common.h>
#include <TPQ.h>
#include <LOBPCG.h>

class TPQ
{
private:
    int Nsites; // 格点数
    double Spin;
    int Sdim;          // 自旋自由度(2*SpinVal + 1)
    int Hdim;          // Hamiltonian的维度
    int NumAve;        // TPQ独立运行的次数
    int Lanczos_Max;   // TPQ每次运行的最大迭代次数
    int ExpecInterval; // 输出关联函数的步长
    double LargeValue; // 虚时演化的最大值

    int Sample_num; // 控制基态还是激发态,0为基态，1为第一激发态

    const string File_theta; // theta的大小，生成文件时的后缀
    const string File_h;     // 磁场的大小，生成文件时的后缀

public:
    TPQ(const string &ParaFilePath, const string file_theta, const string file_h, int sample_num); // 构造函数，从文件中读取参数并赋值给成员变量

    void Read_TPQ_file(const string &pri_file, my_mat &Norm, my_mat &InvTemp, my_mat &Ene, my_mat &Ene2, my_mat &Spc, my_mat &Sz, my_mat &Sz2, my_mat &chi); // 读取TPQ文件

    void mainBasic(const string &TPQ_FilePath); // 主函数，执行基本操作

    void mainPhys(const string &TPQ_FilePath); // 主函数，计算物理量

    void CalcBasic(const my_mat &Norm, const my_mat &Ene, const my_mat &Ene2, const my_mat &Sz, const my_mat &Sz2, const my_mat &InvTemp, my_vec &log_Z);

    void BS_Basic(int &max_BS, const my_vec &log_Z, const my_mat &phys_Z, const my_mat &phys_Ene, const my_mat &phys_Ene2, const my_mat &phys_Sz, const my_mat &phys_Sz2, const my_mat &phys_InvTemp);

    void WpMoments_TPQ(const string &SixBodyPath, my_mat &Wp_mat);                                     // 计算Wp矩
    void S_Corr_TPQ(const string &OneBodyPath, my_mat &SxCorr, my_mat &SyCorr, my_mat &SzCorr);        // 读取TPQ的物理量
    void SS_Corr_TPQ(const string &TwoBodyPath, my_mat &SxSxCorr, my_mat &SySyCorr, my_mat &SzSzCorr); // 读取TPQ的物理量
};

#endif // TPQ_H