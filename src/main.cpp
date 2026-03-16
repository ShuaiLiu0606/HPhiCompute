#include <common.h>
#include <TPQ.h>
#include <LOBPCG.h>
#include <FullDiag.h>

// 第几激发态
#define SAMPLE_NUM 0

int main()
{
     TimeVar ti0 = timeNow();

     cout << "******************  parameters  ******************" << endl;
     cout << "Sample_num: " << SAMPLE_NUM << endl
          << endl;

     if (1)
     {
          cout << "********************  LOBPCG  ********************" << endl;
          cout << "Calculating LOBPCG ... " << flush << endl
               << endl;

          const string FolderPath = "../inputCG/";
          vector<vector<string>> All_files;
          my_vec ValuesTheta;
          my_mat Values_h;

          GetAllFilesWithPrefix(FolderPath, All_files, ValuesTheta, Values_h);
          int NumTheta = ValuesTheta.size();
          cout << endl;

          //******************copy file*******************************

          const string bond_Ofile = FolderPath + All_files[0][0] + "/check_bond.dat";
          const string bond_ifile = "check_bond.dat";
          CopyFile(bond_Ofile, bond_ifile);

          const string Rspace_Ofile = FolderPath + All_files[0][0] + "/new_Rspace_location.dat";
          const string Rspace_ifile = "new_Rspace_location.dat";
          CopyFile(Rspace_Ofile, Rspace_ifile);

          const string Plaq_Ofile = FolderPath + All_files[0][0] + "/check_plaq.dat";
          const string Plaq_ifile = "check_plaq.dat";
          CopyFile(Plaq_Ofile, Plaq_ifile);

          const string FBZLoca_Ofile = FolderPath + All_files[0][0] + "/FBZ_location.dat";
          const string FBZLoca_ifile = "FBZ_location.dat";
          CopyFile(FBZLoca_Ofile, FBZLoca_ifile);

          const string SBZLoca_Ofile = FolderPath + All_files[0][0] + "/SBZ_location.dat";
          const string SBZLoca_ifile = "SBZ_location.dat";
          CopyFile(SBZLoca_Ofile, SBZLoca_ifile);

          const string MLineLoca_Ofile = FolderPath + All_files[0][0] + "/DSSF_Mline_location.dat";
          const string MLineLoca_ifile = "DSSF_Mline_location.dat";
          CopyFile(MLineLoca_Ofile, MLineLoca_ifile);

          const string MSpaceLoca_Ofile = FolderPath + All_files[0][0] + "/DSSF_Mspace_location.dat";
          const string MSpaceLoca_ifile = "DSSF_Mspace_location.dat";
          CopyFile(MSpaceLoca_Ofile, MSpaceLoca_ifile);

          //*****************************************************************

          ofstream ParaFile("ParaFile.dat", ios::app);
          if (!ParaFile.is_open())
          {
               cerr << "Error!!! Failed to open ParaFile.dat" << endl;
          }

          my_cplx_mat Wp_m; // 存储当前所有文件的基态波函数(会非常浪费内存)
          for (int num1 = 0; num1 < NumTheta; num1++)
          {
               const double theta = ValuesTheta[num1];
               ostringstream oss1;
               oss1 << fixed << setprecision(3) << theta; // 设置固定小数点后三位
               const string FileTheta = oss1.str();       // 转换为字符串

               int Num_h = Values_h[num1].size();
               for (int num2 = 0; num2 < Num_h; num2++)
               {
                    TimeVar ti2 = timeNow();

                    const double h = Values_h[num1][num2];
                    ostringstream oss2;
                    oss2 << fixed << setprecision(3) << h; // 设置固定小数点后三位
                    const string File_h = oss2.str();      // 转换为字符串

                    cout << "theta: " << FileTheta << endl;
                    cout << "h: " << File_h << endl;
                    ParaFile << setw(10) << fixed << setprecision(4) << left << theta << setw(10) << fixed << setprecision(4) << right << h << endl;

                    const string filename = All_files[num1][num2];
                    const string CG_Files = FolderPath + filename + "/output";
                    cout << filename << " has been read" << endl;
                    cout << endl;

                    LOBPCG LOBPCG(CG_Files, FileTheta, File_h, SAMPLE_NUM);

                    my_cplx_vec Phi;
                    LOBPCG.Phi_Read(CG_Files, Phi);

                    // Wp_m.push_back(Phi);

                    //   LOBPCG.EnergySpectrum(CG_Files);

                    //     LOBPCG.Fidelity(Wp_m, NumTheta, Num_h);

                    //=====================================================================
                    /*  LOBPCG.S_Moments(Phi.data());
                     LOBPCG.SS_Correlators(Phi.data());
                     LOBPCG.PMoments(Phi.data());
                     LOBPCG.QMoments(Phi.data());
                     LOBPCG.ChiralityMoments(Phi.data());
                     LOBPCG.WpMoments(Phi.data());
                     LOBPCG.PPCorrelators(Phi.data());
                     LOBPCG.QQCorrelators(Phi.data());
                     LOBPCG.SSSS_Correlators(Phi.data()); */

                    // LOBPCG.WpWpCorrelators(Phi.data());

                    //============================================================================
                    //  use hphi file
                    LOBPCG.S_Corr(CG_Files);
                    LOBPCG.SS_Corr(CG_Files);
                    LOBPCG.SS_PCorr(CG_Files);
                    LOBPCG.SS_QCorr(CG_Files);
                    LOBPCG.SSS_Corr(CG_Files, num2);

                    LOBPCG.SSSS_PPCorr(CG_Files);
                    LOBPCG.SSSS_Corr(CG_Files);
                    LOBPCG.SSSSSS_WpCorr(CG_Files);

                    LOBPCG.Magnetization(CG_Files, num2);

                    vector<vector<int>> Wxs;
                    Wxs.resize(1);
                    Wxs[0] = {9, 10, 19, 18, 5, 6, 23, 22, 13, 12, 3, 2};
                    // LOBPCG.WxMoments(Phi.data(), Wxs);

                    vector<vector<int>> Wys;
                    Wys.resize(2);
                    Wys[0] = {1, 22, 21, 8, 9, 14, 13, 4, 5, 18, 17, 0};
                    Wys[1] = {23, 6, 7, 16, 15, 2, 3, 20, 19, 10, 11, 12};
                    // LOBPCG.WyMoments(Phi.data(), Wys);

                    vector<vector<int>> Wzs;
                    Wzs.resize(2);
                    Wzs[0] = {0, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
                    Wzs[1] = {12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
                    // LOBPCG.WzMoments(Phi.data(), Wzs);

                    // LOBPCG.Halfchain_EE(Phi.data());

                    uint32_t sizeA = 2;
                    uint32_t *sitesA = new uint32_t[sizeA]{0, 1};
                    uint32_t sizeB = 2;
                    uint32_t *sitesB = new uint32_t[sizeB]{2, 3};
                    uint32_t sizeC = 2;
                    uint32_t *sitesC = new uint32_t[sizeC]{4, 5};
                    LOBPCG.TEE(Phi.data(), sizeA, sitesA, sizeB, sitesB, sizeC, sitesC);

                    delete[] sitesA;
                    delete[] sitesB;
                    delete[] sitesC;

                    TimeVar tf2 = timeNow();
                    auto diff2 = durationms(tf2 - ti2) / 1000.0;
                    cout << "LOBPCG Cal done. Took " << left << fixed << setprecision(2) << diff2 << " s." << endl
                         << endl;
               }
          }
     }

     if (1)
     {

          cout << "*********************  cTPQ  *********************" << endl;
          cout << "Calculating cTPQ ... " << flush << endl
               << endl;

          const string FolderPath = "../inputTPQ/";
          vector<vector<string>> All_files;
          my_vec ValuesTheta;
          my_mat Values_h;

          GetAllFilesWithPrefix(FolderPath, All_files, ValuesTheta, Values_h);
          int NumTheta = ValuesTheta.size();

          const string Plaq_Ofile = FolderPath + All_files[0][0] + "/check_plaq.dat";
          const string Plaq_ifile = "check_plaq.dat";
          CopyFile(Plaq_Ofile, Plaq_ifile);

          cout << "NumTheta: " << NumTheta << endl;

          for (int num1 = 0; num1 < NumTheta; num1++)
          {

               int Num_h = Values_h[num1].size(); // 读取所有相同前缀的文件夹，以及磁场
               cout << "Num_h: " << Num_h << endl;
               const double theta = ValuesTheta[num1];
               ostringstream oss1;
               oss1 << fixed << setprecision(3) << theta; // 设置固定小数点后三位
               const string FileTheta = oss1.str();       // 转换为字符串

               cout << "theta: " << theta << endl;

               for (int num2 = 0; num2 < Num_h; num2++)
               {
                    TimeVar ti1 = timeNow();

                    const string filename = All_files[num1][num2];

                    const double h = Values_h[num1][num2];
                    ostringstream oss;
                    oss << fixed << setprecision(3) << h; // 设置固定小数点后三位
                    const string File_h = oss.str();      // 转换为字符串

                    const string TPQ_Files = FolderPath + filename + "/output";
                    cout << filename << " has been read" << endl;

                    TPQ TPQ(TPQ_Files, FileTheta, File_h, SAMPLE_NUM); // 构造函数，从文件中读取参数并赋值给成员变量

                    TPQ.mainBasic(TPQ_Files);

                    TPQ.mainPhys(TPQ_Files);

                    TimeVar tf1 = timeNow();
                    auto diff1 = durationms(tf1 - ti1) / 1000.0;
                    cout << "cTPQ done. Took " << left << fixed << setprecision(2) << diff1 << " s." << endl
                         << endl;
               }
          }
     }

     if (0)
     {

          cout << "********************* FullDiag *********************" << endl;
          cout << "Calculating FullDiag ... " << flush << endl
               << endl;
          const string FolderPath = "../inputFull/";
          vector<vector<string>> All_files;
          my_vec ValuesTheta;
          my_mat Values_h;

          GetAllFilesWithPrefix(FolderPath, All_files, ValuesTheta, Values_h);
          int NumTheta = ValuesTheta.size();

          cout << "NumTheta: " << NumTheta << endl;

          for (int num1 = 0; num1 < NumTheta; num1++)
          {

               int Num_h = Values_h[num1].size(); // 读取所有相同前缀的文件夹，以及磁场
               cout << "Num_h: " << Num_h << endl;
               const double theta = ValuesTheta[num1];
               ostringstream oss1;
               oss1 << fixed << setprecision(3) << theta; // 设置固定小数点后三位
               const string FileTheta = oss1.str();       // 转换为字符串

               cout << "theta: " << theta << endl;

               for (int num2 = 0; num2 < Num_h; num2++)
               {
                    TimeVar ti1 = timeNow();

                    const string filename = All_files[num1][num2];

                    const double h = Values_h[num1][num2];
                    ostringstream oss;
                    oss << fixed << setprecision(3) << h; // 设置固定小数点后三位
                    const string File_h = oss.str();      // 转换为字符串

                    const string FullDiag_Files = FolderPath + filename + "/output";
                    cout << filename << " has been read" << endl;

                    FullDiag FullDiag(FullDiag_Files, FileTheta, File_h); // 构造函数，从文件中读取参数并赋值给成员变量

                    double T = 10.0;
                    FullDiag.Fulldiag_cal(FullDiag_Files, T);

                    TimeVar tf1 = timeNow();
                    auto diff1 = durationms(tf1 - ti1) / 1000.0;
                    cout << "FullDiag done. Took " << left << fixed << setprecision(2) << diff1 << " s." << endl
                         << endl;
               }
          }
     }

     TimeVar tf0 = timeNow();
     auto diff0 = durationms(tf0 - ti0) / 1000.0;
     cout << "All Cal done. Took " << left << fixed << setprecision(2) << diff0 << " s." << endl
          << endl;
     return 0;
}
