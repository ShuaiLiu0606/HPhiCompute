#include <common.h>
#include <LOBPCG.h>
#include <filesystem>

LOBPCG::LOBPCG(const string &ParaFilePath, const string file_theta, const string file_h, int sample_num) : Nsites(0), Sdim(0), Hdim(0), Sample_num(sample_num), File_theta(file_theta), File_h(file_h)
{
    // 读取modpara.def文件里面的参数

    const string ParaF = ParaFilePath + "/../modpara.def";

    string AA = "Nsite";        // 关键字
    GetPara(ParaF, AA, Nsites); // 获取格点数

    const string LocspnFile = ParaFilePath + "/../locspn.def";

    string EE = to_string(Nsites - 1); // 关键字
    int S2 = 0;                        // 2S
    GetPara(LocspnFile, EE, S2);

    Spin = S2 / 2.0;
    Sdim = S2 + 1;
    Hdim = pow(Sdim, Nsites);

    cout << "Nsite: " << Nsites << endl;
    cout << "Spin: " << Spin << endl;
    cout << "Sdim: " << Sdim << endl;
    cout << "Hdim: " << Hdim << endl;
    cout << endl;
}


void LOBPCG::Fidelity(my_cplx_mat Phi_Mat, int NumTheta, int Numfile)
{
    if (NumTheta != 1 && Numfile != 1)
    {
        cout << "Function Fidelity; NumTheta != 1 and Numfile != 1" << endl;
        return;
    }
    // 计算保真度
    // Numfile: 当前theta下磁场h的数量
    my_mat Para;
    Mat_Read("ParaFile.dat", Para); // 第一列是 theta, 后面是磁场

    int numWp = Phi_Mat.size(); // 当前theta下的所有的波函数

    int numphi = 0;

    if (NumTheta == 1)
    {
        numphi = Numfile;
    }
    else
    {
        numphi = NumTheta;
    }

    if (numWp == numphi)
    {
        const string path = "fidelity_All_" + to_string((int)Sample_num) + ".dat";

        ofstream fidelity_F(path);

        for (int i = 0; i < (numWp - 1); i++)
        {
            my_cplx TmpFide = Vec_Dot(Hdim, Phi_Mat[i].data(), 1, Phi_Mat[i + 1].data(), 1);
            double resu = norm(TmpFide);

            double theta = (Para[i][0] + Para[i + 1][0]) / 2.0;
            double h = (Para[i][1] + Para[i + 1][1]) / 2.0;
            fidelity_F << setw(10) << setprecision(4) << fixed << left << theta;
            fidelity_F << setw(10) << setprecision(4) << fixed << left << h;
            fidelity_F << setw(24) << setprecision(16) << fixed << right << resu << endl;
        }
    }
    else
    {
        return;
    }
}

//********************************************************************************************************************

/* void LOBPCG::Phi_Read(const string &FilePath, my_cplx_vec &EigenVec)
{
    // 读取波函数
    // num 表示需要读取的波函数的数目,

    for (int rank = 0; rank < 1; rank++)
    {
        const string FileName = FilePath + "/zvo_eigenvec_" + to_string(Sample_num) + "_rank_" + to_string(rank) + ".dat";

        ifstream file(FileName, ios::binary);
        if (!file.is_open())
        {
            cerr << "Function Phi_Read, Failed to open: " << FileName << endl;
            return;
        }

        file.ignore(4);
        uint64_t num_states;
        file.read(reinterpret_cast<char *>(&num_states), sizeof(uint64_t));

        EigenVec.resize(num_states);

        file.ignore(16);

        for (uint64_t i = 0; i < num_states; i++)
        {
            double re, im;
            file.read(reinterpret_cast<char *>(&re), sizeof(double));
            file.read(reinterpret_cast<char *>(&im), sizeof(double));
            EigenVec[i] = complex<double>(re, im);
        }
    }
}
 */

void LOBPCG::Phi_Read(const string &FilePath, my_cplx_vec &EigenVec)
{
    int max_rank = -1;
    string base_pattern = "zvo_eigenvec_" + to_string(Sample_num) + "_rank_";
    string suffix = ".dat";
    vector<string> valid_files;

    // 收集所有有效文件并确定最大rank
    for (const auto &entry : std::filesystem::directory_iterator(FilePath))
    {
        string filename = entry.path().filename().string();
        if (filename.find(base_pattern) != string::npos && filename.find(suffix) != string::npos)
        {
            size_t start_pos = filename.find(base_pattern) + base_pattern.length();
            size_t end_pos = filename.find(suffix);
            string rank_str = filename.substr(start_pos, end_pos - start_pos);

            try
            {
                int rank = stoi(rank_str);
                valid_files.push_back(entry.path().string()); // 保存完整路径
                if (rank > max_rank)
                    max_rank = rank;
            }
            catch (...)
            { /* 忽略转换错误 */
            }
        }
    }

    if (max_rank == -1)
    {
        cerr << "Function Phi_Read, No valid files in: " << FilePath << endl;
        return;
    }

    // 按rank排序文件（确保顺序正确）
    sort(valid_files.begin(), valid_files.end(), [&](const string &a, const string &b)
         { return stoi(a.substr(a.find(base_pattern) + base_pattern.length(),
                                a.find(suffix) - (a.find(base_pattern) + base_pattern.length()))) <
                  stoi(b.substr(b.find(base_pattern) + base_pattern.length(),
                                b.find(suffix) - (b.find(base_pattern) + base_pattern.length()))); });

    // 第一次遍历：计算总数据量并验证一致性
    uint64_t total_size = 0;
    vector<uint64_t> num_states_list;
    for (const auto &filepath : valid_files)
    {
        ifstream file(filepath, ios::binary);
        if (!file.is_open())
        {
            cerr << "Failed to open: " << filepath << endl;
            return;
        }

        file.ignore(4);
        uint64_t num_states;
        file.read(reinterpret_cast<char *>(&num_states), sizeof(num_states));
        num_states_list.push_back(num_states);
        total_size += num_states;
        file.close();
    }

    // 验证总数据量是否符合预期
    if (total_size != Hdim)
    {
        cerr << "ERROR: Hdim of EigenVec (" << total_size << ") != Hdim (" << Hdim << ")\n";
        return;
    }

    // 一次性调整容器大小
    EigenVec.resize(total_size);

    // 第二次遍历：读取数据到正确位置
    uint64_t offset = 0;
    for (size_t i = 0; i < valid_files.size(); ++i)
    {
        ifstream file(valid_files[i], ios::binary);
        file.ignore(4 + sizeof(uint64_t) + 16); // 跳过头部

        // 直接读取到容器内存
        file.read(reinterpret_cast<char *>(&EigenVec[offset]),
                  num_states_list[i] * sizeof(complex<double>));

        offset += num_states_list[i];
        file.close();
    }
}