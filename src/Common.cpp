#include <common.h>
#include <filesystem>
#include <regex>
#include <set>
#include <sstream>

std::random_device rd;                             // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd());                            // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> uni_rd(0.0, 1.0); // call uni_rd to get a random number in [0,1)

//*********************************************************************************************

void GetAllFilesWithPrefix(const string &folder_path, vector<vector<string>> &All_files, vector<double> &extractedValues1, vector<vector<double>> &extractedValues2)
{
    // 读取所有命名为：KG_theta-?PI_h-?的文件名，并将 theta 的值存进 extractedValues1，每个 theta 对应 h 存进 extractedValues2
    // folder_path: 文件夹的路径
    // All_files： 每一行存储 theta 相同的文件夹

    try
    {
        if (!std::filesystem::exists(folder_path) || !std::filesystem::is_directory(folder_path))
        {
            cerr << "Error: Path does not exist or is not a directory.\n";
            return;
        }

        // 使用 map 以 theta 作为 key，存储 (h, 文件名) 的 vector
        map<double, vector<pair<double, string>>> theta_map;

        for (const auto &entry : std::filesystem::directory_iterator(folder_path))
        {
            string filename = entry.path().filename().string();

            if (entry.is_directory())
            {
                // 使用正则表达式提取 theta- 和 h- 的数值
                std::regex pattern(R"(KG_theta([-+]?\d*\.?\d+)PI_h([-+]?\d*\.?\d+))");
                std::smatch match;
                if (std::regex_search(filename, match, pattern))
                {
                    double theta = stod(match[1]);
                    double h = stod(match[2]);

                    // 存入 map
                    theta_map[theta].emplace_back(h, filename);
                }
            }
        }

        // 提取排序后的数据
        for (auto &[theta, file_list] : theta_map)
        {
            // 按 h 排序
            std::sort(file_list.begin(), file_list.end());

            // 存储该 theta 对应的所有文件名
            vector<string> row;
            vector<double> h_values; // 存储当前 theta 对应的所有 h 值
            for (const auto &[h, fname] : file_list)
            {
                row.push_back(fname);
                h_values.push_back(h);
            }

            All_files.push_back(row);
            extractedValues1.push_back(theta);
            extractedValues2.push_back(h_values);

            // 输出当前 theta 对应的 h 的数量
            cout << "Theta: " << setprecision(3) << fixed << theta << " has " << h_values.size() << " h values." << endl;
        }
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        cerr << "Filesystem error: " << e.what() << "\n";
    }
    catch (const std::exception &e)
    {
        cerr << "Error: " << e.what() << "\n";
    }
}

//****************** 通过HPhi计算关联函数 ********************************************************

Operator multiply_terms(const Operator &terms, const Operator &factor)
{
    Operator new_terms;

    for (const auto &[val, coeff] : terms)
    {
        for (const auto &[f_val, f_coeff] : factor)
        {
            new_terms.emplace_back(val + f_val, coeff * f_coeff);
        }
    }
    return new_terms;
}

// -----------------------------------------------------------------------------
// 构造自旋升算符 S+
//
// Bogoliubov 表示:
//
//  S⁺ = Σ_{σ=-S}^{S-1} √[S(S+1) - σ(σ+1)]  c†_{σ+1} c_{σ}
//
// 指标映射关系:
//      σ = i - S
// 其中
//      i = 0 ... 2S
// 例如 (S=1):
//      σ : -1   0   1
//      i :  0   1   2
// 算符:
//      c†_{σ+1} c_{σ}  ->  key = (i, i+1)
// -----------------------------------------------------------------------------

Operator build_Splus(double Spin)
{
    Operator Sp;
    int Sdim = int(2 * Spin + 1);
    for (int i = 0; i < Sdim - 1; ++i)
    {
        double sigma = i - Spin; // 湮灭的 σ 值
        double coeff = sqrt(Spin * (Spin + 1) - sigma * (sigma + 1));
        string key = to_string(i) + to_string(i + 1); // 产生 i，湮灭 i+1
        Sp.emplace_back(key, coeff);
    }
    return Sp;
}
// -----------------------------------------------------------------------------
// 构造自旋降算符 S-
// Bogoliubov 表示:
//  S⁻ = Σ_{σ=-S}^{S-1} √[S(S+1) - σ(σ+1)]  c†_{σ} c_{σ+1}
// 算符映射:
//      c†_{σ} c_{σ+1}  -> key = (i+1, i)
// -----------------------------------------------------------------------------

// 正确的 S⁻：产生指标 i+1（对应 σ），湮灭指标 i（对应 σ+1）
Operator build_Sminus(double Spin)
{
    Operator Sm;
    int Sdim = int(2 * Spin + 1);
    for (int i = 0; i < Sdim - 1; ++i)
    {
        double sigma = i - Spin; // 产生的 σ 值
        double coeff = sqrt(Spin * (Spin + 1) - sigma * (sigma + 1));
        string key = to_string(i + 1) + to_string(i); // 产生 i+1，湮灭 i
        Sm.emplace_back(key, coeff);
    }
    return Sm;
}
// -----------------------------------------------------------------------------
// 构造 Sz 算符
// Bogoliubov 表示:
//  Sz = Σ_{σ=-S}^{S} σ c†_{σ} c_{σ}
// 算符映射:
//      c†_{σ} c_{σ} -> key = (i, i)
// -----------------------------------------------------------------------------

Operator build_Sz(double Spin)
{
    Operator Sz;
    int Sdim = int(2 * Spin + 1);

    // 遍历所有 σ 从 -S 到 S
    for (int i = 0; i < Sdim; i++)
    {
        double sigma = Spin - i;                  // 修改这里// 当前 σ 值，即本征值
        string key = to_string(i) + to_string(i); // 产生和湮灭指标相同，表示粒子数算符
        Sz.emplace_back(key, sigma);
    }

    /*
       cout << "Test: Spin 1/2 " << endl;
       cout << "Sz: " << Sz[0].first << ", " << Sz[0].second << endl;
       cout << "Sz: " << Sz[1].first << ", " << Sz[1].second << endl;
    */


    return Sz;
}

// -----------------------------------------------------------------------------
// 构造 Sx 算符
//      Sx = (S⁺ + S⁻) / 2
// -----------------------------------------------------------------------------
Operator build_Sx(double Spin)
{
    Operator Sx;
    auto Sp = build_Splus(Spin);
    auto Sm = build_Sminus(Spin);

    // 添加 S^+ 项：系数为 1/2
    for (auto &p : Sp)
    {
        Sx.emplace_back(p.first, 0.5 * p.second);
    }

    // 添加 S^- 项：系数也为 1/2
    for (auto &p : Sm)
    {
        Sx.emplace_back(p.first, 0.5 * p.second);
    }

    // cout << "Test: Spin 1/2 " << endl;
    // cout << "Sx: " << Sx[0].first << ", " << Sx[0].second << endl;
    // cout << "Sx: " << Sx[1].first << ", " << Sx[1].second << endl;


    return Sx;
}

// -----------------------------------------------------------------------------
// 构造 Sy 算符
//      Sy = (S⁺ - S⁻) / (2i)
// -----------------------------------------------------------------------------

Operator build_Sy(double Spin)
{
    Operator Sy;
    auto Sp = build_Splus(Spin);
    auto Sm = build_Sminus(Spin);

    // 添加 S^+ 项：系数 -i/2，即 (0, -0.5)
    for (auto &p : Sp)
    {
        Sy.emplace_back(p.first, complex<double>(0, -0.5) * p.second);
    }

    // 添加 S^- 项：系数 i/2，即 (0, 0.5)
    for (auto &p : Sm)
    {
        Sy.emplace_back(p.first, complex<double>(0, 0.5) * p.second);
    }

    //  cout << "Test: Spin 1/2 " << endl;
    //  cout << "Sy: " << Sy[0].first << ", " << Sy[0].second << endl;
    //  cout << "Sy: " << Sy[1].first << ", " << Sy[1].second << endl;

    return Sy;
}

// 提取计算单体关联函数的通用逻辑
void ComS_Corr(int Sdim, const Operator &result, const my_mat &OneBody, int site, double &realPart, double &imagPart)
{
    for (const auto &[value, coeff] : result)
    {
        if (value.size() != 2)
            continue;

        int s1 = value[0] - '0';
        int s2 = value[1] - '0';

        if (s1 < 0 || s1 >= Sdim || s2 < 0 || s2 >= Sdim)
            continue;

        int s = s1 * Sdim + s2;
        int row = site * (Sdim * Sdim) + s;

        if (row < 0 || row >= (int)OneBody.size())
            continue;

        complex<double> ob(OneBody[row][4], OneBody[row][5]);

        complex<double> contrib = coeff * ob;
        realPart += contrib.real();
        imagPart += contrib.imag();
    }
}

void ComSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &TwoBody, int sitei, int sitej, double &realPart, double &imagPart)
{
    for (const auto &[value, coeff] : result)
    {
        if (value.size() != 4)
            continue;

        int s1a = value[0] - '0';
        int s1b = value[1] - '0';
        int s2a = value[2] - '0';
        int s2b = value[3] - '0';

        if (s1a >= Sdim || s1b >= Sdim || s2a >= Sdim || s2b >= Sdim)
            continue;

        int s1 = s1a * Sdim + s1b;
        int s2 = s2a * Sdim + s2b;

        int row = sitei * (Sdim * Sdim * Nsites * Sdim * Sdim) + s1 * (Nsites * Sdim * Sdim) + sitej * Sdim * Sdim + s2;

        if (row < 0 || row >= (int)TwoBody.size())
            continue;

        complex<double> ob(TwoBody[row][8], TwoBody[row][9]);

        complex<double> contrib = coeff * ob;

        realPart += contrib.real();
        imagPart += contrib.imag();
    }
}

void ComSSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &ThreeBody, int sitei, int sitej, int sitek, double &realPart, double &imagPart)
{
    // 必须 i<j<k
    if (!(0 <= sitei && sitei < sitej && sitej < sitek && sitek < Nsites))
        return;

    // 计算组合索引 C_ijk
    int part1 = (Nsites - sitei - 1) * (Nsites - sitei - 2) / 2;
    int part2 = (Nsites - sitej) * (Nsites - sitej - 1) / 2;
    int part3 = sitek - sitej - 1;

    int C_ijk = part1 - part2 + part3;

    const int S2 = Sdim * Sdim;
    const int S4 = S2 * S2;
    const int S6 = S4 * S2;

    for (const auto &[value, coeff] : result)
    {
        if (value.size() != 6)
            continue;

        int s1a = value[0] - '0';
        int s1b = value[1] - '0';
        int s2a = value[2] - '0';
        int s2b = value[3] - '0';
        int s3a = value[4] - '0';
        int s3b = value[5] - '0';

        if (s1a >= Sdim || s1b >= Sdim || s2a >= Sdim || s2b >= Sdim || s3a >= Sdim || s3b >= Sdim)
            continue;

        int s1 = s1a * Sdim + s1b;
        int s2 = s2a * Sdim + s2b;
        int s3 = s3a * Sdim + s3b;

        int offset = s1 * S4 + s2 * S2 + s3;

        int row = C_ijk * S6 + offset;

        if (row < 0 || row >= (int)ThreeBody.size())
            continue;

        std::complex<double> val(ThreeBody[row][12], ThreeBody[row][13]);

        std::complex<double> contrib = coeff * val;

        realPart += contrib.real();
        imagPart += contrib.imag();
    }
}

void ComSSSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &FourBody, int sitei, int sitej, int sitem, int siten, double &realPart, double &imagPart)
{
    realPart = 0.0;
    imagPart = 0.0;

    // 从文件读取键信息
    const string bondF = "check_bond.dat";
    my_mat bonds;
    Mat_Read(bondF, bonds);
    const int NumBond = bonds.size(); // 键总数

    // 查找键索引b1和b2
    int b1 = -1, b2 = -1;
    for (int b = 0; b < NumBond; ++b)
    {
        const uint32_t i = bonds[b][1];
        const uint32_t j = bonds[b][2];
        if ((sitei == i && sitej == j) || (sitei == j && sitej == i))
        {
            b1 = b;
            break;
        }
    }
    for (int b = 0; b < NumBond; ++b)
    {
        const uint32_t k = bonds[b][1];
        const uint32_t l = bonds[b][2];
        if ((sitem == k && siten == l) || (sitem == l && siten == k))
        {
            b2 = b;
            break;
        }
    }

    if (b1 == -1 || b2 == -1)
    {
        cerr << "Error: Invalid bond pair (" << sitei << "," << sitej << ") or (" << sitem << "," << siten << ")\n";
        return;
    }

    const int Sdim2 = Sdim * Sdim;
    const int Sdim4 = Sdim2 * Sdim2;
    const int Sdim8 = Sdim4 * Sdim4; // 每个键对的数据块大小

    const int block_index = (b1 * NumBond + b2) * Sdim8;

    for (const auto &[value, coeff] : result)
    {
        if (value.size() != 8)
        {
            cerr << "Invalid spin format: " << value << endl;
            continue;
        }

        int sigma[8];
        for (int i = 0; i < 8; ++i)
        {
            sigma[i] = value[i] - '0';
            if (sigma[i] < 0 || sigma[i] >= Sdim)
            {
                cerr << "Invalid spin index: " << sigma[i] << endl;
                continue;
            }
        }

        // 计算自旋组合索引，与build_fourBodyG顺序一致
        const int s1 = sigma[0] * Sdim + sigma[1];
        const int s2 = sigma[2] * Sdim + sigma[3];
        const int s3 = sigma[4] * Sdim + sigma[5];
        const int s4 = sigma[6] * Sdim + sigma[7];

        const int s_index = ((s1 * Sdim2 + s2) * Sdim2 + s3) * Sdim2 + s4;
        const int row = block_index + s_index;

        if (row < 0 || row >= (int)FourBody.size())
        {
            cerr << "Row " << row << " out of range" << endl;
            continue;
        }
        if (FourBody[row].size() < 18) // 确保有第16、17列
        {
            cerr << "FourBody row " << row << " has insufficient columns" << endl;
            continue;
        }

        complex<double> ob(FourBody[row][16], FourBody[row][17]);
        complex<double> contrib = coeff * ob;
        realPart += contrib.real();
        imagPart += contrib.imag();
    }
}

void ComSSSSSS_Corr(int Nsites, int Sdim, const Operator &result, const my_mat &SixBody, int plaq_index, double &realPart, double &imagPart)
{
    realPart = 0.0;
    imagPart = 0.0;

    const int Sdim2 = Sdim * Sdim;
    const int block = (int)pow(Sdim2, 6); // rows per plaquette

    for (const auto &[value, coeff] : result)
    {
        // 必须是 12 个自旋指标
        if (value.size() != 12)
        {
            cerr << "Invalid spin string length: " << value << endl;
            continue;
        }

        bool valid = true;
        int sigma[12];

        // 解析 12 个 σ
        for (int i = 0; i < 12; i++)
        {
            sigma[i] = value[i] - '0';
            if (sigma[i] < 0 || sigma[i] >= Sdim)
            {
                cerr << "Invalid spin index in " << value << endl;
                valid = false;
                break;
            }
        }
        if (!valid)
            continue;

        // 6 个复合自旋
        int s1 = sigma[0] * Sdim + sigma[1];
        int s2 = sigma[2] * Sdim + sigma[3];
        int s3 = sigma[4] * Sdim + sigma[5];
        int s4 = sigma[6] * Sdim + sigma[7];
        int s5 = sigma[8] * Sdim + sigma[9];
        int s6 = sigma[10] * Sdim + sigma[11];

        // 行索引（与 build_sixBodyG 一致）
        int row = plaq_index * block +
                  s1 * (int)pow(Sdim2, 5) +
                  s2 * (int)pow(Sdim2, 4) +
                  s3 * (int)pow(Sdim2, 3) +
                  s4 * (int)pow(Sdim2, 2) +
                  s5 * Sdim2 +
                  s6;

        if (row < 0 || row >= (int)SixBody.size())
        {
            cerr << "Row out of range: " << row << endl;
            continue;
        }
        if (SixBody[row].size() < 26) // 确保有第 24、25 列
        {
            cerr << "SixBody row " << row << " has insufficient columns" << endl;
            continue;
        }

        complex<double> ob(SixBody[row][24], SixBody[row][25]);
        complex<double> contrib = coeff * ob;
        realPart += contrib.real();
        imagPart += contrib.imag();
    }
}
//****************** 随机向量归一 *********************

void rd_wf(const int &dim, double *wf)
{
    // 生成随机向量
    for (int i = 0; i < dim; i++)
    {
        wf[i] = (uni_rd(gen) - 0.5);
    }
    // normalize
    cblas_dscal(dim, 1.0 / cblas_dnrm2(dim, wf, 1), wf, 1);
}

void rd_wf(const int &dim, my_cplx *wf)
{
    // 生成随机向量
    my_cplx II(0, 1);
    for (int i = 0; i < dim; i++)
    {
        wf[i] = (uni_rd(gen) - 0.5) + (uni_rd(gen) - 0.5) * II;
    }
    // normalize
    cblas_zdscal(dim, 1.0 / cblas_dznrm2(dim, wf, 1), wf, 1);
}

void rd_wf_product(const int &dim, double *wf)
{
    // 生成随机product state
    for (int i = 0; i < dim; i++)
    {
        wf[i] = 0;
    }
    int idx = uni_rd(gen) * dim;
    wf[idx] = 1.0;
    // cout << "norm of random product states"  << cblas_dnrm2(dim, wf, 1) << endl;
}

void uni_wf(const int &dim, double *wf)
{
    // 生成均匀向量
    for (int i = 0; i < dim; i++)
    {
        wf[i] = 1;
    }
    // normalize
    cblas_dscal(dim, 1.0 / cblas_dnrm2(dim, wf, 1), wf, 1);
}

//********************** bitwise operations ********************
int Bits_Count1(const int &b)
{
    // 计算整数b的二进制表示中1的个数
    int a = b;
    int cnt = 0;
    while (a != 0)
    {
        ++cnt;
        a &= (a - 1);
    }
    return cnt;
}

int Bits_CycleRight(const int &s, const int &n)
{
    // 将整数s的二进制表示循环右移n位

    // return (s >> 1) + ((s & 1) << (n-1));
    return (s >> 1) | (s << (n - 1));
}

int Bits_CycleLeft(const int &s, const int &n)
{
    // 将整数s的二进制表示循环左移n位

    // return ((s << 1) | (s >> (n - 1))) - ((s & (1 << (n - 1))) << 1);
    return ((s << 1) | (s >> (n - 1))) - ((s >> (n - 1)) << n);
}

int Bits_Reflect(const int &s, const int &n)
{
    // 将整数s的二进制表示左右翻转
    int aux = 0;
    for (int i = 0; i < n; i++)
    {
        aux |= (aux | (((s >> i) & 1) << (n - 1 - i)));
    }
    return aux;
}

int Bits_Invert(const int &s, const int &n)
{
    // 将整数s的二进制表示取反
    int aux = (1 << n) - 1;
    return s ^ aux;
}

int Bits_HammingDis(const int &a, const int &b)
{
    // 计算两个整数a和b的汉明距离
    return Bits_Count1(a ^ b);
}

void Bits_Decomposition(const uint32_t &s, const uint32_t &n, const uint32_t &size_A, uint32_t *sites_A, uint32_t &a, uint32_t &b)
{
    // 其主要作用是将一个整数 s 按照一组指定的位索引 sites_A 进行分解，分解成两部分 a 和 b。
    // 用n是格点数Nsites
    // size_A表示sites_A中有几个个位置
    // sites_A = {1, 3} (表示要将 s 的第 1 和第 3 位放入 a 中)
    uint32_t aux_B = 0;
    //
    a = 0;
    b = 0;
    for (uint32_t ia = 0; ia < size_A; ia++)
    {
        // 将 s 右移 sites_A[ia] 位； 取出第 sites_A[ia] 位的值(0或1），将这个值放到 a 的第 ia 位，
        // 右移：0101->0010(要将最右边的数字丢掉，最左边补0)，左移一样
        a += (((s >> sites_A[ia]) & 1) << ia); // 取出s的第sites_A[ia])位置。&：同1为1，否则为0
        aux_B += (1 << sites_A[ia]);
    }
    // aux_B 用来记录 sites_A 中指定的位,例如(00101)，1表示指定的位，之后通过取反操作来得到那些未被指定的位。（通过二进制的方式）
    // 例如 n=4, size_A = 2， sites_A = {1, 3}， aux_B首先变为2(0010)，然后变为2+8=10(1010)
    // 1010表示系统A中标记了第1，3位

    uint32_t ib = 0;
    aux_B = (~aux_B);
    // 取反操作，将二进制中的0变为1，1变为0
    for (uint32_t i = 0; i < n; i++)
    {
        if ((aux_B >> i) & 1)
        { // 检查 aux_B 的第 i 位是否为 1。如果是，则表示 s 的第 i 位应该被分配到 b 中。
          // 这里的i应该是前面ia的补集
            b += ((s >> i) & 1) << ib;
            // 将 s 右移 i 位； 取出第 i 位的值(1或0)，将这个值放到 b 的第 ib 位
            ib++;
        }
    }

    // 例如：s = 011010, sites_A = {0,4}
    // 则：a = {1,0}(2), b = {0,1,0,1}(5)
}

void Bits_Print(const int &a, const int &n)
{
    // 输出整数a的二进制表示

    cout << "  ";
    for (int ix = 0; ix < n; ix++)
        cout << ((a >> (n - ix - 1)) & 1);
    cout << "  ";
}

int nchoosek(const int &n, const int &_k)
{
    // 计算组合数C(n,k)

    if (_k > n)
    {
        return 0;
    }

    int k = _k < (n - _k) ? _k : (n - _k);
    if (0 == k)
    {
        return 1;
    }
    if (1 == k)
    {
        return n;
    }
    double aux = 1.0;
    for (int i = 0; i < k; i++)
    {
        aux *= double(n - i) / double(k - i);
    }

    return (int)(aux + 1e-2);
}

void Debug_Check_Conj_dense(double *mat, int dim)
{
    // 检查矩阵是否是厄米矩阵
    double aux = 0;
    for (int i = 0; i < dim; i++)
    {
        for (int j = i + 1; j < dim; j++)
        {
            aux += abs(mat[i * dim + j] - mat[j * dim + i]);
        }
    }
    cout << "abs(H - H^T) = " << aux << endl;
}

//***************** 文件读取 **********************

void GetPara(const string &fname, const string &string_match, int &para)
{
    // 从文件fname中读取string_match开头的行，并将int类型的参数赋值给para
    ifstream f_in(fname);
    if (!f_in)
    {
        cerr << "Function GetPara; Failed to open: " << fname << endl;

        return;
    }
    string line;
    while (getline(f_in, line))
    {
        if (line.find(string_match) == 0)
        {
            // 判断行是否以 string_match 开头
            string::size_type pos = string_match.length();
            while (pos < line.length() && isspace(line[pos]))
            {
                pos++; // 跳过空格
            }
            if (pos < line.length())
            {
                stringstream ss(line.substr(pos)); // 提取参数部分
                ss >> para;                        // 将字符串转换为 int 类型并赋值给 para
            }
            break; // 找到匹配行后退出循环
        }
    }

    f_in.close(); // 关闭文件
}

void GetPara(const string &fname, const string &string_match, double &para)
{
    // 从文件fname中读取string_match开头的行，并将double类型的参数赋值给para
    ifstream f_in(fname);
    if (!f_in)
    {
        cerr << "Failed to open: " << fname << endl;
        return;
    }
    string line;
    while (getline(f_in, line))
    {
        if (line.find(string_match) == 0)
        {
            string::size_type pos = string_match.length();
            while (pos < line.length() && isspace(line[pos]))
            {
                pos++; // 跳过空格
            }
            if (pos < line.length())
            {
                stringstream ss(line.substr(pos)); // 提取参数部分
                ss >> para;                        // 将字符串转换为 double 类型并赋值给 para
            }
            break; // 找到匹配行后退出循环
        }
    }

    f_in.close(); // 关闭文件
}

void GetPara(const string &fname, const string &string_match, char &para)
{
    // 从文件fname中读取string_match开头的行，并将char类型的参数赋值给para
    ifstream f_in(fname);
    if (!f_in)
    {
        cerr << "Failed to open: " << fname << endl;
        return;
    }
    string line;
    while (getline(f_in, line))
    {
        if (line.find(string_match) == 0)
        {
            string::size_type pos = string_match.length();
            while (pos < line.length() && isspace(line[pos]))
            {
                pos++; // 跳过空格
            }
            if (pos < line.length())
            {
                para = line[pos]; // 提取非空格字符作为 char 参数
            }
            break; // 找到匹配行后退出循环
        }
    }

    f_in.close(); // 关闭文件
}

//****************** 向量,矩阵的读取，输出 **************************
void Vec_Read(const string &fname, vector<complex<double>> &data)
{
    // 从文件fname中读取复数向量

    ifstream f_in(fname);
    if (!f_in)
    {
        cerr << "Failed to open: " << fname << endl;
        return;
    }
    data.clear();

    double real, imag;
    while (f_in >> real >> imag)
    {
        data.emplace_back(real, imag);
    }
    f_in.close();
}

void Vec_Read(const string &fname, vector<double> &data)
{
    // 从文件fname中读取实数向量

    ifstream f_in(fname);
    if (!f_in)
    {
        cerr << "Failed to open: " << fname << endl;
    }
    data.clear();

    double aa;
    while (f_in >> aa)
    {
        data.emplace_back(aa);
    }

    f_in.close();
}

void Vec_Print(const complex<double> *data, const int dsize)
{
    // 输出复数向量
    for (int i = 0; i < dsize; ++i)
    {
        double realPart = data[i].real();
        double imagPart = data[i].imag();

        cout << setprecision(8) << "(" << right << setw(11) << realPart << "," << right << setw(11) << imagPart << ")  ";
    }
}

void Vec_Print(const double *data, const int dsize)
{
    // 输出实数向量
    for (int i = 0; i < dsize; i++)
    {
        cout << setprecision(8) << setw(11) << right << data[i] << "  ";
    }
}

// mat
void Mat_Read(const string &filename, vector<vector<int>> &matrix)
{
    // 从文件filename中读取实数矩阵
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed open to file: " << filename << endl;
        return;
    }

    matrix.clear();

    string line;
    while (getline(file, line))
    {
        vector<int> row;
        istringstream iss(line);
        int value;
        while (iss >> value)
        {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    file.close();
}

void Mat_Read(const string &filename, vector<vector<double>> &matrix)
{
    // 从文件filename中读取实数矩阵
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed open to file: " << filename << endl;
        return;
    }

    matrix.clear();

    string line;
    while (getline(file, line))
    {
        vector<double> row;
        istringstream iss(line);
        double value;
        while (iss >> value)
        {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    file.close();
}

void Mat_Read(const string &filename, vector<vector<complex<double>>> &matrix)
{
    // 从文件filename中读取复数矩阵

    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed open to file: " << filename << endl;
        return;
    }

    matrix.clear(); // 关键：清空旧数据，防止内存累积

    string line;
    while (getline(file, line))
    {
        vector<complex<double>> row;
        istringstream iss(line);
        string value_str;

        // 从行中读取每个复数值
        while (iss >> value_str)
        {
            // 去除括号
            value_str.erase(remove(value_str.begin(), value_str.end(), '('), value_str.end());
            value_str.erase(remove(value_str.begin(), value_str.end(), ')'), value_str.end());

            // 找到逗号的位置
            size_t comma_pos = value_str.find(',');

            if (comma_pos != string::npos)
            {
                try
                {
                    // 提取实部和虚部字符串
                    string real_str = value_str.substr(0, comma_pos);
                    string imag_str = value_str.substr(comma_pos + 1);

                    // 转换为 double 并创建复数
                    double real = stod(real_str);
                    double imag = stod(imag_str);

                    complex<double> value(real, imag);
                    row.push_back(value);
                }
                catch (const std::invalid_argument &e)
                {
                    cerr << "Invalid number format in: " << value_str << endl;
                }
                catch (const std::out_of_range &e)
                {
                    cerr << "Number out of range in: " << value_str << endl;
                }
            }
        }

        matrix.push_back(row);
    }
    file.close();
}

void Mat_Print(const complex<double> *data, const int rows, const int cols)
{
    // 输出复数矩阵
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {

            cout << setprecision(8) << fixed << right << data[i * cols + j] << "  ";
        }
        cout << endl;
    }
}

void Mat_Print(const double *data, const int rows, const int cols)
{
    // 输出实数矩阵
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << setprecision(8) << fixed << right << data[i * cols + j] << "  ";
        }
        cout << endl;
    }
}

void Mat_Pointer(const vector<vector<double>> &matrix, double *flat_array)
{
    // 将二维矩阵转换为一维数组
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            flat_array[i * cols + j] = matrix[i][j];
        }
    }
}

void Mat_Pointer(const vector<vector<complex<double>>> &matrix, complex<double> *flat_array)
{
    // 将二维复数矩阵转换为一维复数数组
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            flat_array[i * cols + j] = matrix[i][j];
        }
    }
}

//*******************二进制向量的写入，读取，输出***************************

// 二进制写入
void Vec_Fwrite(const string &fname, const double *data, const int dsize)
{
    // 将实数向量写入二进制文件
    ofstream f_out(fname.c_str(), ios::binary);
    if (!f_out)
    {
        cerr << "Failed to write: " << fname << endl;
        return;
    }
    f_out.write(reinterpret_cast<const char *>(data), dsize * sizeof(double));
    f_out.close();
}

void Vec_Fwrite(const string &fname, const complex<double> *data, const int dsize)
{
    // 将复数向量写入二进制文件
    ofstream f_out(fname.c_str(), ios::binary);
    if (!f_out)
    {
        cerr << "Failed to write: " << fname << endl;
        return;
    }
    f_out.write(reinterpret_cast<const char *>(data), dsize * sizeof(complex<double>));
    f_out.close();
}

// 二进制读取
void Vec_Fread(const string &fname, double *data, const int dsize)
{
    // 从二进制文件中读取实数向量
    ifstream f_in(fname.c_str(), ios::binary);
    if (!f_in)
    {
        cerr << "Failed to open: " << fname << endl;
        return;
    }
    f_in.read(reinterpret_cast<char *>(data), dsize * sizeof(double));
    f_in.close();
}

void Vec_Fread(const string &fname, complex<double> *data, const int dsize)
{
    // 从二进制文件中读取复数向量
    ifstream f_in(fname.c_str(), ios::binary);
    if (!f_in)
    {
        cerr << "Failed to open: " << fname << endl;
        return;
    }
    f_in.read(reinterpret_cast<char *>(data), dsize * sizeof(complex<double>));
    f_in.close();
}

//*****************************向量矩阵运算****************************
double Vec_Dot(const int &dim, double *x, const int incx, double *y, const int incy)
{
    // 计算两个实数向量的点积
    return cblas_ddot(dim, x, incx, y, incy);
}

my_cplx Vec_Dot(const int &dim, my_cplx *x, const int incx, my_cplx *y, const int incy)
{
    // 计算两个复数向量的点积
    my_cplx aux;
    cblas_zdotc_sub(dim, x, incx, y, incy, &aux);
    return aux;
}

void Out_Dot(int row, const double *vec1, int col, const double *vec2, double *result)
{
    // 计算外积
    double alpha = 1.0;

    cblas_dger(CblasRowMajor, row, col, alpha, vec1, 1, vec2, 2, result, col);
}

void Out_Dot(int row, const complex<double> *vec1, int col, const complex<double> *vec2, complex<double> *result)
{
    // 计算外积
    const complex<double> alpha(1.0, 0.0);
    cblas_zgerc(CblasRowMajor, row, col, &alpha, vec1, 1, vec2, 1, result, col);
}

void Mat_Vec_Dot(const double *matrix, const double *vec, int rows, int cols, double *result)
{
    // 计算矩阵和向量的乘积
    const double alpha = 1.0;
    const double beta = 0.0;

    cblas_dgemv(CblasRowMajor, CblasNoTrans, rows, cols, alpha, matrix, cols, vec, 1, beta, result, 1);
}

void Mat_Vec_Dot(const my_cplx *matrix, const my_cplx *vec, int rows, int cols, my_cplx *result)
{
    // 计算矩阵和向量的乘积
    const complex<double> alpha(1.0, 0.0);
    const complex<double> beta(0.0, 0.0);

    cblas_zgemv(CblasRowMajor, CblasNoTrans, rows, cols, &alpha, matrix, cols, vec, 1, &beta, result, 1);
}

void Mat_Dot(CBLAS_LAYOUT matrix_layout, CBLAS_TRANSPOSE A_transpose, CBLAS_TRANSPOSE B_transpose, int M, int N, int K, my_cplx *A, int lda, my_cplx *B, int ldb, my_cplx *C, int ldc)
{
    // 计算矩阵乘法

    my_cplx alpha = {1.0, 0.0}; // 乘数因子
    my_cplx beta = {0.0, 0.0};  // 加法因子
    cblas_zgemm(matrix_layout, A_transpose, B_transpose, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc);
    // M：结果矩阵 C(A) 的行数。N：结果矩阵 C(B) 的列数。
    // K：矩阵 A 的列数（或矩阵 B 的行数）。
    // lda：矩阵 A 的前导维度。ldb：矩阵 B 的前导维度。ldc：矩阵 C 的前导维度。(列数)
}

void Mat_Dot(CBLAS_LAYOUT matrix_layout, CBLAS_TRANSPOSE A_transpose, CBLAS_TRANSPOSE B_transpose, int M, int N, int K, double *A, int lda, double *B, int ldb, double *C, int ldc)
{
    // 计算矩阵乘法
    double alpha = 1.0;
    double beta = 0.0;
    cblas_dgemm(matrix_layout, A_transpose, B_transpose, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

void Vec_Cross(const double *vec1, const double *vec2, double *result)
{
    // 计算叉积
    fill(result, result + 3, 0.0);
    if (vec1 == nullptr || vec2 == nullptr || result == nullptr)
    {
        cerr << "Error: Null pointer passed to Vec_cross." << std::endl;
        return;
    }

    // 计算叉积
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

void Vec_Cross(const my_cplx *vec1, const my_cplx *vec2, my_cplx *result)
{
    // 计算叉积
    fill(result, result + 3, complex<double>(0.0, 0.0));
    if (vec1 == nullptr || vec2 == nullptr || result == nullptr)
    {
        cerr << "Error: Null pointer passed to Vec_cross." << endl;
        return;
    }

    // 计算叉积
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

double *Kronecker(const double *matA, const double *matB, int A_rows, int A_cols, int B_rows, int B_cols)
{
    // 计算Kronecker积

    // Kronecker积的矩阵大小
    int C_rows = A_rows * B_rows;
    int C_cols = A_cols * B_cols;

    // 分配Kronecker积矩阵空间
    double *matC = new double[C_rows * C_cols];

    // 计算Kronecker积
    for (int i = 0; i < A_rows; i++)
    {
        for (int j = 0; j < A_cols; j++)
        {
            for (int k = 0; k < B_rows; k++)
            {
                for (int l = 0; l < B_cols; l++)
                {
                    // 一维数组模拟二维数组的索引计算
                    matC[(i * B_rows + k) * C_cols + (j * B_cols + l)] = matA[i * A_cols + j] * matB[k * B_cols + l];
                }
            }
        }
    }

    return matC;
}

// 复数矩阵的直乘函数
my_cplx *Kronecker(const my_cplx *matA, const my_cplx *matB, int A_rows, int A_cols, int B_rows, int B_cols)
{
    // 计算Kronecker积

    // Kronecker积的矩阵大小
    int C_rows = A_rows * B_rows;
    int C_cols = A_cols * B_cols;

    // 分配Kronecker积矩阵空间
    my_cplx *matC = new my_cplx[C_rows * C_cols];

    // 计算Kronecker积
    for (int i = 0; i < A_rows; i++)
    {
        for (int j = 0; j < A_cols; j++)
        {
            for (int k = 0; k < B_rows; k++)
            {
                for (int l = 0; l < B_cols; l++)
                {
                    // 一维数组模拟二维数组的索引计算
                    matC[(i * B_rows + k) * C_cols + (j * B_cols + l)] = matA[i * A_cols + j] * matB[k * B_cols + l];
                }
            }
        }
    }

    return matC;
}

//******************************others************************
double VN_Entropy(const int &dim, double *p)
{
    // 计算诺伊曼熵
    double aux = 0;
    int n = 0;
    for (int i = 0; i < dim; i++)
    {
        if (p[i] > 0)
        {
            aux += -p[i] * log(p[i]);
        }
    }

    return aux;
}

pair<double, double> Rotate_Coord(double x, double y, double angle)
{
    // 旋转坐标

    // 正角度表示逆时针旋转
    double rad = angle * PI / 180.0; // 转换为弧度
    double cosTheta = cos(rad);
    double sinTheta = sin(rad);
    return {x * cosTheta - y * sinTheta, x * sinTheta + y * cosTheta};
}

void CopyFile(const string &sourcePath, const string &destinationPath)
{
    // 复制文件

    ifstream source(sourcePath);
    ofstream destination(destinationPath);

    if (!source)
    {
        cerr << "Failed to open: " << sourcePath << endl;
    }

    if (!destination)
    {
        cerr << "Failed to open: " << destinationPath << endl;
    }

    destination << source.rdbuf(); // 将源文件的内容复制到目标文件
}

//****************************有限温cTPQ**********************

// Function to calculate mean along a specific axis for 2D data
my_vec cla_mean(const my_mat &data, int axis)
{
    // 计算2D数据沿特定轴的均值
    int rows = data.size();
    if (rows == 0)
    {
        return {};
    }

    int cols = data[0].size();

    my_vec result;

    if (axis == 0)
    {
        // Compute mean along rows
        result.resize(cols, 0.0);
        for (int j = 0; j < cols; j++)
        {
            double sum = 0.0;
            for (int i = 0; i < rows; i++)
            {
                sum += data[i][j];
            }
            result[j] = sum / rows;
        }
    }
    else if (axis == 1)
    {
        // Compute mean along columns
        result.resize(rows, 0.0);
        for (int i = 0; i < rows; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < cols; j++)
            {
                sum += data[i][j];
            }
            result[i] = sum / cols;
        }
    }
    else
    {
        throw std::invalid_argument("Invalid axis value. Axis must be 0 or 1.");
    }

    return result;
}

// Function to calculate standard deviation along a specific axis for 2D data
my_vec cal_standard_deviation(const my_mat &data, int axis, bool ddof)
{
    // 计算2D数据沿特定轴的标准差
    //  s^2 = (\sum (x_i - x(平均值))^2) / (N-1)
    //  ddof = true;无偏估计假设样本是从总体中随机抽取的，因此调整了分母以纠正偏差
    my_vec means = cla_mean(data, axis);
    my_vec result;

    int rows = data.size();
    if (rows == 0)
    {
        return {};
    }

    int cols = data[0].size();

    if (axis == 0)
    {
        // 0 表示按列计算（即对每列的值求标准差）

        result.resize(cols, 0.0);
        for (int j = 0; j < cols; j++)
        {
            double sum_sq_diff = 0.0;
            for (int i = 0; i < rows; i++)
            {
                sum_sq_diff += pow(data[i][j] - means[j], 2);
            }
            result[j] = sqrt(sum_sq_diff / (ddof ? rows - 1 : rows));
        }
    }
    else if (axis == 1)
    {
        // 1 表示按行计算（即对每行的值求标准差）
        result.resize(rows, 0.0);
        for (int i = 0; i < rows; i++)
        {
            double sum_sq_diff = 0.0;
            for (int j = 0; j < cols; j++)
            {
                sum_sq_diff += pow(data[i][j] - means[i], 2);
            }
            result[i] = sqrt(sum_sq_diff / (ddof ? cols - 1 : cols));
        }
    }
    else
    {
        throw std::invalid_argument("Invalid axis value. Axis must be 0 or 1.");
    }

    return result;
}

// Function to calculate mean along a specific axis for 3D data
my_mat cla_mean(const vector<my_mat> &data, int axis)
{
    // 计算3D数据沿特定轴的均值

    int depth = data.size();
    if (depth == 0)
        return {};
    int rows = data[0].size();
    if (rows == 0)
        return {};
    int cols = data[0][0].size();

    my_mat result;

    if (axis == 0)
    {
        // Calculate mean along the depth
        result.resize(rows, my_vec(cols, 0.0));
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < depth; k++)
                {
                    sum += data[k][i][j];
                }
                result[i][j] = sum / depth;
            }
        }
    }
    else if (axis == 1)
    {
        // Calculate mean along the rows
        result.resize(depth, my_vec(cols, 0.0));
        for (int k = 0; k < depth; k++)
        {
            for (int j = 0; j < cols; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < rows; i++)
                {
                    sum += data[k][i][j];
                }
                result[k][j] = sum / rows;
            }
        }
    }
    else if (axis == 2)
    {
        // Calculate mean along the columns
        result.resize(depth, my_vec(rows, 0.0));
        for (int k = 0; k < depth; k++)
        {
            for (int i = 0; i < rows; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < cols; j++)
                {
                    sum += data[k][i][j];
                }
                result[k][i] = sum / cols;
            }
        }
    }
    else
    {
        throw std::invalid_argument("Invalid axis value. Axis must be 0, 1, or 2.");
    }

    return result;
}

// Function to calculate standard deviation along a specific axis for 3D data
my_mat cal_standard_deviation(const vector<my_mat> &data, int axis, bool ddof)
{
    // 计算3D数据沿特定轴的标准差

    my_mat means = cla_mean(data, axis);
    my_mat result;

    int depth = data.size();
    if (depth == 0)
        return {};
    int rows = data[0].size();
    if (rows == 0)
        return {};
    int cols = data[0][0].size();

    if (axis == 0)
    {
        // Calculate std along the depth
        result.resize(rows, my_vec(cols, 0.0));
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                double sum_sq_diff = 0.0;
                for (int k = 0; k < depth; k++)
                {
                    sum_sq_diff += pow(data[k][i][j] - means[i][j], 2);
                }
                result[i][j] = sqrt(sum_sq_diff / (ddof ? depth - 1 : depth));
            }
        }
    }
    else if (axis == 1)
    {
        // Calculate std along the rows
        result.resize(depth, my_vec(cols, 0.0));
        for (int k = 0; k < depth; k++)
        {
            for (int j = 0; j < cols; j++)
            {
                double sum_sq_diff = 0.0;
                for (int i = 0; i < rows; i++)
                {
                    sum_sq_diff += pow(data[k][i][j] - means[k][j], 2);
                }
                result[k][j] = sqrt(sum_sq_diff / (ddof ? rows - 1 : rows));
            }
        }
    }
    else if (axis == 2)
    {
        // Calculate std along the columns
        result.resize(depth, my_vec(rows, 0.0));
        for (int k = 0; k < depth; k++)
        {
            for (int i = 0; i < rows; i++)
            {
                double sum_sq_diff = 0.0;
                for (int j = 0; j < cols; j++)
                {
                    sum_sq_diff += pow(data[k][i][j] - means[k][i], 2);
                }
                result[k][i] = sqrt(sum_sq_diff / (ddof ? cols - 1 : cols));
            }
        }
    }
    else
    {
        throw std::invalid_argument("Invalid axis value. Axis must be 0, 1, or 2.");
    }

    return result;
}

//******************************Matrix evd (use mkl lapacke function)**********************************

void Mat_Svd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double *a, lapack_int lda, double *s, double *u, lapack_int ldu, double *vt, lapack_int ldvt)
{
    // 使用mkl lapacke函数计算矩阵的奇异值分解
    vector<double> superb(std::min(m, n));
    lapack_int info = LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb.data());
    if (info > 0)
    {
        // 错误处理：LAPACKE_zgesvd未能收敛
        cerr << "SVD failed to converge, error code: " << info << endl;
    }
}

void Mat_Svd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, my_cplx *a, lapack_int lda, double *s, my_cplx *u, lapack_int ldu, my_cplx *vt, lapack_int ldvt)
{
    // 使用mkl lapacke函数计算矩阵的奇异值分解
    vector<double> superb(std::min(m, n));
    lapack_int info = LAPACKE_zgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb.data());
    if (info > 0)
    {
        // 错误处理：LAPACKE_zgesvd未能收敛
        cerr << "SVD failed to converge, error code: " << info << endl;
    }
}

// ===============================================================================
// Matrix evd (use mkl lapacke function)

void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, double *a, lapack_int lda, double *w)
{
    // 使用mkl lapacke函数计算实数对称矩阵的特征值和特征向量
    int eigtype = 1; // 0, dsyevd; 1, dsyev; 2 dsyevr;
    if (0 == eigtype)
    {
        LAPACKE_dsyev(matrix_layout, jobz, uplo, n, a, lda, w);
    }
    if (1 == eigtype)
    {
        LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, a, lda, w);
    }
    if (2 == eigtype)
    {
        char range = 'A'; // all eigenvalues
        lapack_int *isuppz = new lapack_int[2 * n];
        double abstol = 0;
        LAPACKE_dsyevr(matrix_layout, jobz, range, uplo,
                       n, a, lda, NULL, NULL, NULL, NULL,
                       abstol, &n, w, a, n, isuppz);
    }
}

void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, my_cplx *a, lapack_int lda, double *w)
{
    // 使用mkl lapacke函数计算复数对称矩阵的特征值和特征向量

    my_int eigtype = 1; // 0, dsyevd; 1, dsyev; 2 dsyevr;
    if (0 == eigtype)
    {
        LAPACKE_zheev(matrix_layout, jobz, uplo, n, a, lda, w);
    }
    if (1 == eigtype)
    {
        LAPACKE_zheevd(matrix_layout, jobz, uplo, n, a, lda, w);
    }
    // jobz: 'N'：仅计算特征值。'V'：同时计算特征值和特征向量。
    if (2 == eigtype)
    {
        char range = 'A'; // all eigenvalues
        lapack_int *isuppz = new lapack_int[2 * n];
        double abstol = 0;
        LAPACKE_zheevr(matrix_layout, jobz, range, uplo,
                       n, a, lda, NULL, NULL, NULL, NULL,
                       abstol, &n, w, a, n, isuppz);
        free(isuppz);
    }
}

void Lanczos(const int &dim, double *vals, int *cols, int *PointerBE, const int &M, char job_vec, double *alpha, double *l_vecs, double *l_eigvecs)
{
    // Lanczos算法

    sparse_matrix_t A;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, dim, dim, PointerBE, PointerBE + 1, cols, vals);
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    mkl_sparse_optimize(A);
    //
    int steps = M;
    // Lanczos matrix elements
    double *beta = new double[steps]; // beta[0] is unused

    // if we do not need lanzcos vecters, only 3 vectors are used
    if ('N' == job_vec)
    {
        double *phi0 = new double[dim];
        double *phi1 = new double[dim];
        double *phi2 = new double[dim];
        // copy initial vector form 1st dim elements form l_vecs
        cblas_dcopy(dim, l_vecs, 1, phi0, 1);

        // alpha[0], beta[1] -- beta[0] is not required
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi0, 0, phi1);
        alpha[0] = cblas_ddot(dim, phi0, 1, phi1, 1);
        cblas_daxpy(dim, -alpha[0], phi0, 1, phi1, 1);
        beta[1] = cblas_dnrm2(dim, phi1, 1);
        cblas_dscal(dim, 1.0 / beta[1], phi1, 1);

        // alpha[m], beta[m+1], m = [0,steps-2]
        for (int m = 1; m < steps - 1; m++)
        {
            // alpha[m], m = [0,steps-2]
            mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
            alpha[m] = cblas_ddot(dim, phi1, 1, phi2, 1);
            // beta[m+1], m = [0,steps-2]
            cblas_daxpy(dim, -alpha[m], phi1, 1, phi2, 1);
            cblas_daxpy(dim, -beta[m], phi0, 1, phi2, 1);
            beta[m + 1] = cblas_dnrm2(dim, phi2, 1);
            cblas_dscal(dim, 1.0 / beta[m + 1], phi2, 1);
            cblas_dcopy(dim, phi1, 1, phi0, 1);
            cblas_dcopy(dim, phi2, 1, phi1, 1);
        }
        // alpha[m], m = steps-1
        int m = steps - 1;
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
        alpha[m] = cblas_ddot(dim, phi1, 1, phi2, 1);

        LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'N', steps, alpha, &beta[0] + 1, NULL, steps);

        delete[] phi0;
        delete[] phi1;
        delete[] phi2;
    }

    // lanzcos vecters are required
    if ('V' == job_vec)
    {
        // alpha[0], beta[1] -- beta[0] is not required
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[0], 0, &l_vecs[dim]);
        alpha[0] = cblas_ddot(dim, &l_vecs[0], 1, &l_vecs[dim], 1);
        cblas_daxpy(dim, -alpha[0], &l_vecs[0], 1, &l_vecs[dim], 1);
        beta[1] = cblas_dnrm2(dim, &l_vecs[dim], 1);
        cblas_dscal(dim, 1.0 / beta[1], &l_vecs[dim], 1);

        // alpha[0], beta[1] -- beta[0] is not required
        for (int m = 1; m < steps - 1; m++)
        {
            // alpha[m], m = [0,steps-2]
            mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, &l_vecs[(m + 1) * dim]);
            alpha[m] = cblas_ddot(dim, &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1);
            // beta[m+1], m = [0,steps-2]
            cblas_daxpy(dim, -alpha[m], &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1);
            cblas_daxpy(dim, -beta[m], &l_vecs[(m - 1) * dim], 1, &l_vecs[(m + 1) * dim], 1);
            beta[m + 1] = cblas_dnrm2(dim, &l_vecs[(m + 1) * dim], 1);
            cblas_dscal(dim, 1.0 / beta[m + 1], &l_vecs[(m + 1) * dim], 1);
        }
        // alpha[m], m = steps-1
        int m = steps - 1;
        double *tmpvec = new double[dim];
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, tmpvec);
        alpha[m] = cblas_ddot(dim, &l_vecs[m * dim], 1, tmpvec, 1);
        delete[] tmpvec;

        LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);

        // test wf0
        /*
        double* wf0 = new double[dim];
        double* wf1 = new double[dim];
        for (int i = 0; i < dim; i++)
        {
            wf0[i] = Vec_Dot(M, &l_eigvecs[0], M, &l_vecs[i], dim);
        }
        double aux0 = cblas_dnrm2(dim, wf0, 1);
        cout << "aux0: " << aux0 << endl;
        cblas_dscal(dim, 1.0 / aux0, wf0, 1);
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf0, 0, wf1);
        cout << "<wf0|H|wf0> = " << setprecision(14) << cblas_ddot(dim, wf0, 1, wf1, 1) << endl;
        delete[]wf0;
        delete[]wf1;
        */
    }
    delete[] beta;
    mkl_sparse_destroy(A);
}

void Lanczos(const int &dim, my_cplx *vals, int *cols, int *PointerBE, const int &M, char job_vec,
             double *alpha, my_cplx *l_vecs, double *l_eigvecs)
{
    // Lanczos算法
    sparse_matrix_t A;
    mkl_sparse_z_create_csr(&A, SPARSE_INDEX_BASE_ZERO, dim, dim, PointerBE, PointerBE + 1, cols, vals);
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_HERMITIAN;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    mkl_sparse_optimize(A);
    //
    my_int steps = M;
    // Lanczos matrix elements
    double *beta = new double[steps]; // beta[0] is unused

    // if we do not need lanzcos vecters, only 3 vectors are used
    if ('N' == job_vec)
    {
        my_cplx *phi0 = new my_cplx[dim];
        my_cplx *phi1 = new my_cplx[dim];
        my_cplx *phi2 = new my_cplx[dim];
        // copy initial vector form 1st dim elements form l_vecs
        cblas_zcopy(dim, l_vecs, 1, phi0, 1);

        // alpha[0], beta[1] -- beta[0] is not required
        mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi0, 0, phi1);
        my_cplx aux;
        cblas_zdotc_sub(dim, phi0, 1, phi1, 1, &aux);
        alpha[0] = aux.real();
        my_cplx minus_aux = -aux;
        cblas_zaxpy(dim, &minus_aux, phi0, 1, phi1, 1);
        beta[1] = cblas_dznrm2(dim, phi1, 1);
        cblas_zdscal(dim, 1.0 / beta[1], phi1, 1);

        // alpha[m], beta[m+1], m = [0,steps-2]
        for (my_int m = 1; m < steps - 1; m++)
        {
            // alpha[m], m = [0,steps-2]
            mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
            cblas_zdotc_sub(dim, phi1, 1, phi2, 1, &aux);
            alpha[m] = aux.real();
            // beta[m+1], m = [0,steps-2]
            minus_aux = -aux;
            cblas_zaxpy(dim, &minus_aux, phi1, 1, phi2, 1);
            my_cplx minus_beta = -beta[m];
            cblas_zaxpy(dim, &minus_beta, phi0, 1, phi2, 1);
            beta[m + 1] = cblas_dznrm2(dim, phi2, 1);
            cblas_zdscal(dim, 1.0 / beta[m + 1], phi2, 1);
            cblas_zcopy(dim, phi1, 1, phi0, 1);
            cblas_zcopy(dim, phi2, 1, phi1, 1);
        }
        // alpha[m], m = steps-1
        my_int m = steps - 1;
        mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
        cblas_zdotc_sub(dim, phi1, 1, phi2, 1, &aux);
        alpha[m] = aux.real();

        LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'N', steps, alpha, &beta[0] + 1, NULL, steps);

        delete[] phi0;
        delete[] phi1;
        delete[] phi2;
    }

    // lanzcos vecters are required
    // my_cplx* alpha_cmplx = new my_cplx[dim];
    if ('V' == job_vec)
    {
        // alpha[0], beta[1] -- beta[0] is not required
        my_cplx aux;
        mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[0], 0, &l_vecs[dim]);
        cblas_zdotc_sub(dim, &l_vecs[0], 1, &l_vecs[dim], 1, &aux);
        alpha[0] = aux.real();
        my_cplx minus_aux = -aux;
        cblas_zaxpy(dim, &minus_aux, &l_vecs[0], 1, &l_vecs[dim], 1);
        beta[1] = cblas_dznrm2(dim, &l_vecs[dim], 1);
        cblas_zdscal(dim, 1.0 / beta[1], &l_vecs[dim], 1);

        // alpha[m], beta[m+1], m = [0,steps-2]
        for (my_int m = 1; m < steps - 1; m++)
        {
            // alpha[m], m = [0,steps-2]
            mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, &l_vecs[(m + 1) * dim]);
            cblas_zdotc_sub(dim, &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1, &aux);
            alpha[m] = aux.real();
            // alpha[m+1], m = [0,steps-2]
            minus_aux = -aux;
            cblas_zaxpy(dim, &minus_aux, &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1);
            my_cplx minus_beta = -beta[m];
            cblas_zaxpy(dim, &minus_beta, &l_vecs[(m - 1) * dim], 1, &l_vecs[(m + 1) * dim], 1);
            beta[m + 1] = cblas_dznrm2(dim, &l_vecs[(m + 1) * dim], 1);
            cblas_zdscal(dim, 1.0 / beta[m + 1], &l_vecs[(m + 1) * dim], 1);
        }
        // alpha[m], m = steps-1
        my_int m = steps - 1;
        my_cplx *tmpvec = new my_cplx[dim];
        mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, tmpvec);
        cblas_zdotc_sub(dim, &l_vecs[m * dim], 1, tmpvec, 1, &aux);
        delete[] tmpvec;
        alpha[m] = aux.real();

        // LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);
        // for (my_int m = 0; m < steps; m++) alpha[m] = alpha_cmplx[m].real();
        // LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);
        LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);

        // test wf0
        /*
        my_cplx* wf0 = new my_cplx[dim];
        my_cplx* wf1 = new my_cplx[dim];
        for (my_int i = 0; i < dim; i++)
        {
            //wf0[i] = my_vec_dot(M, &l_eigvecs[0], M, &l_vecs[i], dim);
            wf0[i] = 0;
            for (my_int m = 0; m < M; m++)
            {
                wf0[i] += l_eigvecs[m * M] * l_vecs[m * dim + i];
            }
        }
        double aux0 = cblas_dznrm2(dim, wf0, 1);
        cout << "aux0: " << endl;
        cblas_zdscal(dim, 1.0 / aux0, wf0, 1);
        mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf0, 0, wf1);
        cout << "<wf0|H|wf0> = " << setprecision(14) << my_vec_dot(dim, wf0, 1, wf1, 1) << endl;
        delete[]wf0;
        delete[]wf1;
        */
    }

    delete[] beta;
    mkl_sparse_destroy(A);
}

// extra vector operations
void Vec_ax(const int &dim, const my_cplx &a, my_cplx *x, my_cplx *y)
{
    // y = a*x
    for (int i = 0; i < dim; i++)
    {
        y[i] = a * x[i];
    }
}
void Vec_ax(const int &dim, const double &a, double *x, double *y)
{
    // y = a*x
    for (int i = 0; i < dim; i++)
    {
        y[i] = a * x[i];
    }
}

void Vec_real(const int &dim, my_cplx *x, double *a)
{
    // a = real(x)
    for (int i = 0; i < dim; i++)
    {
        a[i] = x[i].real();
    }
}

void Vec_imag(const int &dim, my_cplx *x, double *a)
{
    // a = imag(x)
    for (int i = 0; i < dim; i++)
    {
        a[i] = x[i].imag();
    }
}

void Vec_xpiy(const int &dim, double *x, double *y, my_cplx *a)
{
    // a = x + i*y, x and y are real vectors

    my_cplx II(0, 1);
    for (int i = 0; i < dim; i++)
    {
        a[i] = x[i] + II * y[i];
    }
}
