#include <common.h>
#include <LOBPCG.h>

void LOBPCG::EnergySpectrum(const string &EnPath)
{

    vector<vector<string>> Matrix;
    string in_file = EnPath + "/zvo_energy.dat";
    ifstream E_file(in_file);
    if (!E_file.is_open())
    {
        cerr << "Function EnergySpectrum, Failed to open file: " << in_file << endl;
        return;
    }

    string line;
    while (getline(E_file, line))
    {
        vector<string> row;
        istringstream iss(line);
        string value;
        while (iss >> value)
        {
            row.push_back(value);
        }
        Matrix.push_back(row);
    }
    E_file.close();
    vector<double> Energies;
    vector<double> Szs;

    for (auto vec : Matrix)
    {
        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (vec[i] == "Energy")
            {
                // 提取能量值
                if (i + 1 < vec.size())
                {
                    double energy_value = stod(vec[i + 1]);
                    Energies.push_back(energy_value);
                }
                break;
            }
            if (vec[i] == "Sz")
            {
                // 提取Sz
                if (i + 1 < vec.size())
                {
                    double Sz = stod(vec[i + 1]);
                    Szs.push_back(Sz);
                }
                break;
            }
        }
    }

    const string EnFile = "Energy_All.dat";
    ofstream EnOutFile(EnFile, ios::app); // 第一列是theta， 第二列是磁场，后面是能量
    EnOutFile << File_theta << "    " << File_h << "    ";

    const string SpecFile = "EnSpectrum_All.dat";
    ofstream SpectOutFile(SpecFile, ios::app);
    SpectOutFile << File_theta << "    " << File_h << "    ";

    const string SzFile = "SzbyEn_All.dat";
    ofstream SzOutFile(SzFile, ios::app);
    SzOutFile << File_theta << "    " << File_h << "    ";

    double FirstEnergy = Energies[0];

    for (size_t i = 0; i < Energies.size(); i++)
    {
        EnOutFile << setprecision(16) << setw(26) << fixed << right << Energies[i]; // 每一行存入同一个态的所有能量

        if (i < 3)
        {
            cout << "Energy" << i << ":" << right << setw(12) << setprecision(8) << fixed << Energies[i] << endl;
            cout << "    Sz" << i << ":" << right << setw(12) << setprecision(8) << fixed << Szs[i] << endl;
        }

        SpectOutFile << setw(24) << fixed << right << setprecision(16) << Energies[i] - FirstEnergy;
        SzOutFile << setw(24) << fixed << right << setprecision(16) << Szs[i];
    }
    cout << endl;

    EnOutFile << endl;
    SpectOutFile << endl;
    SzOutFile << endl;

    EnOutFile.close();
    SpectOutFile.close();
}
