#include <SFML/Graphics.hpp>

#include <Python.h>

#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

void query_files(std::vector<std::string> &to_carve);

int main(int argc, char** argv) {
    std::vector<std::string> to_carve;
    query_files(to_carve);

    for (int i = 0; i < to_carve.size(); i++) {
        std::cout << to_carve[i] << "\n";
    }
}

void query_files(std::vector<std::string> &to_carve) {
    const std::string py_file = "query_file.py";

    FILE *f = fopen(py_file.c_str(), "r");

    Py_Initialize();

    std::string convert = "./seamcarver";
    std::wstring widestr(convert.begin(), convert.end());
    wchar_t * argv[1];
    argv[0] = const_cast<wchar_t*>(widestr.c_str());
    PySys_SetArgv(1, argv);
    
    PyRun_SimpleFile(f, py_file.c_str());
    Py_Finalize();

    std::ifstream fin("toseamcarve.txt");
    std::string line;
    while (getline(fin, line)) {
        to_carve.push_back(line);
    }
}
