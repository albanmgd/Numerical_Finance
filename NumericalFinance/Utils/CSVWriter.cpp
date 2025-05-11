#include "CSVWriter.h"
#include <fstream>
#include <iostream>
#include <vector>

void WriteCSV(const std::vector<std::vector<double> > &data, const std::string &fileName) {
    std::ofstream file(fileName, std::ios::trunc);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << fileName << "the file is already open" << std::endl;
        return;
    } // if the file is already open we raise an error and cancel
    // otherwise we write the header
    file << "meanPrice, var, nbSimul\n";
    // we then write the rows
    for (const auto&row : data) {
        if (row.size() != 3) {
            std::cerr << "wrong number of rows, skipping" << std::endl;
            continue;
        }
        file << row[0] << "," << row[1] << "," << row[2] << "\n";
    }
    file.close();
    std::cout << "Done writing csv, file location: " << fileName << std::endl;
}

