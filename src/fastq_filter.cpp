#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <pybind11/pybind11.h>


bool containsSubsequence(const std::string& sequence, const std::string& subsequence, int start, int end) {
    if (sequence.length() < end) return false;
    for (int i = start - 1; i <= end; ++i) {
        if (sequence.substr(i, subsequence.length()) == subsequence) {
            return true;
        }
    }
    return false;
}

void processFastqFiles(const char* inputFile1, const char* inputFile2, const std::string& output_dir, 
                       const std::string& barcode1, int start1, int end1, 
                       const std::string& barcode2, int start2, int end2) {
    std::string output_mate1 = output_dir + "/" + barcode1 + "_" + barcode2 + "_mate1.fastq.gz";
    std::string output_mate2 = output_dir + "/" + barcode1 + "_" + barcode2 + "_mate2.fastq.gz";
    gzFile inFile1 = gzopen(inputFile1, "rb");
    gzFile inFile2 = gzopen(inputFile2, "rb");
    gzFile outFile1 = gzopen(output_mate1.c_str(), "wb");
    gzFile outFile2 = gzopen(output_mate2.c_str(), "wb");

    if (!inFile1 || !inFile2 || !outFile1 || !outFile2) {
        std::cerr << "Error opening files." << std::endl;
        if (inFile1) gzclose(inFile1);
        if (inFile2) gzclose(inFile2);
        if (outFile1) gzclose(outFile1);
        if (outFile2) gzclose(outFile2);
        return;
    }

    char buffer1[1024], buffer2[1024];
    std::vector<std::string> read1(4), read2(4);
    int lineIndex = 0;
    int count = 0;

    while (gzgets(inFile1, buffer1, sizeof(buffer1)) != Z_NULL && gzgets(inFile2, buffer2, sizeof(buffer2)) != Z_NULL) {
        count += 1;
        if(count % 4000000 == 0) {
            std::cout << "Processed " << count / 4 << " reads" << std::endl;
        }
        read1[lineIndex] = buffer1;
        read2[lineIndex] = buffer2;
        if (++lineIndex == 4) {
            if (containsSubsequence(read1[1], barcode1, start1, end1) && containsSubsequence(read2[1], barcode2, start2, end2)) {
                for (const auto& readLine : read1) {
                    gzputs(outFile1, readLine.c_str());
                }
                for (const auto& readLine : read2) {
                    gzputs(outFile2, readLine.c_str());
                }
            }
            lineIndex = 0;
        }
    }

    gzclose(inFile1);
    gzclose(inFile2);
    gzclose(outFile1);
    gzclose(outFile2);
}


namespace py = pybind11;

PYBIND11_MODULE(cpp, m) {
    m.def("process_fastq_files", &processFastqFiles, "A function to process FASTQ files");
}