#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


/**
 * Checks if a given sequence contains a subsequence within a specified range.
 * Uses exact matching (no mismatches allowed).
 *
 * The entire subsequence must fit within [start, end] to match seqkit behavior.
 *
 * @param sequence The sequence to search within.
 * @param subsequence The subsequence to search for.
 * @param start The starting index of the range (1-indexed, inclusive).
 * @param end The ending index of the range (1-indexed, inclusive).
 * @return True if the subsequence is found entirely within the range, false otherwise.
 */
bool containsSubsequence(const std::string& sequence, const std::string& subsequence, int start, int end) {
    if (static_cast<int>(sequence.length()) < end) return false;
    int subseq_len = static_cast<int>(subsequence.length());
    // Convert to 0-indexed. Pattern at position i ends at i+subseq_len-1 (0-indexed).
    // For pattern to fit in [start, end] (1-indexed), we need:
    //   i+1 >= start  =>  i >= start-1
    //   i+subseq_len <= end  =>  i <= end-subseq_len
    for (int i = start - 1; i <= end - subseq_len && i >= 0; ++i) {
        if (sequence.substr(i, subseq_len) == subsequence) {
            return true;
        }
    }
    return false;
}


/**
 * Struct to hold barcode information for batch processing.
 */
struct BarcodeInfo {
    std::string barcode_5p;      // 5' barcode (searches in R2)
    std::string barcode_3p_rc;   // 3' barcode reverse complement (searches in R1)
    int start_5p;                // Start position for 5' barcode (1-indexed)
    int end_5p;                  // End position for 5' barcode (1-indexed)
    int start_3p;                // Start position for 3' barcode (1-indexed)
    int end_3p;                  // End position for 3' barcode (1-indexed)
    std::string output_name;     // Output filename prefix (e.g., "BARCODE1_BARCODE2")
};


/**
 * Processes two FASTQ files and filters the reads based on specified criteria.
 * Original single-barcode version for backward compatibility.
 */
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


/**
 * Two-pass batch processing of FASTQ files for internal demultiplexing.
 *
 * Pass 1: Read all pairs, match against all barcodes, store match index in memory.
 * Pass 2: Write outputs in batches to limit open file handles.
 *
 * @param r1_path Path to R1 FASTQ file (contains 3' barcode, reverse complemented).
 * @param r2_path Path to R2 FASTQ file (contains 5' barcode, forward).
 * @param output_dir Directory for output files.
 * @param barcodes_5p Vector of 5' barcodes.
 * @param barcodes_3p_rc Vector of 3' barcodes (reverse complement).
 * @param starts_5p Vector of start positions for 5' barcodes.
 * @param ends_5p Vector of end positions for 5' barcodes.
 * @param starts_3p Vector of start positions for 3' barcodes.
 * @param ends_3p Vector of end positions for 3' barcodes.
 * @param output_names Vector of output name prefixes.
 * @param max_open_files Maximum number of file handles to open at once.
 * @return Vector of match counts per barcode.
 */
std::vector<int> processFastqFilesBatch(
    const std::string& r1_path,
    const std::string& r2_path,
    const std::string& output_dir,
    const std::vector<std::string>& barcodes_5p,
    const std::vector<std::string>& barcodes_3p_rc,
    const std::vector<int>& starts_5p,
    const std::vector<int>& ends_5p,
    const std::vector<int>& starts_3p,
    const std::vector<int>& ends_3p,
    const std::vector<std::string>& output_names,
    int max_open_files = 200
) {
    size_t num_barcodes = barcodes_5p.size();
    std::vector<int> match_counts(num_barcodes, 0);

    if (num_barcodes == 0) {
        std::cerr << "No barcodes provided" << std::endl;
        return match_counts;
    }

    // === PASS 1: Read all pairs and find matches ===
    std::cout << "Pass 1: Finding barcode matches..." << std::endl;

    gzFile inR1 = gzopen(r1_path.c_str(), "rb");
    gzFile inR2 = gzopen(r2_path.c_str(), "rb");

    if (!inR1 || !inR2) {
        std::cerr << "Error opening input files" << std::endl;
        if (inR1) gzclose(inR1);
        if (inR2) gzclose(inR2);
        return match_counts;
    }

    // Store match index for each read (-1 = no match)
    std::vector<int32_t> matches;
    matches.reserve(1000000);  // Pre-allocate for ~1M reads

    char buffer1[1024], buffer2[1024];
    std::vector<std::string> read1(4), read2(4);
    int lineIndex = 0;
    int total_reads = 0;

    while (gzgets(inR1, buffer1, sizeof(buffer1)) != Z_NULL &&
           gzgets(inR2, buffer2, sizeof(buffer2)) != Z_NULL) {
        read1[lineIndex] = buffer1;
        read2[lineIndex] = buffer2;

        if (++lineIndex == 4) {
            total_reads++;
            if (total_reads % 1000000 == 0) {
                std::cout << "  Pass 1: Processed " << total_reads << " reads" << std::endl;
            }

            // Find matching barcode
            int match_idx = -1;
            for (size_t i = 0; i < num_barcodes; ++i) {
                // 5' barcode in R2, 3' barcode (RC) in R1
                if (containsSubsequence(read2[1], barcodes_5p[i], starts_5p[i], ends_5p[i]) &&
                    containsSubsequence(read1[1], barcodes_3p_rc[i], starts_3p[i], ends_3p[i])) {
                    match_idx = static_cast<int>(i);
                    match_counts[i]++;
                    break;  // First match wins
                }
            }
            matches.push_back(match_idx);
            lineIndex = 0;
        }
    }

    gzclose(inR1);
    gzclose(inR2);

    std::cout << "Pass 1 complete: " << total_reads << " reads, "
              << matches.size() << " stored" << std::endl;

    // Count total matches
    int total_matched = 0;
    for (size_t i = 0; i < num_barcodes; ++i) {
        total_matched += match_counts[i];
    }
    std::cout << "Total matched: " << total_matched << " ("
              << (100.0 * total_matched / total_reads) << "%)" << std::endl;

    // === PASS 2: Write outputs in batches ===
    std::cout << "Pass 2: Writing output files..." << std::endl;

    // Calculate batch size based on max_open_files (2 files per barcode)
    int batch_size = max_open_files / 2;
    if (batch_size < 1) batch_size = 1;

    int num_batches = (static_cast<int>(num_barcodes) + batch_size - 1) / batch_size;

    for (int batch = 0; batch < num_batches; ++batch) {
        int batch_start = batch * batch_size;
        int batch_end = std::min(batch_start + batch_size, static_cast<int>(num_barcodes));

        std::cout << "  Batch " << (batch + 1) << "/" << num_batches
                  << ": barcodes " << batch_start << "-" << (batch_end - 1) << std::endl;

        // Open output files for this batch
        std::vector<gzFile> outR1(batch_end - batch_start);
        std::vector<gzFile> outR2(batch_end - batch_start);

        for (int i = batch_start; i < batch_end; ++i) {
            std::string out_r1 = output_dir + "/" + output_names[i] + "_mate1.fastq.gz";
            std::string out_r2 = output_dir + "/" + output_names[i] + "_mate2.fastq.gz";
            outR1[i - batch_start] = gzopen(out_r1.c_str(), "wb");
            outR2[i - batch_start] = gzopen(out_r2.c_str(), "wb");

            if (!outR1[i - batch_start] || !outR2[i - batch_start]) {
                std::cerr << "Error opening output files for barcode " << i << std::endl;
            }
        }

        // Re-read input files and write matching reads
        gzFile inR1_pass2 = gzopen(r1_path.c_str(), "rb");
        gzFile inR2_pass2 = gzopen(r2_path.c_str(), "rb");

        if (!inR1_pass2 || !inR2_pass2) {
            std::cerr << "Error reopening input files for pass 2" << std::endl;
            continue;
        }

        lineIndex = 0;
        int read_idx = 0;

        while (gzgets(inR1_pass2, buffer1, sizeof(buffer1)) != Z_NULL &&
               gzgets(inR2_pass2, buffer2, sizeof(buffer2)) != Z_NULL) {
            read1[lineIndex] = buffer1;
            read2[lineIndex] = buffer2;

            if (++lineIndex == 4) {
                int match_idx = matches[read_idx];

                // Check if this read's match is in current batch
                if (match_idx >= batch_start && match_idx < batch_end) {
                    int local_idx = match_idx - batch_start;
                    for (const auto& line : read1) {
                        gzputs(outR1[local_idx], line.c_str());
                    }
                    for (const auto& line : read2) {
                        gzputs(outR2[local_idx], line.c_str());
                    }
                }

                read_idx++;
                lineIndex = 0;
            }
        }

        gzclose(inR1_pass2);
        gzclose(inR2_pass2);

        // Close output files for this batch
        for (int i = 0; i < batch_end - batch_start; ++i) {
            if (outR1[i]) gzclose(outR1[i]);
            if (outR2[i]) gzclose(outR2[i]);
        }
    }

    std::cout << "Pass 2 complete" << std::endl;

    return match_counts;
}


namespace py = pybind11;

PYBIND11_MODULE(cpp, m) {
    m.doc() = "C++ module for fast FASTQ processing and demultiplexing";

    m.def("process_fastq_files", &processFastqFiles,
          "Process FASTQ files with single barcode pair (legacy)",
          py::arg("input_file1"),
          py::arg("input_file2"),
          py::arg("output_dir"),
          py::arg("barcode1"),
          py::arg("start1"),
          py::arg("end1"),
          py::arg("barcode2"),
          py::arg("start2"),
          py::arg("end2"));

    m.def("process_fastq_files_batch", &processFastqFilesBatch,
          "Process FASTQ files with multiple barcode pairs (batch mode)",
          py::arg("r1_path"),
          py::arg("r2_path"),
          py::arg("output_dir"),
          py::arg("barcodes_5p"),
          py::arg("barcodes_3p_rc"),
          py::arg("starts_5p"),
          py::arg("ends_5p"),
          py::arg("starts_3p"),
          py::arg("ends_3p"),
          py::arg("output_names"),
          py::arg("max_open_files") = 200);
}
