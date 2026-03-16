#include "writer.hpp"
#include "htslib/hts.h"

//#define DEBUG_WRITER

namespace stoat {

Writer::Writer(const std::string output_file_path) : file_path(output_file_path) {};

std::string Writer::get_file_path() const {
    return file_path;
}

// Uncompressed writer using a standard output file stream
StdWriter::StdWriter(const std::string output_file_path) : Writer(output_file_path) {
    // do nothing if the file path is null
    if (output_file_path == "") {
        return;
    }
    // open steam
    file_stream.open(file_path);
}

bool StdWriter::write(const std::string out_content) {
    #pragma omp critical (writer)
    {
        file_stream << out_content;
    }
    return true;
}

void StdWriter::close() {
    // close steam
    file_stream.close();
}

// Bgzipped writer using a HTSlib
BgzWriter::BgzWriter(const std::string output_file_path) : Writer(output_file_path) {
    // do nothing if the file path is null
    if (output_file_path == "") {
        return;
    }
    // open hts connections
    file_p = bgzf_open(output_file_path.c_str(), "w6");
    if (!file_p) {
        stoat::LOG_ERROR("Error opening file for writing: " + output_file_path + "\n");
    }
}

bool BgzWriter::write(const std::string out_content) {
    if (bgzf_write(file_p, out_content.c_str(), out_content.size()) < 0) {
        stoat::LOG_ERROR("Error writing to BGZ output: " + file_path + "\n");
        bgzf_close(file_p);
        return 0;
    }
    return 1;
}

void BgzWriter::close() {
    // close steam
    bgzf_close(file_p);
}

// Reader management
Reader::Reader(const std::string input_file_path): file_path(input_file_path) {}

// Standard reader using infile stream
StdReader::StdReader(const std::string input_file_path): Reader(input_file_path) {
    // do nothing if the file path is null
    if (input_file_path == "") {
        return;
    }
    // open hts connections
    file_stream.open(input_file_path);
}

bool StdReader::getline(std::string& line) {
    if(std::getline(file_stream, line)) {
        return true;
    } else {
        return false;
    }
}

void StdReader::close() {
    // close steam
    file_stream.close();
}

// Bgzipped TSV writer using a HTSlib
BgzReader::BgzReader(const std::string input_file_path): Reader(input_file_path) {
    // do nothing if the file path is null
    if (input_file_path == "") {
        return;
    }
    // open hts connections
    file_p = bgzf_open(input_file_path.c_str(), "r");
    if (!file_p) {
        stoat::LOG_ERROR("Error opening file for reading: " + input_file_path + "\n");
    }
}

bool BgzReader::getline(std::string& line) {
    kstring_t tmp_kline = KS_INITIALIZE;
    int ret = bgzf_getline(file_p, '\n', &tmp_kline);
    if (ret < 0) {
        return false;
    }
    std::string tmp_line(ks_str(&tmp_kline), ks_len(&tmp_kline));
    line = tmp_line;
    return ret >= 0;
}

void BgzReader::close() {
    // close steam
    bgzf_close(file_p);
}

void Writer::write_stoat_output_header(stoat::phenotype_type_t phenotype_type) {

    if (phenotype_type == stoat::BINARY) {
        write("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tGROUP_PATHS\tDEPTH\n");
    } else if (phenotype_type == stoat::BINARY_COVAR) {
        write("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP\tALLELE_PATHS\tDEPTH\n");
    } else if (phenotype_type == stoat::QUANTITATIVE) {
        write("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP\tALLELE_PATHS\tDEPTH\n");
    } else if (phenotype_type == stoat::EQTL) {
        write("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tGENE\tP\tALLELE_PATHS\tDEPTH\n");
    }
}

void Writer::write_binary(const stoat::snarl_info_t& snarl_data, const stoat::test_result_t& test_result) {
    // snarl information
    std::string outl = snarl_data.ref_path;
    outl += "\t" + std::to_string(snarl_data.start_position) + "\t" + std::to_string(snarl_data.end_position);
    outl += "\t" + stoat::pairToString(std::make_pair(snarl_data.start_node.get_node_id(), snarl_data.end_node.get_node_id()));
    // allele counts
    if (snarl_data.walks_by_allele.empty()) {
        //TODO: This writes "NA" if there are no paths/path lengths for the alleles. idk if this is what we want to do
        outl += "\tNA";
    } else {
        outl += "\t" + stoat::vectorPathToString(snarl_data.walks_by_allele, true);
    }
    // pvalues and contingency table
    outl += "\t" + test_result.pv + "\t" + test_result.second_pv + "\t" + test_result.group_paths + "\t" +
        std::to_string(snarl_data.depth) + "\n";
    // write
    write(outl);
}

void Writer::write_quantitative(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result) {
    // snarl information
    std::string outl = snarl_data.ref_path;
    outl += "\t" + std::to_string(snarl_data.start_position) + "\t" + std::to_string(snarl_data.end_position);
    outl += "\t" + stoat::pairToString(std::make_pair(snarl_data.start_node.get_node_id(), snarl_data.end_node.get_node_id()));
    // allele counts
    if (snarl_data.walks_by_allele.empty()) {
        //TODO: This writes "NA" if there are no paths/path lengths for the alleles. idk if this is what we want to do
        outl += "\tNA";
    } else {
        outl += "\t" + stoat::vectorPathToString(snarl_data.walks_by_allele, true);
    }
    // pvalues and contingency table
    outl += "\t" + test_result.pv + "\t" + test_result.allele_paths + "\t" + std::to_string(snarl_data.depth) + "\n";
    // write
    write(outl);
}

void Writer::write_binary_covar(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result) {
    // now it's the same output
    write_quantitative(snarl_data, test_result);
}

void Writer::write_eqtl(const snarl_info_t& snarl_data, const std::string& gene_name, const stoat::test_result_t& test_result) {
    // snarl information
    std::string outl = snarl_data.ref_path;
    outl += "\t" + std::to_string(snarl_data.start_position) + "\t" + std::to_string(snarl_data.end_position);
    outl += "\t" + stoat::pairToString(std::make_pair(snarl_data.start_node.get_node_id(), snarl_data.end_node.get_node_id()));
    // allele counts
    if (snarl_data.walks_by_allele.empty()) {
        //TODO: This writes "NA" if there are no paths/path lengths for the alleles. idk if this is what we want to do
        outl += "\tNA";
    } else {
        outl += "\t" + stoat::vectorPathToString(snarl_data.walks_by_allele, true);
    }
    // pvalues and contingency table
    outl += "\t" + gene_name + "\t" + test_result.pv + "\t" + test_result.allele_paths + "\t" + std::to_string(snarl_data.depth) + "\n";
    // write
    write(outl);
}

}//end namespace

