#include "writer.hpp"
#include "htslib/hts.h"

//#define DEBUG_WRITER

namespace stoat {

Writer::Writer(const std::string output_file_path, size_t thread_count, size_t max_buffer_length) : 
    file_path(output_file_path),
    max_buffer_length(max_buffer_length) {

        // Initialize the buffer with empty strings and allocate memory for them
        buffers.resize(thread_count+1);
        for (size_t i = 0 ; i < buffers.size(); i++) {
            buffers.at(i).reserve(max_buffer_length);
        }
    };

std::string Writer::get_file_path() const {
    return file_path;
}

bool Writer::write(const std::string out_content) {
    bool success = true;
    buffers.at(omp_get_thread_num()).append(out_content);
    if (buffers.at(omp_get_thread_num()).size() >= max_buffer_length) {
        success = flush();
    }
    return success;
}

bool Writer::flush(size_t thread_num) {
    // Call the virtual function to write the buffer
    if (thread_num == std::numeric_limits<size_t>::max()) {
        thread_num = omp_get_thread_num();
    }

    bool written;
    #pragma omp critical (writer)
    {
        written = write_string_to_file(buffers.at(thread_num));
    }
    // clear the buffer 
    buffers.at(omp_get_thread_num()).clear();
    return written;
}

void Writer::close() {
    // Flush the buffer
    for (size_t i = 0 ; i < buffers.size() ; i++) {
        flush(i);
    }
    // Call the virtual function to close the file
    close_file();
}

// Uncompressed writer using a standard output file stream
StdWriter::StdWriter(const std::string output_file_path, size_t max_buffer_length) : Writer(output_file_path, max_buffer_length) {
    // do nothing if the file path is null
    if (output_file_path == "") {
        return;
    }
    // open steam
    file_stream.open(file_path);
}

bool StdWriter::write_string_to_file(const std::string& out_content) {
    // Write to the output file. Since this is always called within the omp guards from write(), this doesn't need its own guards (I think)
    file_stream << out_content;
    return true;
}



void StdWriter::close_file() {
    // close steam
    file_stream.close();
}

// Bgzipped writer using a HTSlib
BgzWriter::BgzWriter(const std::string output_file_path, size_t max_buffer_length) : Writer(output_file_path, max_buffer_length) {
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

bool BgzWriter::write_string_to_file(const std::string& out_content) {
    // Write to the output file. Since this is always called within the omp guards from write(), this doesn't need its own guards
    if (bgzf_write(file_p, out_content.c_str(), out_content.size()) < 0) {
        stoat::LOG_ERROR("Error writing to BGZ output: " + file_path + "\n");
        bgzf_close(file_p);
        return false;
    }
    return true;
}

void BgzWriter::close_file() {
    // Flush the rest of the buffer
    flush();
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
        write("#CHR\tSTART_OFFSET\tEND_OFFSET\tSTART_NODE\tEND_NODE\tALLELE_LENGTHS\tP_FISHER\tP_CHI2\tALLELE_COUNT_PER_PHENO\tDEPTH\n");
    } else if (phenotype_type == stoat::BINARY_COVAR) {
        write("#CHR\tSTART_OFFSET\tEND_OFFSET\tSTART_NODE\tEND_NODE\tALLELE_LENGTHS\tP\tALLEL_COUNT\tDEPTH\n");
    } else if (phenotype_type == stoat::QUANTITATIVE) {
        write("#CHR\tSTART_OFFSET\tEND_OFFSET\tSTART_NODE\tEND_NODE\tALLELE_LENGTHS\tP\tALLEL_COUNT\tDEPTH\n");
    } else if (phenotype_type == stoat::EQTL) {
        write("#CHR\tSTART_OFFSET\tEND_OFFSET\tSTART_NODE\tEND_NODE\tALLELE_LENGTHS\tGENE\tP\tALLEL_COUNT\tDEPTH\n");
    }
}

void Writer::write_binary(const stoat::snarl_info_t& snarl_data, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles) {
    // snarl information
    std::string outl = snarl_data.ref_path;
    outl += "\t" + std::to_string(snarl_data.start_position) + "\t" + std::to_string(snarl_data.end_position);
    outl += "\t" + snarl_data.start_node.to_string() + "\t" + snarl_data.end_node.to_string();
    // allele counts
    if (snarl_data.walks_by_allele.empty()) {
        //TODO: This writes "NA" if there are no paths/path lengths for the alleles. idk if this is what we want to do
        outl += "\tNA";
    } else {
        outl += "\t" + stoat::vectorPathToString(snarl_data.walks_by_allele, true, active_alleles);
    }
    // pvalues and contingency table
    outl += "\t" + stoat::set_precision(test_result.pv) + "\t" + stoat::set_precision(test_result.second_pv) + "\t" + test_result.group_paths + "\t" +
        std::to_string(snarl_data.depth) + "\n";
    // write
    write(outl);
}

void Writer::write_quantitative(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles) {
    // snarl information
    std::string outl = snarl_data.ref_path;
    outl += "\t" + std::to_string(snarl_data.start_position) + "\t" + std::to_string(snarl_data.end_position);
    outl += "\t" + snarl_data.start_node.to_string() + "\t" + snarl_data.end_node.to_string();
    // allele counts
    if (snarl_data.walks_by_allele.empty()) {
        //TODO: This writes "NA" if there are no paths/path lengths for the alleles. idk if this is what we want to do
        outl += "\tNA";
    } else {
        outl += "\t" + stoat::vectorPathToString(snarl_data.walks_by_allele, true, active_alleles);
    }
    // pvalues and contingency table
    outl += "\t" + stoat::set_precision(test_result.pv) + "\t" + test_result.allele_paths + "\t" + std::to_string(snarl_data.depth) + "\n";
    // write
    write(outl);
}

void Writer::write_binary_covar(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles) {
    // now it's the same output
    write_quantitative(snarl_data, test_result, active_alleles);
}

void Writer::write_eqtl(const snarl_info_t& snarl_data, const std::string& gene_name, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles) {
    // snarl information
    std::string outl = snarl_data.ref_path;
    outl += "\t" + std::to_string(snarl_data.start_position) + "\t" + std::to_string(snarl_data.end_position);
    outl += "\t" + snarl_data.start_node.to_string() + "\t" + snarl_data.end_node.to_string();
    // allele counts
    if (snarl_data.walks_by_allele.empty()) {
        //TODO: This writes "NA" if there are no paths/path lengths for the alleles. idk if this is what we want to do
        outl += "\tNA";
    } else {
        outl += "\t" + stoat::vectorPathToString(snarl_data.walks_by_allele, true, active_alleles);
    }
    // pvalues and contingency table
    outl += "\t" + gene_name + "\t" + stoat::set_precision(test_result.pv) + "\t" + test_result.allele_paths + "\t" + std::to_string(snarl_data.depth) + "\n";
    // write
    write(outl);
}

}//end namespace

