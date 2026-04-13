#ifndef WRITER_INCLUDED
#define WRITER_INCLUDED

#include <fstream>
#include <htslib/bgzf.h>

#include "types_and_structs.hpp"
#include "stats_test.hpp"

// Classes to abstract the writing (and reading) of files in stoat. Standard write/read can be done with StdWriter/StdReader,
// while bgzipped content can be written/read using BgzWriter/BgzReader. For Writers, some specialized methods are implemented
// to write specific output. Or the generic *write* function can be used. For Readers, it implements the basic getline function.

namespace stoat {

// Generic Writer class with the specialized output writer functions
class Writer {
public:
    Writer(const std::string output_file_path);

    std::string get_file_path() const;

    // write the given string
    virtual bool write(const std::string out_content) = 0;

    // cleanly close the writer
    virtual void close() = 0;

    // header writer
    void write_node_stoat_output_header(phenotype_type_t phenotype_type);
    void write_snarl_stoat_output_header(phenotype_type_t phenotype_type);

    // stoat snarl output writers
    void write_snarl_binary(const stoat::snarl_info_t& snarl_data, const stoat::test_result_t& test_result);
    void write_snarl_binary_covar(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result);
    void write_snarl_quantitative(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result);
    void write_snarl_eqtl(const snarl_info_t& snarl_data, const std::string& gene_name, const stoat::test_result_t& test_result);

    // stoat node output writers
    void write_node_binary(const stoat::node_info_t& node_data, const stoat::test_result_t& test_result);
    void write_node_binary_covar(const node_info_t& node_data, const stoat::test_result_t& test_result);
    void write_node_quantitative(const node_info_t& node_data, const stoat::test_result_t& test_result);
    void write_node_eqtl(const node_info_t& node_data, const std::string& gene_name, const stoat::test_result_t& test_result);
        
protected:
    const std::string file_path;
    
};

// Standart Writer: uncompressed using a standard output file stream
class StdWriter: public Writer {
public:
    StdWriter(const std::string output_file_path);

    bool write(const std::string out_content);

    void close();
protected:
    std::ofstream file_stream;
    
};

// Bgzipped writer using a HTSlib
class BgzWriter: public Writer {
public:
    BgzWriter(const std::string output_file_path);

    bool write(const std::string out_content);

    void close();
protected:
    BGZF *file_p;    
};

// Generic Reader class
class Reader {
public:
    Reader(const std::string input_file_path);

    // read a line and update the *line* argument. Returns false if EOF.
    virtual bool getline(std::string& line) = 0;

    virtual void close() = 0;
protected:
    const std::string file_path;
};

// Standart Reader for uncompressed content using a standard output file stream
class StdReader: public Reader {
public:
    StdReader(const std::string input_file_path);

    // read a line and update the *line* argument. Returns false if EOF.
    bool getline(std::string& line);

    void close();
protected:
    std::ifstream file_stream;
};

// Bgzipped reader using a HTSlib
class BgzReader: public Reader {
public:
    BgzReader(const std::string input_file_path);

    // read a line and update the *line* argument. Returns false if EOF.
    bool getline(std::string& line);

    void close();
protected:
    BGZF *file_p;    
};

} //end namespace

#endif
