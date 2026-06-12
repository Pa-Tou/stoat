#ifndef WRITER_INCLUDED
#define WRITER_INCLUDED

#include <fstream>
#include <htslib/bgzf.h>

#include "types_and_structs.hpp"
#include "stats_test.hpp"

// Classes to abstract the writing (and reading) of files in stoat. Standard write/read can be done with StdWriter/StdReader,
// while bgzipped content can be written/read using BgzWriter/BgzReader. For Writers, some specialized methods are implemented
// to write specific output. Or the generic *write* function can be used. For Readers, it implements the basic getline function.
// To facilitate parallelization, the actual writing to a file happens in batches. When write() is called, it writes to a buffer,
// which is written to the file (flush()) when it becomes big enough. The buffer is also flushed when the writer is closed.
// Note that if the writer is never closed, then whatever is in the buffer will be lost.

namespace stoat {

// Generic Writer class with the specialized output writer functions
class Writer {
public:
    Writer(const std::string output_file_path, size_t thread_count, size_t max_buffer_length = 1000000);

    std::string get_file_path() const;

    /// Write the given thread to the output file.
    /// Under the hood, this will write the given string to the buffer. 
    /// If this causes the buffer to exceed the maximum length, then flush it and write the whole buffer to the output file.
    bool write(const std::string out_content);


    // flush the output and cleanly close the writer
    void close();

    // header writer
    void write_stoat_output_header(phenotype_type_t phenotype_type);

    // stoat output writers
    void write_binary(const stoat::snarl_info_t& snarl_data, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles);

    void write_binary_covar(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles);

    void write_quantitative(const snarl_info_t& snarl_data, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles);

    void write_eqtl(const snarl_info_t& snarl_data, const std::string& gene_name, const stoat::test_result_t& test_result, const std::vector<bool>& active_alleles);
    
protected:

    //////////////////////// Helper functions
    // Write the buffer to the output file and clear the buffer
    bool flush();

    // The virtual function to write the buffer to a file
    virtual bool write_string_to_file(const std::string& out_content) = 0;
    // The virtual function to close the file
    virtual void close_file() = 0;



    const std::string file_path;
    std::string buffer;
    // Once a buffer exceeds this length, flush it
    size_t max_buffer_length;
    
};

// Standart Writer: uncompressed using a standard output file stream
class StdWriter: public Writer {
public:
    StdWriter(const std::string output_file_path, size_t max_buffer_length = 1000000);

protected:
    bool write_string_to_file(const std::string& out_content);

    void close_file();

    std::ofstream file_stream;
    
};

// Bgzipped writer using a HTSlib
class BgzWriter: public Writer {
public:
    BgzWriter(const std::string output_file_path, size_t max_buffer_length = 1000000);

protected:
    bool write_string_to_file(const std::string& out_content);

    void close_file();

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
