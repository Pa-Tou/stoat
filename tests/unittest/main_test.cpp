#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <omp.h>

int main(int argc, char* argv[]) {
    // Existing parser tests assert record order; production callbacks are unordered.
    omp_set_num_threads(1);
    return Catch::Session().run(argc, argv);
}
