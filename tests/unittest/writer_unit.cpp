#include <catch.hpp>
#include <filesystem>
#include <fstream>
#include "../../src/writer.hpp"  // adjust path as needed
#include "../../src/log.hpp"  // adjust path as needed

TEST_CASE("Write trivial file", "[writer]") {
    SECTION("StdWriter") {
        StdWriter writer("temp.txt", 5);
        writer.write("ABC\n");
        writer.write("DEFGH\n");
        writer.write("I");
        writer.close();

        SECTION("StdReader") {
            StdReader reader("temp.txt");
            std::string line
            reader.getline(line);
            REQUIRE(line == "ABC");
            reader.getline(line);
            REQUIRE(line == "DEFGH");
            reader.getline(line);
            REQUIRE(line == "I");
            REQUIRE(!reader.getline(line));
            reader.close();
        }
    }
    SECTION("BgzWriter") {
        BgzWriter writer("temp.txt", 5);
        writer.write("ABC\n");
        writer.write("DEFGH\n");
        writer.write("I");
        writer.close();

        SECTION("BgzReader") {
            BgzReader reader("temp.txt");
            std::string line
            reader.getline(line);
            REQUIRE(line == "ABC");
            reader.getline(line);
            REQUIRE(line == "DEFGH");
            reader.getline(line);
            REQUIRE(line == "I");
            REQUIRE(!reader.getline(line));
            reader.close();
        }
    }
}
