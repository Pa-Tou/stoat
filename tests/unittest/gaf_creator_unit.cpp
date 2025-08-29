#include <catch.hpp>
#include "../../src/gaf_creator.hpp"
#include <fstream>
#include <sstream>

using namespace stoat_vcf;

TEST_CASE("Calcul des proportions significatives", "[calcul_proportion_signi]") {
    SECTION("Cas normal avec groupes égaux") {
        auto result = stoat_vcf::calcul_proportion_signi(10, 10, 0.01);
        REQUIRE(result.first >= 0.0);
        REQUIRE(result.first <= 60.0);
        REQUIRE(result.second >= 0.0);
        REQUIRE(result.second <= 60.0);
        REQUIRE((result.first + result.second) <= 60.0);
    }

    SECTION("Cas avec groupes de tailles différentes") {
        auto result = stoat_vcf::calcul_proportion_signi(15, 5, 0.01);
        REQUIRE(result.first >= 0.0);
        REQUIRE(result.first <= 60.0);
        REQUIRE(result.second >= 0.0);
        REQUIRE(result.second <= 60.0);
        REQUIRE((result.first + result.second) <= 60.0);
    }

    SECTION("Cas limite avec petits groupes") {
        auto result = stoat_vcf::calcul_proportion_signi(1, 1, 0.01);
        REQUIRE(result.first >= 0.0);
        REQUIRE(result.first <= 60.0);
        REQUIRE(result.second >= 0.0);
        REQUIRE(result.second <= 60.0);
        REQUIRE((result.first + result.second) <= 60.0);
    }

    SECTION("Cas avec groupe vide") {
        auto result = stoat_vcf::calcul_proportion_signi(0, 0, 0.01);
        REQUIRE(result.first == 0.0);
        REQUIRE(result.second == 0.0);
    }

    SECTION("Cas avec p-value extrême") {
        auto result = stoat_vcf::calcul_proportion_signi(10, 10, 1e-10);
        REQUIRE(result.first >= 0.0);
        REQUIRE(result.first <= 60.0);
        REQUIRE(result.second >= 0.0);
        REQUIRE(result.second <= 60.0);
        REQUIRE((result.first + result.second) <= 60.0);
    }
}

TEST_CASE("Ajout de suffixe au nom de fichier", "[stoat_vcf::addSuffixToFilename]") {
    SECTION("Ajout de suffixe à un nom simple") {
        REQUIRE(stoat_vcf::addSuffixToFilename("test.txt", "_suffix") == "test_suffix.txt");
    }

    SECTION("Ajout de suffixe à un nom avec chemin") {
        REQUIRE(stoat_vcf::addSuffixToFilename("/path/to/test.txt", "_suffix") == "/path/to/test_suffix.txt");
    }

    SECTION("Ajout de suffixe à un nom sans extension") {
        REQUIRE(stoat_vcf::addSuffixToFilename("test", "_suffix") == "test_suffix");
    }
}

TEST_CASE("Décomposition de snarl", "[stoat_vcf::decompose_snarl]") {
    SECTION("Snarl simple") {
        std::string snarl = "1>2>3";
        auto result = stoat_vcf::decompose_snarl(snarl);
        REQUIRE(result.size() == 3);
        REQUIRE(result[0] == 1);
        REQUIRE(result[1] == 2);
        REQUIRE(result[2] == 3);
    }

    SECTION("Snarl vide") {
        std::string snarl = "";
        auto result = stoat_vcf::decompose_snarl(snarl);
        REQUIRE(result.empty());
    }

    SECTION("Snarl avec un seul nœud") {
        std::string snarl = "42";
        auto result = stoat_vcf::decompose_snarl(snarl);
        REQUIRE(result.size() == 1);
        REQUIRE(result[0] == 42);
    }
}

TEST_CASE("Écriture des lignes GAF", "[stoat_vcf::write_gaf_lines]") {
    SECTION("Test d'écriture de base") {
        std::ofstream outfile("test_gaf.txt");
        stoat_vcf::write_gaf_lines("seq1", "1>2>3", 100, 0.75, outfile);
        outfile.close();

        std::ifstream infile("test_gaf.txt");
        std::string line;
        REQUIRE(std::getline(infile, line));
        REQUIRE(line.find("seq1") != std::string::npos);
        REQUIRE(line.find("1>2>3") != std::string::npos);
        infile.close();
        std::remove("test_gaf.txt");
    }
}
