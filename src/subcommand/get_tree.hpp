// This file is part of STOAT 0.0.1, copyright (C) 2024-2025 
// Authors : Matis Alias-Bagarre, Xian Hui Chang & Jean Monlong.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "../log.hpp"
#include <bdsg/snarl_distance_index.hpp>
#include <handlegraph/path_position_handle_graph.hpp>

namespace stoat_command {

void print_help_get_tree();

int main_stoat_get_tree(int argc, char* argv[]);

void print_tree_tsv(const bdsg::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, std::string& ref_path, size_t start, size_t end);
} // end namespace
