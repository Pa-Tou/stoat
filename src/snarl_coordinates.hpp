#ifndef STOAT_SNARL_COORDINATES_HPP_INCLUDED
#define STOAT_SNARL_COORDINATES_HPP_INCLUDED

#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "types_and_structs.hpp"

namespace stoat {

/// A class for holding the reference coordinates of snarls
/// This will first try to find coordinates on the given sample. If the given sample does not pass through the snarl, then report the coordinates of the lowest
/// ancestor of the snarl with reference coordinates. 
/// If none of the ancestors have reference coordinates, then try reference-sense paths then any path.
/// To use specific references, call add_reference_path() with the full path name. For consistency, add_reference_path() should be called all at once before looking for coordinates
/// Since a path may traverse a snarl multiple times, report the largest range of coordinates of the given path
class SnarlCoordinates {

    //////////////////////////////////////////////////////// Public interface
    public:
        SnarlCoordinates() {};

        /// Tell the SnarlCoordinates class about a reference, as its path name (the full name which includes the sample name, chr name, etc)
        /// This will remove duplicates
        void add_reference_path(std::string path_name);

        /// Return the reference coordinates as a tuple of reference name, start offset, and end offset
        /// The offsets do not include the boundary nodes 
        /// This firs tries to find a reference path given by add_reference_path, then a reference path given by add_reference_path on an ancestor snarl,
        /// then a reference-sense path on this snarl, then a reference-sense path on an ancestor snarl, then any path on this snarl, then any path on an
        /// ancestor snarl
        /// Returns "NA", max(), max() if no coordinates were found
        std::tuple<std::string, size_t, size_t> get_reference_coordinates_as_string(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                    const bdsg::SnarlDistanceIndex& distance_index,
                                                                                    net_handle_t snarl);

        /// As above, but instead of returning the reference name, return an index into reference_names_as_vector() 
        /// Returns max(), max(), max() if no coordinates were found
        std::tuple<size_t, size_t, size_t> get_reference_coordinates_as_index(const handlegraph::PathPositionHandleGraph& graph,
                                                                              const bdsg::SnarlDistanceIndex& distance_index,
                                                                              net_handle_t snarl);

        /// Get a vector of reference path names, ordered by the order in which the references added in add_reference(), then by the order in which new paths were found
        std::vector<std::string> reference_names_as_vector() const;

        /// Get the path name as a string
        std::string get_path_name_from_index(size_t i) const;

        /// Have we seen this reference before?
        bool has_reference(const std::string& pathname) const { return path_to_index.count(pathname) != 0; }

        /// Clear any previously found data, including previously found paths
        void clear();


    //////////////////////////////////////////////////// Private data members
    private:
        enum PathType {REF, REF_SENSE, OTHER};

        // Keep track of which reference paths we have and an index for them
        // The bool indicates whether it was a reference. This is so that we can prioritize references
        std::unordered_map<std::string, std::pair<size_t, PathType>> path_to_index;

        // The inverse of path_to_index. This is duplicative but it lets us store the paths as indices in the snarl to coordinates map
        std::vector<std::pair<std::string, PathType>> path_by_index;

        // Map a snarl (as a start_end traversal of the snarl net handle) to its coordinates, as a tuple of <reference index, start offset, end offset> 
        std::unordered_map<net_handle_t, std::tuple<size_t, size_t, size_t>> snarl_to_coordinates; 


    /////////////////////////////////////////////////// Private helper functions

    private:
    /// Get a list of all traversals (which may or may not traverse both boundary nodes) as path handle, start coordinate, end coordinate
    std::vector<std::tuple<handlegraph::path_handle_t, size_t, size_t>> get_traversals_of_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                                const bdsg::SnarlDistanceIndex& distance_index,
                                                                                                net_handle_t snarl);

};
}

#endif
