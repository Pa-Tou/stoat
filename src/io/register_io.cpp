/**
 * \file register_io.hpp
 * Includes calls to register all libvg types with libvgio.
 * Copied from vg
 */

// Keep these includes in alphabetical order.

#include "register_loader_saver_gbz.hpp"
#include "register_loader_saver_gbzgraph.hpp"
#include "register_loader_saver_packed_graph.hpp"
#include "register_loader_saver_hash_graph.hpp"

#include "register_io.hpp"

namespace stoat {

namespace io {

bool register_libvg_io() {
    register_loader_saver_gbz();
    register_loader_saver_gbzgraph();
    register_loader_saver_packed_graph();
    register_loader_saver_hash_graph();
    return true;
}
    
} // end io

} // end stoat
