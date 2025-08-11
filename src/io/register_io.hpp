#ifndef STOAT_IO_REGISTER_IO_HPP_INCLUDED
#define STOAT_IO_REGISTER_IO_HPP_INCLUDED

/**
 * \file register_io.hpp
 * Includes a function to call to register IO handlers for vg types.
 * Copied from vg
 */

namespace stoat {

namespace io {

/**
 * Register libvg types with libvgio.
 * Must be called by library users before doing IO.
 * Does not magically statically call itself.
 * Returns true on success.
 */
bool register_libvg_io();

}

}

#endif
