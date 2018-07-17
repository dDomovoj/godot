/*************************************************************************/
/*  stream_peer_ssl.cpp                                                  */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2018 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2018 Godot Engine contributors (cf. AUTHORS.md)    */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "stream_peer_ssl.h"
#include "os/file_access.h"
#include "project_settings.h"

StreamPeerSSL *(*StreamPeerSSL::_create)() = NULL;

StreamPeerSSL *StreamPeerSSL::create() {

	return _create();
}

StreamPeerSSL::LoadCertsFromMemory StreamPeerSSL::load_certs_func = NULL;
bool StreamPeerSSL::available = false;
bool StreamPeerSSL::initialize_certs = true;

void StreamPeerSSL::load_certs_from_memory(const PoolByteArray &p_memory) {
	if (load_certs_func)
		load_certs_func(p_memory);
}

bool StreamPeerSSL::is_available() {
	return available;
}

void StreamPeerSSL::set_blocking_handshake_enabled(bool p_enabled) {
	blocking_handshake = p_enabled;
}

bool StreamPeerSSL::is_blocking_handshake_enabled() const {
	return blocking_handshake;
}

PoolByteArray StreamPeerSSL::get_project_cert_array() {

	PoolByteArray out;
	String certs_path = GLOBAL_DEF("network/ssl/certificates", "");
	ProjectSettings::get_singleton()->set_custom_property_info("network/ssl/certificates", PropertyInfo(Variant::STRING, "network/ssl/certificates", PROPERTY_HINT_FILE, "*.crt"));

	if (certs_path != "") {

		FileAccess *f = FileAccess::open(certs_path, FileAccess::READ);
		if (f) {
			int flen = f->get_len();
			out.resize(flen + 1);
			{
				PoolByteArray::Write w = out.write();
				f->get_buffer(w.ptr(), flen);
				w[flen] = 0; //end f string
			}

			memdelete(f);

#ifdef DEBUG_ENABLED
			print_line("Loaded certs from '" + certs_path);
#endif
		}
	}

	return out;
}

void StreamPeerSSL::_bind_methods() {

	ClassDB::bind_method(D_METHOD("poll"), &StreamPeerSSL::poll);
	ClassDB::bind_method(D_METHOD("accept_stream"), &StreamPeerSSL::accept_stream);
	ClassDB::bind_method(D_METHOD("connect_to_stream", "stream", "validate_certs", "for_hostname"), &StreamPeerSSL::connect_to_stream, DEFVAL(false), DEFVAL(String()));
	ClassDB::bind_method(D_METHOD("get_status"), &StreamPeerSSL::get_status);
	ClassDB::bind_method(D_METHOD("disconnect_from_stream"), &StreamPeerSSL::disconnect_from_stream);
	ClassDB::bind_method(D_METHOD("set_blocking_handshake_enabled", "enabled"), &StreamPeerSSL::set_blocking_handshake_enabled);
	ClassDB::bind_method(D_METHOD("is_blocking_handshake_enabled"), &StreamPeerSSL::is_blocking_handshake_enabled);

	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "blocking_handshake"), "set_blocking_handshake_enabled", "is_blocking_handshake_enabled");

	BIND_ENUM_CONSTANT(STATUS_DISCONNECTED);
	BIND_ENUM_CONSTANT(STATUS_CONNECTED);
	BIND_ENUM_CONSTANT(STATUS_ERROR);
	BIND_ENUM_CONSTANT(STATUS_ERROR_HOSTNAME_MISMATCH);
}

StreamPeerSSL::StreamPeerSSL() {
	blocking_handshake = true;
}
