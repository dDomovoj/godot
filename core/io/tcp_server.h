/*************************************************************************/
/*  tcp_server.h                                                         */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2019 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2019 Godot Engine contributors (cf. AUTHORS.md)    */
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

#ifndef TCP_SERVER_H
#define TCP_SERVER_H

#include "core/io/ip.h"
#include "core/io/net_socket.h"
#include "core/io/stream_peer.h"
#include "core/io/stream_peer_tcp.h"

class TCP_Server : public Reference {

	GDCLASS(TCP_Server, Reference);

protected:
	enum {
		MAX_PENDING_CONNECTIONS = 8
	};

	Ref<NetSocket> _sock;
	static void _bind_methods();

public:
	Error listen(uint16_t p_port, const IP_Address &p_bind_address = IP_Address("*"));
	bool is_connection_available() const;
	Ref<StreamPeerTCP> take_connection();

	void stop(); // Stop listening

	TCP_Server();
	~TCP_Server();
};

#endif // TCP_SERVER_H
