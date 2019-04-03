/*************************************************************************/
/*  multiplayer_api.cpp                                                  */
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

#include "multiplayer_api.h"

#include "core/io/marshalls.h"
#include "scene/main/node.h"

_FORCE_INLINE_ bool _should_call_local(MultiplayerAPI::RPCMode mode, bool is_master, bool &r_skip_rpc) {

	switch (mode) {

		case MultiplayerAPI::RPC_MODE_DISABLED: {
			// Do nothing.
		} break;
		case MultiplayerAPI::RPC_MODE_REMOTE: {
			// Do nothing also. Remote cannot produce a local call.
		} break;
		case MultiplayerAPI::RPC_MODE_MASTERSYNC: {
			if (is_master)
				r_skip_rpc = true; // I am the master, so skip remote call.
		} // Do not break, fall over to other sync.
		case MultiplayerAPI::RPC_MODE_REMOTESYNC:
		case MultiplayerAPI::RPC_MODE_PUPPETSYNC: {
			// Call it, sync always results in a local call.
			return true;
		} break;
		case MultiplayerAPI::RPC_MODE_MASTER: {
			if (is_master)
				r_skip_rpc = true; // I am the master, so skip remote call.
			return is_master;
		} break;
		case MultiplayerAPI::RPC_MODE_PUPPET: {
			return !is_master;
		} break;
	}
	return false;
}

_FORCE_INLINE_ bool _can_call_mode(MultiplayerAPI::RPCMode mode, int p_node_master, int p_remote_id) {
	switch (mode) {

		case MultiplayerAPI::RPC_MODE_DISABLED: {
			return false;
		} break;
		case MultiplayerAPI::RPC_MODE_REMOTE:
		case MultiplayerAPI::RPC_MODE_REMOTESYNC: {
			return true;
		} break;
		case MultiplayerAPI::RPC_MODE_MASTERSYNC:
		case MultiplayerAPI::RPC_MODE_MASTER: {
			return p_node_master == NetworkedMultiplayerPeer::TARGET_PEER_SERVER;
		} break;
		case MultiplayerAPI::RPC_MODE_PUPPETSYNC:
		case MultiplayerAPI::RPC_MODE_PUPPET: {
			return p_node_master != NetworkedMultiplayerPeer::TARGET_PEER_SERVER && p_remote_id == p_node_master;
		} break;
	}

	return false;
}

void MultiplayerAPI::poll() {

	if (!network_peer.is_valid() || network_peer->get_connection_status() == NetworkedMultiplayerPeer::CONNECTION_DISCONNECTED)
		return;

	network_peer->poll();

	if (!network_peer.is_valid()) // It's possible that polling might have resulted in a disconnection, so check here.
		return;

	while (network_peer->get_available_packet_count()) {

		int sender = network_peer->get_packet_peer();
		const uint8_t *packet;
		int len;

		Error err = network_peer->get_packet(&packet, len);
		if (err != OK) {
			ERR_PRINT("Error getting packet!");
		}

		rpc_sender_id = sender;
		_process_packet(sender, packet, len);
		rpc_sender_id = 0;

		if (!network_peer.is_valid()) {
			break; // It's also possible that a packet or RPC caused a disconnection, so also check here.
		}
	}
}

void MultiplayerAPI::clear() {
	connected_peers.clear();
	path_get_cache.clear();
	path_send_cache.clear();
	packet_cache.clear();
	last_send_cache_id = 1;
}

void MultiplayerAPI::set_root_node(Node *p_node) {
	root_node = p_node;
}

void MultiplayerAPI::set_network_peer(const Ref<NetworkedMultiplayerPeer> &p_peer) {

	if (p_peer == network_peer) return; // Nothing to do

	if (network_peer.is_valid()) {
		network_peer->disconnect("peer_connected", this, "_add_peer");
		network_peer->disconnect("peer_disconnected", this, "_del_peer");
		network_peer->disconnect("connection_succeeded", this, "_connected_to_server");
		network_peer->disconnect("connection_failed", this, "_connection_failed");
		network_peer->disconnect("server_disconnected", this, "_server_disconnected");
		clear();
	}

	network_peer = p_peer;

	ERR_EXPLAIN("Supplied NetworkedNetworkPeer must be connecting or connected.");
	ERR_FAIL_COND(p_peer.is_valid() && p_peer->get_connection_status() == NetworkedMultiplayerPeer::CONNECTION_DISCONNECTED);

	if (network_peer.is_valid()) {
		network_peer->connect("peer_connected", this, "_add_peer");
		network_peer->connect("peer_disconnected", this, "_del_peer");
		network_peer->connect("connection_succeeded", this, "_connected_to_server");
		network_peer->connect("connection_failed", this, "_connection_failed");
		network_peer->connect("server_disconnected", this, "_server_disconnected");
	}
}

Ref<NetworkedMultiplayerPeer> MultiplayerAPI::get_network_peer() const {
	return network_peer;
}

void MultiplayerAPI::_process_packet(int p_from, const uint8_t *p_packet, int p_packet_len) {

	ERR_EXPLAIN("Multiplayer root node was not initialized. If you are using custom multiplayer, remember to set the root node via MultiplayerAPI.set_root_node before using it");
	ERR_FAIL_COND(root_node == NULL);
	ERR_EXPLAIN("Invalid packet received. Size too small.");
	ERR_FAIL_COND(p_packet_len < 1);

	uint8_t packet_type = p_packet[0];

	switch (packet_type) {

		case NETWORK_COMMAND_SIMPLIFY_PATH: {

			_process_simplify_path(p_from, p_packet, p_packet_len);
		} break;

		case NETWORK_COMMAND_CONFIRM_PATH: {

			_process_confirm_path(p_from, p_packet, p_packet_len);
		} break;

		case NETWORK_COMMAND_REMOTE_CALL:
		case NETWORK_COMMAND_REMOTE_SET: {

			ERR_EXPLAIN("Invalid packet received. Size too small.");
			ERR_FAIL_COND(p_packet_len < 6);

			Node *node = _process_get_node(p_from, p_packet, p_packet_len);

			ERR_EXPLAIN("Invalid packet received. Requested node was not found.");
			ERR_FAIL_COND(node == NULL);

			// Detect cstring end.
			int len_end = 5;
			for (; len_end < p_packet_len; len_end++) {
				if (p_packet[len_end] == 0) {
					break;
				}
			}

			ERR_EXPLAIN("Invalid packet received. Size too small.");
			ERR_FAIL_COND(len_end >= p_packet_len);

			StringName name = String::utf8((const char *)&p_packet[5]);

			if (packet_type == NETWORK_COMMAND_REMOTE_CALL) {

				_process_rpc(node, name, p_from, p_packet, p_packet_len, len_end + 1);

			} else {

				_process_rset(node, name, p_from, p_packet, p_packet_len, len_end + 1);
			}

		} break;

		case NETWORK_COMMAND_RAW: {

			_process_raw(p_from, p_packet, p_packet_len);
		} break;
	}
}

Node *MultiplayerAPI::_process_get_node(int p_from, const uint8_t *p_packet, int p_packet_len) {

	uint32_t target = decode_uint32(&p_packet[1]);
	Node *node = NULL;

	if (target & 0x80000000) {
		// Use full path (not cached yet).

		int ofs = target & 0x7FFFFFFF;

		ERR_EXPLAIN("Invalid packet received. Size smaller than declared.");
		ERR_FAIL_COND_V(ofs >= p_packet_len, NULL);

		String paths;
		paths.parse_utf8((const char *)&p_packet[ofs], p_packet_len - ofs);

		NodePath np = paths;

		node = root_node->get_node(np);

		if (!node)
			ERR_PRINTS("Failed to get path from RPC: " + String(np));
	} else {
		// Use cached path.
		int id = target;

		Map<int, PathGetCache>::Element *E = path_get_cache.find(p_from);
		ERR_EXPLAIN("Invalid packet received. Requests invalid peer cache.");
		ERR_FAIL_COND_V(!E, NULL);

		Map<int, PathGetCache::NodeInfo>::Element *F = E->get().nodes.find(id);
		ERR_EXPLAIN("Invalid packet received. Unabled to find requested cached node.");
		ERR_FAIL_COND_V(!F, NULL);

		PathGetCache::NodeInfo *ni = &F->get();
		// Do proper caching later.

		node = root_node->get_node(ni->path);
		if (!node)
			ERR_PRINTS("Failed to get cached path from RPC: " + String(ni->path));
	}
	return node;
}

void MultiplayerAPI::_process_rpc(Node *p_node, const StringName &p_name, int p_from, const uint8_t *p_packet, int p_packet_len, int p_offset) {

	ERR_EXPLAIN("Invalid packet received. Size too small.");
	ERR_FAIL_COND(p_offset >= p_packet_len);

	// Check that remote can call the RPC on this node.
	RPCMode rpc_mode = RPC_MODE_DISABLED;
	const Map<StringName, RPCMode>::Element *E = p_node->get_node_rpc_mode(p_name);
	if (E) {
		rpc_mode = E->get();
	} else if (p_node->get_script_instance()) {
		rpc_mode = p_node->get_script_instance()->get_rpc_mode(p_name);
	}

	int node_master_id = p_node->get_network_master();
	ERR_EXPLAIN("RPC '" + String(p_name) + "' is not allowed on node " + p_node->get_path() + " from: " + itos(p_from) + ". Mode is " + itos((int)rpc_mode) + ", master is " + itos(node_master_id) + ".");
	ERR_FAIL_COND(!_can_call_mode(rpc_mode, p_from, node_master_id));

	int argc = p_packet[p_offset];
	Vector<Variant> args;
	Vector<const Variant *> argp;
	args.resize(argc);
	argp.resize(argc);

	p_offset++;

	for (int i = 0; i < argc; i++) {

		ERR_EXPLAIN("Invalid packet received. Size too small.");
		ERR_FAIL_COND(p_offset >= p_packet_len);

		int vlen;
		Error err = decode_variant(args.write[i], &p_packet[p_offset], p_packet_len - p_offset, &vlen, allow_object_decoding || network_peer->is_object_decoding_allowed());
		ERR_EXPLAIN("Invalid packet received. Unable to decode RPC argument.");
		ERR_FAIL_COND(err != OK);

		argp.write[i] = &args[i];
		p_offset += vlen;
	}

	Variant::CallError ce;

	p_node->call(p_name, (const Variant **)argp.ptr(), argc, ce);
	if (ce.error != Variant::CallError::CALL_OK) {
		String error = Variant::get_call_error_text(p_node, p_name, (const Variant **)argp.ptr(), argc, ce);
		error = "RPC - " + error;
		ERR_PRINTS(error);
	}
}

void MultiplayerAPI::_process_rset(Node *p_node, const StringName &p_name, int p_from, const uint8_t *p_packet, int p_packet_len, int p_offset) {

	ERR_EXPLAIN("Invalid packet received. Size too small.");
	ERR_FAIL_COND(p_offset >= p_packet_len);

	// Check that remote can call the RSET on this node.
	RPCMode rset_mode = RPC_MODE_DISABLED;
	const Map<StringName, RPCMode>::Element *E = p_node->get_node_rset_mode(p_name);
	if (E) {
		rset_mode = E->get();
	} else if (p_node->get_script_instance()) {
		rset_mode = p_node->get_script_instance()->get_rset_mode(p_name);
	}

	int node_master_id = p_node->get_network_master();
	ERR_EXPLAIN("RSET '" + String(p_name) + "' is not allowed on node " + p_node->get_path() + " from: " + itos(p_from) + ". Mode is " + itos((int)rset_mode) + ", master is " + itos(node_master_id) + ".");
	ERR_FAIL_COND(!_can_call_mode(rset_mode, p_from, node_master_id));

	Variant value;
	Error err = decode_variant(value, &p_packet[p_offset], p_packet_len - p_offset, NULL, allow_object_decoding || network_peer->is_object_decoding_allowed());

	ERR_EXPLAIN("Invalid packet received. Unable to decode RSET value.");
	ERR_FAIL_COND(err != OK);

	bool valid;

	p_node->set(p_name, value, &valid);
	if (!valid) {
		String error = "Error setting remote property '" + String(p_name) + "', not found in object of type " + p_node->get_class();
		ERR_PRINTS(error);
	}
}

void MultiplayerAPI::_process_simplify_path(int p_from, const uint8_t *p_packet, int p_packet_len) {

	ERR_EXPLAIN("Invalid packet received. Size too small.");
	ERR_FAIL_COND(p_packet_len < 5);
	int id = decode_uint32(&p_packet[1]);

	String paths;
	paths.parse_utf8((const char *)&p_packet[5], p_packet_len - 5);

	NodePath path = paths;

	if (!path_get_cache.has(p_from)) {
		path_get_cache[p_from] = PathGetCache();
	}

	PathGetCache::NodeInfo ni;
	ni.path = path;
	ni.instance = 0;

	path_get_cache[p_from].nodes[id] = ni;

	// Encode path to send ack.
	CharString pname = String(path).utf8();
	int len = encode_cstring(pname.get_data(), NULL);

	Vector<uint8_t> packet;

	packet.resize(1 + len);
	packet.write[0] = NETWORK_COMMAND_CONFIRM_PATH;
	encode_cstring(pname.get_data(), &packet.write[1]);

	network_peer->set_transfer_mode(NetworkedMultiplayerPeer::TRANSFER_MODE_RELIABLE);
	network_peer->set_target_peer(p_from);
	network_peer->put_packet(packet.ptr(), packet.size());
}

void MultiplayerAPI::_process_confirm_path(int p_from, const uint8_t *p_packet, int p_packet_len) {

	ERR_EXPLAIN("Invalid packet received. Size too small.");
	ERR_FAIL_COND(p_packet_len < 2);

	String paths;
	paths.parse_utf8((const char *)&p_packet[1], p_packet_len - 1);

	NodePath path = paths;

	PathSentCache *psc = path_send_cache.getptr(path);
	ERR_EXPLAIN("Invalid packet received. Tries to confirm a path which was not found in cache.");
	ERR_FAIL_COND(!psc);

	Map<int, bool>::Element *E = psc->confirmed_peers.find(p_from);
	ERR_EXPLAIN("Invalid packet received. Source peer was not found in cache for the given path.");
	ERR_FAIL_COND(!E);
	E->get() = true;
}

bool MultiplayerAPI::_send_confirm_path(NodePath p_path, PathSentCache *psc, int p_target) {
	bool has_all_peers = true;
	List<int> peers_to_add; // If one is missing, take note to add it.

	for (Set<int>::Element *E = connected_peers.front(); E; E = E->next()) {

		if (p_target < 0 && E->get() == -p_target)
			continue; // Continue, excluded.

		if (p_target > 0 && E->get() != p_target)
			continue; // Continue, not for this peer.

		Map<int, bool>::Element *F = psc->confirmed_peers.find(E->get());

		if (!F || !F->get()) {
			// Path was not cached, or was cached but is unconfirmed.
			if (!F) {
				// Not cached at all, take note.
				peers_to_add.push_back(E->get());
			}

			has_all_peers = false;
		}
	}

	// Those that need to be added, send a message for this.

	for (List<int>::Element *E = peers_to_add.front(); E; E = E->next()) {

		// Encode function name.
		CharString pname = String(p_path).utf8();
		int len = encode_cstring(pname.get_data(), NULL);

		Vector<uint8_t> packet;

		packet.resize(1 + 4 + len);
		packet.write[0] = NETWORK_COMMAND_SIMPLIFY_PATH;
		encode_uint32(psc->id, &packet.write[1]);
		encode_cstring(pname.get_data(), &packet.write[5]);

		network_peer->set_target_peer(E->get()); // To all of you.
		network_peer->set_transfer_mode(NetworkedMultiplayerPeer::TRANSFER_MODE_RELIABLE);
		network_peer->put_packet(packet.ptr(), packet.size());

		psc->confirmed_peers.insert(E->get(), false); // Insert into confirmed, but as false since it was not confirmed.
	}

	return has_all_peers;
}

void MultiplayerAPI::_send_rpc(Node *p_from, int p_to, bool p_unreliable, bool p_set, const StringName &p_name, const Variant **p_arg, int p_argcount) {

	if (network_peer.is_null()) {
		ERR_EXPLAIN("Attempt to remote call/set when networking is not active in SceneTree.");
		ERR_FAIL();
	}

	if (network_peer->get_connection_status() == NetworkedMultiplayerPeer::CONNECTION_CONNECTING) {
		ERR_EXPLAIN("Attempt to remote call/set when networking is not connected yet in SceneTree.");
		ERR_FAIL();
	}

	if (network_peer->get_connection_status() == NetworkedMultiplayerPeer::CONNECTION_DISCONNECTED) {
		ERR_EXPLAIN("Attempt to remote call/set when networking is disconnected.");
		ERR_FAIL();
	}

	if (p_argcount > 255) {
		ERR_EXPLAIN("Too many arguments >255.");
		ERR_FAIL();
	}

	if (p_to != 0 && !connected_peers.has(ABS(p_to))) {
		if (p_to == network_peer->get_unique_id()) {
			ERR_EXPLAIN("Attempt to remote call/set yourself! unique ID: " + itos(network_peer->get_unique_id()));
		} else {
			ERR_EXPLAIN("Attempt to remote call unexisting ID: " + itos(p_to));
		}

		ERR_FAIL();
	}

	NodePath from_path = (root_node->get_path()).rel_path_to(p_from->get_path());
	ERR_EXPLAIN("Unable to send RPC. Relative path is empty. THIS IS LIKELY A BUG IN THE ENGINE!");
	ERR_FAIL_COND(from_path.is_empty());

	// See if the path is cached.
	PathSentCache *psc = path_send_cache.getptr(from_path);
	if (!psc) {
		// Path is not cached, create.
		path_send_cache[from_path] = PathSentCache();
		psc = path_send_cache.getptr(from_path);
		psc->id = last_send_cache_id++;
	}

	// Create base packet, lots of hardcode because it must be tight.

	int ofs = 0;

#define MAKE_ROOM(m_amount) \
	if (packet_cache.size() < m_amount) packet_cache.resize(m_amount);

	// Encode type.
	MAKE_ROOM(1);
	packet_cache.write[0] = p_set ? NETWORK_COMMAND_REMOTE_SET : NETWORK_COMMAND_REMOTE_CALL;
	ofs += 1;

	// Encode ID.
	MAKE_ROOM(ofs + 4);
	encode_uint32(psc->id, &(packet_cache.write[ofs]));
	ofs += 4;

	// Encode function name.
	CharString name = String(p_name).utf8();
	int len = encode_cstring(name.get_data(), NULL);
	MAKE_ROOM(ofs + len);
	encode_cstring(name.get_data(), &(packet_cache.write[ofs]));
	ofs += len;

	if (p_set) {
		// Set argument.
		Error err = encode_variant(*p_arg[0], NULL, len, allow_object_decoding || network_peer->is_object_decoding_allowed());
		ERR_EXPLAIN("Unable to encode RSET value. THIS IS LIKELY A BUG IN THE ENGINE!");
		ERR_FAIL_COND(err != OK);
		MAKE_ROOM(ofs + len);
		encode_variant(*p_arg[0], &(packet_cache.write[ofs]), len, allow_object_decoding || network_peer->is_object_decoding_allowed());
		ofs += len;

	} else {
		// Call arguments.
		MAKE_ROOM(ofs + 1);
		packet_cache.write[ofs] = p_argcount;
		ofs += 1;
		for (int i = 0; i < p_argcount; i++) {
			Error err = encode_variant(*p_arg[i], NULL, len, allow_object_decoding || network_peer->is_object_decoding_allowed());
			ERR_EXPLAIN("Unable to encode RPC argument. THIS IS LIKELY A BUG IN THE ENGINE!");
			ERR_FAIL_COND(err != OK);
			MAKE_ROOM(ofs + len);
			encode_variant(*p_arg[i], &(packet_cache.write[ofs]), len, allow_object_decoding || network_peer->is_object_decoding_allowed());
			ofs += len;
		}
	}

	// See if all peers have cached path (is so, call can be fast).
	bool has_all_peers = _send_confirm_path(from_path, psc, p_to);

	// Take chance and set transfer mode, since all send methods will use it.
	network_peer->set_transfer_mode(p_unreliable ? NetworkedMultiplayerPeer::TRANSFER_MODE_UNRELIABLE : NetworkedMultiplayerPeer::TRANSFER_MODE_RELIABLE);

	if (has_all_peers) {

		// They all have verified paths, so send fast.
		network_peer->set_target_peer(p_to); // To all of you.
		network_peer->put_packet(packet_cache.ptr(), ofs); // A message with love.
	} else {
		// Not all verified path, so send one by one.

		// Append path at the end, since we will need it for some packets.
		CharString pname = String(from_path).utf8();
		int path_len = encode_cstring(pname.get_data(), NULL);
		MAKE_ROOM(ofs + path_len);
		encode_cstring(pname.get_data(), &(packet_cache.write[ofs]));

		for (Set<int>::Element *E = connected_peers.front(); E; E = E->next()) {

			if (p_to < 0 && E->get() == -p_to)
				continue; // Continue, excluded.

			if (p_to > 0 && E->get() != p_to)
				continue; // Continue, not for this peer.

			Map<int, bool>::Element *F = psc->confirmed_peers.find(E->get());
			ERR_CONTINUE(!F); // Should never happen.

			network_peer->set_target_peer(E->get()); // To this one specifically.

			if (F->get()) {
				// This one confirmed path, so use id.
				encode_uint32(psc->id, &(packet_cache.write[1]));
				network_peer->put_packet(packet_cache.ptr(), ofs);
			} else {
				// This one did not confirm path yet, so use entire path (sorry!).
				encode_uint32(0x80000000 | ofs, &(packet_cache.write[1])); // Offset to path and flag.
				network_peer->put_packet(packet_cache.ptr(), ofs + path_len);
			}
		}
	}
}

void MultiplayerAPI::_add_peer(int p_id) {
	connected_peers.insert(p_id);
	path_get_cache.insert(p_id, PathGetCache());
	emit_signal("network_peer_connected", p_id);
}

void MultiplayerAPI::_del_peer(int p_id) {
	connected_peers.erase(p_id);
	path_get_cache.erase(p_id); // I no longer need your cache, sorry.
	emit_signal("network_peer_disconnected", p_id);
}

void MultiplayerAPI::_connected_to_server() {

	emit_signal("connected_to_server");
}

void MultiplayerAPI::_connection_failed() {

	emit_signal("connection_failed");
}

void MultiplayerAPI::_server_disconnected() {

	emit_signal("server_disconnected");
}

void MultiplayerAPI::rpcp(Node *p_node, int p_peer_id, bool p_unreliable, const StringName &p_method, const Variant **p_arg, int p_argcount) {

	ERR_EXPLAIN("Trying to call an RPC while no network peer is active.");
	ERR_FAIL_COND(!network_peer.is_valid());
	ERR_EXPLAIN("Trying to call an RPC on a node which is not inside SceneTree.");
	ERR_FAIL_COND(!p_node->is_inside_tree());
	ERR_EXPLAIN("Trying to call an RPC via a network peer which is not connected.");
	ERR_FAIL_COND(network_peer->get_connection_status() != NetworkedMultiplayerPeer::CONNECTION_CONNECTED);

	int node_id = network_peer->get_unique_id();
	bool skip_rpc = false;
	bool call_local_native = false;
	bool call_local_script = false;
	bool is_master = p_node->is_network_master();

	if (p_peer_id == 0 || p_peer_id == node_id || (p_peer_id < 0 && p_peer_id != -node_id)) {
		// Check that send mode can use local call.

		const Map<StringName, RPCMode>::Element *E = p_node->get_node_rpc_mode(p_method);
		if (E) {
			call_local_native = _should_call_local(E->get(), is_master, skip_rpc);
		}

		if (call_local_native) {
			// Done below.
		} else if (p_node->get_script_instance()) {
			// Attempt with script.
			RPCMode rpc_mode = p_node->get_script_instance()->get_rpc_mode(p_method);
			call_local_script = _should_call_local(rpc_mode, is_master, skip_rpc);
		}
	}

	if (!skip_rpc) {
		_send_rpc(p_node, p_peer_id, p_unreliable, false, p_method, p_arg, p_argcount);
	}

	if (call_local_native) {
		Variant::CallError ce;
		p_node->call(p_method, p_arg, p_argcount, ce);
		if (ce.error != Variant::CallError::CALL_OK) {
			String error = Variant::get_call_error_text(p_node, p_method, p_arg, p_argcount, ce);
			error = "rpc() aborted in local call:  - " + error;
			ERR_PRINTS(error);
			return;
		}
	}

	if (call_local_script) {
		Variant::CallError ce;
		ce.error = Variant::CallError::CALL_OK;
		p_node->get_script_instance()->call(p_method, p_arg, p_argcount, ce);
		if (ce.error != Variant::CallError::CALL_OK) {
			String error = Variant::get_call_error_text(p_node, p_method, p_arg, p_argcount, ce);
			error = "rpc() aborted in script local call:  - " + error;
			ERR_PRINTS(error);
			return;
		}
	}
}

void MultiplayerAPI::rsetp(Node *p_node, int p_peer_id, bool p_unreliable, const StringName &p_property, const Variant &p_value) {

	ERR_EXPLAIN("Trying to RSET while no network peer is active.");
	ERR_FAIL_COND(!network_peer.is_valid());
	ERR_EXPLAIN("Trying to RSET on a node which is not inside SceneTree.");
	ERR_FAIL_COND(!p_node->is_inside_tree());
	ERR_EXPLAIN("Trying to send an RSET via a network peer which is not connected.");
	ERR_FAIL_COND(network_peer->get_connection_status() != NetworkedMultiplayerPeer::CONNECTION_CONNECTED);

	int node_id = network_peer->get_unique_id();
	bool is_master = p_node->is_network_master();
	bool skip_rset = false;

	if (p_peer_id == 0 || p_peer_id == node_id || (p_peer_id < 0 && p_peer_id != -node_id)) {
		// Check that send mode can use local call.

		bool set_local = false;

		const Map<StringName, RPCMode>::Element *E = p_node->get_node_rset_mode(p_property);
		if (E) {

			set_local = _should_call_local(E->get(), is_master, skip_rset);
		}

		if (set_local) {
			bool valid;
			p_node->set(p_property, p_value, &valid);

			if (!valid) {
				String error = "rset() aborted in local set, property not found:  - " + String(p_property);
				ERR_PRINTS(error);
				return;
			}
		} else if (p_node->get_script_instance()) {
			// Attempt with script.
			RPCMode rpc_mode = p_node->get_script_instance()->get_rset_mode(p_property);

			set_local = _should_call_local(rpc_mode, is_master, skip_rset);

			if (set_local) {

				bool valid = p_node->get_script_instance()->set(p_property, p_value);

				if (!valid) {
					String error = "rset() aborted in local script set, property not found:  - " + String(p_property);
					ERR_PRINTS(error);
					return;
				}
			}
		}
	}

	if (skip_rset)
		return;

	const Variant *vptr = &p_value;

	_send_rpc(p_node, p_peer_id, p_unreliable, true, p_property, &vptr, 1);
}

Error MultiplayerAPI::send_bytes(PoolVector<uint8_t> p_data, int p_to, NetworkedMultiplayerPeer::TransferMode p_mode) {

	ERR_EXPLAIN("Trying to send an empty raw packet.");
	ERR_FAIL_COND_V(p_data.size() < 1, ERR_INVALID_DATA);
	ERR_EXPLAIN("Trying to send a raw packet while no network peer is active.");
	ERR_FAIL_COND_V(!network_peer.is_valid(), ERR_UNCONFIGURED);
	ERR_EXPLAIN("Trying to send a raw packet via a network peer which is not connected.");
	ERR_FAIL_COND_V(network_peer->get_connection_status() != NetworkedMultiplayerPeer::CONNECTION_CONNECTED, ERR_UNCONFIGURED);

	MAKE_ROOM(p_data.size() + 1);
	PoolVector<uint8_t>::Read r = p_data.read();
	packet_cache.write[0] = NETWORK_COMMAND_RAW;
	memcpy(&packet_cache.write[1], &r[0], p_data.size());

	network_peer->set_target_peer(p_to);
	network_peer->set_transfer_mode(p_mode);

	return network_peer->put_packet(packet_cache.ptr(), p_data.size() + 1);
}

void MultiplayerAPI::_process_raw(int p_from, const uint8_t *p_packet, int p_packet_len) {

	ERR_EXPLAIN("Invalid packet received. Size too small.");
	ERR_FAIL_COND(p_packet_len < 2);

	PoolVector<uint8_t> out;
	int len = p_packet_len - 1;
	out.resize(len);
	{
		PoolVector<uint8_t>::Write w = out.write();
		memcpy(&w[0], &p_packet[1], len);
	}
	emit_signal("network_peer_packet", p_from, out);
}

int MultiplayerAPI::get_network_unique_id() const {

	ERR_EXPLAIN("No network peer is assigned. Unable to get unique network ID.");
	ERR_FAIL_COND_V(!network_peer.is_valid(), 0);
	return network_peer->get_unique_id();
}

bool MultiplayerAPI::is_network_server() const {

	// XXX Maybe fail silently? Maybe should actually return true to make development of both local and online multiplayer easier?
	ERR_EXPLAIN("No network peer is assigned. I can't be a server.");
	ERR_FAIL_COND_V(!network_peer.is_valid(), false);
	return network_peer->is_server();
}

void MultiplayerAPI::set_refuse_new_network_connections(bool p_refuse) {

	ERR_EXPLAIN("No network peer is assigned. Unable to set 'refuse_new_connections'.");
	ERR_FAIL_COND(!network_peer.is_valid());
	network_peer->set_refuse_new_connections(p_refuse);
}

bool MultiplayerAPI::is_refusing_new_network_connections() const {

	ERR_EXPLAIN("No network peer is assigned. Unable to get 'refuse_new_connections'.");
	ERR_FAIL_COND_V(!network_peer.is_valid(), false);
	return network_peer->is_refusing_new_connections();
}

Vector<int> MultiplayerAPI::get_network_connected_peers() const {

	ERR_EXPLAIN("No network peer is assigned. Assume no peers are connected.");
	ERR_FAIL_COND_V(!network_peer.is_valid(), Vector<int>());

	Vector<int> ret;
	for (Set<int>::Element *E = connected_peers.front(); E; E = E->next()) {
		ret.push_back(E->get());
	}

	return ret;
}

void MultiplayerAPI::set_allow_object_decoding(bool p_enable) {

	allow_object_decoding = p_enable;
}

bool MultiplayerAPI::is_object_decoding_allowed() const {

	return allow_object_decoding;
}

void MultiplayerAPI::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_root_node", "node"), &MultiplayerAPI::set_root_node);
	ClassDB::bind_method(D_METHOD("send_bytes", "bytes", "id", "mode"), &MultiplayerAPI::send_bytes, DEFVAL(NetworkedMultiplayerPeer::TARGET_PEER_BROADCAST), DEFVAL(NetworkedMultiplayerPeer::TRANSFER_MODE_RELIABLE));
	ClassDB::bind_method(D_METHOD("has_network_peer"), &MultiplayerAPI::has_network_peer);
	ClassDB::bind_method(D_METHOD("get_network_peer"), &MultiplayerAPI::get_network_peer);
	ClassDB::bind_method(D_METHOD("get_network_unique_id"), &MultiplayerAPI::get_network_unique_id);
	ClassDB::bind_method(D_METHOD("is_network_server"), &MultiplayerAPI::is_network_server);
	ClassDB::bind_method(D_METHOD("get_rpc_sender_id"), &MultiplayerAPI::get_rpc_sender_id);
	ClassDB::bind_method(D_METHOD("_add_peer", "id"), &MultiplayerAPI::_add_peer);
	ClassDB::bind_method(D_METHOD("_del_peer", "id"), &MultiplayerAPI::_del_peer);
	ClassDB::bind_method(D_METHOD("set_network_peer", "peer"), &MultiplayerAPI::set_network_peer);
	ClassDB::bind_method(D_METHOD("poll"), &MultiplayerAPI::poll);
	ClassDB::bind_method(D_METHOD("clear"), &MultiplayerAPI::clear);

	ClassDB::bind_method(D_METHOD("_connected_to_server"), &MultiplayerAPI::_connected_to_server);
	ClassDB::bind_method(D_METHOD("_connection_failed"), &MultiplayerAPI::_connection_failed);
	ClassDB::bind_method(D_METHOD("_server_disconnected"), &MultiplayerAPI::_server_disconnected);
	ClassDB::bind_method(D_METHOD("get_network_connected_peers"), &MultiplayerAPI::get_network_connected_peers);
	ClassDB::bind_method(D_METHOD("set_refuse_new_network_connections", "refuse"), &MultiplayerAPI::set_refuse_new_network_connections);
	ClassDB::bind_method(D_METHOD("is_refusing_new_network_connections"), &MultiplayerAPI::is_refusing_new_network_connections);
	ClassDB::bind_method(D_METHOD("set_allow_object_decoding", "enable"), &MultiplayerAPI::set_allow_object_decoding);
	ClassDB::bind_method(D_METHOD("is_object_decoding_allowed"), &MultiplayerAPI::is_object_decoding_allowed);

	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "allow_object_decoding"), "set_allow_object_decoding", "is_object_decoding_allowed");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "refuse_new_network_connections"), "set_refuse_new_network_connections", "is_refusing_new_network_connections");
	ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "network_peer", PROPERTY_HINT_RESOURCE_TYPE, "NetworkedMultiplayerPeer", 0), "set_network_peer", "get_network_peer");

	ADD_SIGNAL(MethodInfo("network_peer_connected", PropertyInfo(Variant::INT, "id")));
	ADD_SIGNAL(MethodInfo("network_peer_disconnected", PropertyInfo(Variant::INT, "id")));
	ADD_SIGNAL(MethodInfo("network_peer_packet", PropertyInfo(Variant::INT, "id"), PropertyInfo(Variant::POOL_BYTE_ARRAY, "packet")));
	ADD_SIGNAL(MethodInfo("connected_to_server"));
	ADD_SIGNAL(MethodInfo("connection_failed"));
	ADD_SIGNAL(MethodInfo("server_disconnected"));

	BIND_ENUM_CONSTANT(RPC_MODE_DISABLED);
	BIND_ENUM_CONSTANT(RPC_MODE_REMOTE);
	BIND_ENUM_CONSTANT(RPC_MODE_MASTER);
	BIND_ENUM_CONSTANT(RPC_MODE_PUPPET);
	BIND_ENUM_CONSTANT(RPC_MODE_SLAVE); // Deprecated.
	BIND_ENUM_CONSTANT(RPC_MODE_REMOTESYNC);
	BIND_ENUM_CONSTANT(RPC_MODE_SYNC); // Deprecated.
	BIND_ENUM_CONSTANT(RPC_MODE_MASTERSYNC);
	BIND_ENUM_CONSTANT(RPC_MODE_PUPPETSYNC);
}

MultiplayerAPI::MultiplayerAPI() :
		allow_object_decoding(false) {
	rpc_sender_id = 0;
	root_node = NULL;
	clear();
}

MultiplayerAPI::~MultiplayerAPI() {
	clear();
}
