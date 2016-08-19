#ifndef NETWORKED_MULTIPLAYER_ENET_H
#define NETWORKED_MULTIPLAYER_ENET_H

#include "io/networked_multiplayer_peer.h"
#include "enet/enet.h"

class NetworkedMultiplayerENet : public NetworkedMultiplayerPeer {

	OBJ_TYPE(NetworkedMultiplayerENet,NetworkedMultiplayerPeer)

	enum {
		SYSMSG_ADD_PEER,
		SYSMSG_REMOVE_PEER
	};

	bool active;
	bool server;

	uint32_t unique_id;

	int target_peer;
	TransferMode transfer_mode;

	ENetEvent event;
	ENetPeer *peer;
	ENetHost *host;

	bool refuse_connections;

	ConnectionStatus connection_status;

	Map<int,ENetPeer*> peer_map;

	struct Packet {

		ENetPacket *packet;
		int from;
	};

	mutable List<Packet> incoming_packets;

	mutable Packet current_packet;

	uint32_t _gen_unique_id() const;
	void _pop_current_packet() const;

protected:
	static void _bind_methods();
public:

	virtual void set_transfer_mode(TransferMode p_mode);
	virtual void set_target_peer(int p_peer);


	virtual int get_packet_peer() const;


	Error create_server(int p_port, int p_max_peers=32, int p_in_bandwidth=0, int p_out_bandwidth=0);
	Error create_client(const IP_Address& p_ip, int p_port, int p_in_bandwidth=0, int p_out_bandwidth=0);

	void close_connection();

	virtual void poll();

	virtual bool is_server() const;

	virtual int get_available_packet_count() const;
	virtual Error get_packet(const uint8_t **r_buffer,int &r_buffer_size) const; ///< buffer is GONE after next get_packet
	virtual Error put_packet(const uint8_t *p_buffer,int p_buffer_size);

	virtual int get_max_packet_size() const;

	virtual ConnectionStatus get_connection_status() const;

	virtual void set_refuse_new_connections(bool p_enable);
	virtual bool is_refusing_new_connections() const;

	virtual int get_unique_id() const;

	NetworkedMultiplayerENet();
	~NetworkedMultiplayerENet();
};


#endif // NETWORKED_MULTIPLAYER_ENET_H
