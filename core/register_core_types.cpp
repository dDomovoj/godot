/*************************************************************************/
/*  register_core_types.cpp                                              */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                    http://www.godotengine.org                         */
/*************************************************************************/
/* Copyright (c) 2007-2017 Juan Linietsky, Ariel Manzur.                 */
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
#include "register_core_types.h"

#include "io/tcp_server.h"
#include "io/packet_peer_udp.h"
#include "io/config_file.h"
#include "os/main_loop.h"
#include "io/packet_peer.h"
#include "math/a_star.h"
#include "math/triangle_mesh.h"
#include "globals.h"
#include "object_type_db.h"
#include "geometry.h"
#include "bind/core_bind.h"
#include "core_string_names.h"
#include "path_remap.h"
#include "translation.h"
#include "compressed_translation.h"
#include "io/translation_loader_po.h"
#include "io/resource_format_xml.h"
#include "io/resource_format_binary.h"
#include "io/stream_peer_ssl.h"
#include "os/input.h"
#include "core/io/xml_parser.h"
#include "io/http_client.h"
#include "io/pck_packer.h"
#include "packed_data_container.h"
#include "func_ref.h"
#include "input_map.h"
#include "undo_redo.h"

#ifdef XML_ENABLED
static ResourceFormatSaverXML *resource_saver_xml=NULL;
static ResourceFormatLoaderXML *resource_loader_xml=NULL;
#endif
static ResourceFormatSaverBinary *resource_saver_binary=NULL;
static ResourceFormatLoaderBinary *resource_loader_binary=NULL;


static _ResourceLoader *_resource_loader=NULL;
static _ResourceSaver *_resource_saver=NULL;
static _OS *_os=NULL;
static _ClassDB *_classdb=NULL;
static _Marshalls *_marshalls = NULL;
static TranslationLoaderPO *resource_format_po=NULL;

static IP* ip = NULL;

static _Geometry *_geometry=NULL;

extern Mutex *_global_mutex;



extern void register_variant_methods();
extern void unregister_variant_methods();


void register_core_types() {


	_global_mutex=Mutex::create();


	StringName::setup();


	register_variant_methods();


	CoreStringNames::create();

	resource_format_po = memnew( TranslationLoaderPO );
	ResourceLoader::add_resource_format_loader( resource_format_po );


	resource_saver_binary = memnew( ResourceFormatSaverBinary );
	ResourceSaver::add_resource_format_saver(resource_saver_binary);
	resource_loader_binary = memnew( ResourceFormatLoaderBinary );
	ResourceLoader::add_resource_format_loader(resource_loader_binary);

#ifdef XML_ENABLED
	resource_saver_xml = memnew( ResourceFormatSaverXML );
	ResourceSaver::add_resource_format_saver(resource_saver_xml);
	resource_loader_xml = memnew( ResourceFormatLoaderXML );
	ResourceLoader::add_resource_format_loader(resource_loader_xml);
#endif

	ClassDB::register_class<Object>();


	ClassDB::register_class<Reference>();
	ClassDB::register_class<WeakRef>();
	ClassDB::register_class<ResourceImportMetadata>();
	ClassDB::register_class<Resource>();
	ClassDB::register_class<FuncRef>();
	ClassDB::register_virtual_class<StreamPeer>();
	ClassDB::register_class<StreamPeerBuffer>();
	ClassDB::register_custom_instance_class<StreamPeerTCP>();
	ClassDB::register_custom_instance_class<TCP_Server>();
	ClassDB::register_custom_instance_class<PacketPeerUDP>();
	ClassDB::register_custom_instance_class<StreamPeerSSL>();
	ClassDB::register_virtual_class<IP>();
	ClassDB::register_virtual_class<PacketPeer>();
	ClassDB::register_class<PacketPeerStream>();
	ClassDB::register_class<MainLoop>();
//	ClassDB::register_type<OptimizedSaver>();
	ClassDB::register_class<Translation>();
	ClassDB::register_class<PHashTranslation>();
	ClassDB::register_class<UndoRedo>();
	ClassDB::register_class<HTTPClient>();
	ClassDB::register_class<TriangleMesh>();

	ClassDB::register_virtual_class<ResourceInteractiveLoader>();

	ClassDB::register_class<_File>();
	ClassDB::register_class<_Directory>();
	ClassDB::register_class<_Thread>();
	ClassDB::register_class<_Mutex>();
	ClassDB::register_class<_Semaphore>();

	ClassDB::register_class<XMLParser>();

	ClassDB::register_class<ConfigFile>();

	ClassDB::register_class<PCKPacker>();

	ClassDB::register_class<PackedDataContainer>();
	ClassDB::register_virtual_class<PackedDataContainerRef>();
	ClassDB::register_class<AStar>();

	ip = IP::create();


	_geometry = memnew(_Geometry);


	_resource_loader=memnew(_ResourceLoader);
	_resource_saver=memnew(_ResourceSaver);
	_os=memnew(_OS);
	_classdb=memnew(_ClassDB);
	_marshalls = memnew(_Marshalls);



}

void register_core_singletons() {

	Globals::get_singleton()->add_singleton( Globals::Singleton("Globals",Globals::get_singleton()) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("IP",IP::get_singleton()) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("Geometry",_Geometry::get_singleton()) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("ResourceLoader",_ResourceLoader::get_singleton()) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("ResourceSaver",_ResourceSaver::get_singleton()) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("PathRemap",PathRemap::get_singleton() ) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("OS",_OS::get_singleton() ) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("ClassDB",_classdb ) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("Marshalls",_Marshalls::get_singleton() ) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("TranslationServer",TranslationServer::get_singleton() ) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("TS",TranslationServer::get_singleton() ) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("Input",Input::get_singleton() ) );
	Globals::get_singleton()->add_singleton( Globals::Singleton("InputMap",InputMap::get_singleton() )  );


}

void unregister_core_types() {



	memdelete( _resource_loader );
	memdelete( _resource_saver );
	memdelete( _os);
	memdelete( _classdb );
	memdelete( _marshalls );

	memdelete( _geometry );
#ifdef XML_ENABLED
	if (resource_saver_xml)
		memdelete(resource_saver_xml);
	if (resource_loader_xml)
		memdelete(resource_loader_xml);
#endif

	if (resource_saver_binary)
		memdelete(resource_saver_binary);
	if (resource_loader_binary)
		memdelete(resource_loader_binary);


	memdelete( resource_format_po );

	if (ip)
		memdelete(ip);


	ObjectDB::cleanup();

	unregister_variant_methods();

	ClassDB::cleanup();
	ResourceCache::clear();
	CoreStringNames::free();
	StringName::cleanup();

	if (_global_mutex) {
		memdelete(_global_mutex);
		_global_mutex=NULL; //still needed at a few places
	};
}
