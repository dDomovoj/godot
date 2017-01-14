/*************************************************************************/
/*  resource.cpp                                                         */
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
#include "resource.h"

#include "core_string_names.h"
#include "os/file_access.h"
#include "io/resource_loader.h"
#include "script_language.h"

#include <stdio.h>

void ResourceImportMetadata::set_editor(const String& p_editor) {

	editor=p_editor;
}

String ResourceImportMetadata::get_editor() const{

	return editor;
}

void ResourceImportMetadata::add_source(const String& p_path,const String& p_md5) {

	Source s;
	s.md5=p_md5;
	s.path=p_path;
	sources.push_back(s);
}

String ResourceImportMetadata::get_source_path(int p_idx) const{
	ERR_FAIL_INDEX_V(p_idx,sources.size(),String());
	return sources[p_idx].path;
}
String ResourceImportMetadata::get_source_md5(int p_idx) const{
	ERR_FAIL_INDEX_V(p_idx,sources.size(),String());
	return sources[p_idx].md5;
}

void ResourceImportMetadata::set_source_md5(int p_idx,const String& p_md5) {

	ERR_FAIL_INDEX(p_idx,sources.size());
	sources[p_idx].md5=p_md5;

}

void ResourceImportMetadata::remove_source(int p_idx){

	ERR_FAIL_INDEX(p_idx,sources.size());
	sources.remove(p_idx);

}

int ResourceImportMetadata::get_source_count() const {

	return sources.size();
}
void ResourceImportMetadata::set_option(const String& p_key, const Variant& p_value) {

	if (p_value.get_type()==Variant::NIL) {
		options.erase(p_key);
		return;
	}

	ERR_FAIL_COND(p_value.get_type() == Variant::OBJECT);
	ERR_FAIL_COND(p_value.get_type() == Variant::_RID);

	options[p_key]=p_value;

}

bool ResourceImportMetadata::has_option(const String& p_key) const {

	return options.has(p_key);
}

Variant ResourceImportMetadata::get_option(const String& p_key) const {

	ERR_FAIL_COND_V(!options.has(p_key),Variant());

	return options[p_key];
}

void ResourceImportMetadata::get_options(List<String> *r_options) const {

	for(Map<String,Variant>::Element *E=options.front();E;E=E->next()) {

		r_options->push_back(E->key());
	}

}

PoolStringArray ResourceImportMetadata::_get_options() const {

	PoolStringArray option_names;
	option_names.resize(options.size());
	int i=0;
	for(Map<String,Variant>::Element *E=options.front();E;E=E->next()) {

		option_names.set(i++,E->key());
	}

	return option_names;
}

void ResourceImportMetadata::_bind_methods() {

	ClassDB::bind_method(_MD("set_editor","name"),&ResourceImportMetadata::set_editor);
	ClassDB::bind_method(_MD("get_editor"),&ResourceImportMetadata::get_editor);
	ClassDB::bind_method(_MD("add_source","path","md5"),&ResourceImportMetadata::add_source, "");
	ClassDB::bind_method(_MD("get_source_path","idx"),&ResourceImportMetadata::get_source_path);
	ClassDB::bind_method(_MD("get_source_md5","idx"),&ResourceImportMetadata::get_source_md5);
	ClassDB::bind_method(_MD("set_source_md5","idx", "md5"),&ResourceImportMetadata::set_source_md5);
	ClassDB::bind_method(_MD("remove_source","idx"),&ResourceImportMetadata::remove_source);
	ClassDB::bind_method(_MD("get_source_count"),&ResourceImportMetadata::get_source_count);
	ClassDB::bind_method(_MD("set_option","key","value"),&ResourceImportMetadata::set_option);
	ClassDB::bind_method(_MD("get_option","key"),&ResourceImportMetadata::get_option);
	ClassDB::bind_method(_MD("get_options"),&ResourceImportMetadata::_get_options);
}

ResourceImportMetadata::ResourceImportMetadata() {


}


void Resource::emit_changed() {

	emit_signal(CoreStringNames::get_singleton()->changed);
}


void Resource::_resource_path_changed() {


}

void Resource::set_path(const String& p_path, bool p_take_over) {

	if (path_cache==p_path)
		return;

	if (path_cache!="") {

		ResourceCache::lock->write_lock();
		ResourceCache::resources.erase(path_cache);
		ResourceCache::lock->write_unlock();
	}

	path_cache="";

	ResourceCache::lock->read_lock();
	bool has_path = ResourceCache::resources.has( p_path );
	ResourceCache::lock->read_unlock();

	if (has_path) {
		if (p_take_over) {

			ResourceCache::lock->write_lock();
			ResourceCache::resources.get(p_path)->set_name("");
			ResourceCache::lock->write_unlock();
		} else {
			ERR_EXPLAIN("Another resource is loaded from path: "+p_path);

			ResourceCache::lock->read_lock();
			bool exists = ResourceCache::resources.has( p_path );
			ResourceCache::lock->read_unlock();

			ERR_FAIL_COND( exists );
		}

	}
	path_cache=p_path;

	if (path_cache!="") {

		ResourceCache::lock->write_lock();
		ResourceCache::resources[path_cache]=this;
		ResourceCache::lock->write_unlock();
	}

	_change_notify("resource_path");
	_resource_path_changed();

}

String Resource::get_path() const {

	return path_cache;
}

void Resource::set_subindex(int p_sub_index) {

	subindex=p_sub_index;
}

int Resource::get_subindex() const{

	return subindex;
}


void Resource::set_name(const String& p_name) {

	name=p_name;
	_change_notify("resource_name");

}
String Resource::get_name() const {

	return name;
}

bool Resource::editor_can_reload_from_file() {

	return true; //by default yes
}

void Resource::reload_from_file() {


	String path=get_path();
	if (!path.is_resource_file())
		return;

	Ref<Resource> s = ResourceLoader::load(path,get_class(),true);

	if (!s.is_valid())
		return;

	List<PropertyInfo> pi;
	s->get_property_list(&pi);

	for (List<PropertyInfo>::Element *E=pi.front();E;E=E->next()) {

		if (!(E->get().usage&PROPERTY_USAGE_STORAGE))
			continue;
		if (E->get().name=="resource_path")
			continue; //do not change path

		set(E->get().name,s->get(E->get().name));

	}
}


Ref<Resource> Resource::duplicate_for_local_scene(Node *p_for_scene, Map<Ref<Resource>, Ref<Resource> > &remap_cache) {


	List<PropertyInfo> plist;
	get_property_list(&plist);


	Resource *r = (Resource*)ClassDB::instance(get_class());
	ERR_FAIL_COND_V(!r,Ref<Resource>());

	r->local_scene=p_for_scene;

	for(List<PropertyInfo>::Element *E=plist.front();E;E=E->next()) {

		if (!(E->get().usage&PROPERTY_USAGE_STORAGE))
			continue;
		Variant p = get(E->get().name);
		if (p.get_type()==Variant::OBJECT) {

			RES sr = p;
			if (sr.is_valid()) {

				if (sr->is_local_to_scene()) {
					if (remap_cache.has(sr)) {
						p=remap_cache[sr];
					} else {


						RES dupe = sr->duplicate_for_local_scene(p_for_scene,remap_cache);
						p=dupe;
						remap_cache[sr]=dupe;
					}
				}
			}
		}

		r->set(E->get().name,p);
	}

	return Ref<Resource>(r);
}

Ref<Resource> Resource::duplicate(bool p_subresources) {


	List<PropertyInfo> plist;
	get_property_list(&plist);


	Resource *r = (Resource*)ClassDB::instance(get_class());
	ERR_FAIL_COND_V(!r,Ref<Resource>());

	for(List<PropertyInfo>::Element *E=plist.front();E;E=E->next()) {

		if (!(E->get().usage&PROPERTY_USAGE_STORAGE))
			continue;
		Variant p = get(E->get().name);
		if (p.get_type()==Variant::OBJECT && p_subresources) {

			RES sr = p;
			if (sr.is_valid())
				p=sr->duplicate(true);
		}

		r->set(E->get().name,p);
	}

	return Ref<Resource>(r);
}


void Resource::_set_path(const String& p_path) {

	set_path(p_path,false);
}

void Resource::_take_over_path(const String& p_path) {

	set_path(p_path,true);
}


RID Resource::get_rid() const {

	return RID();
}


void Resource::register_owner(Object *p_owner) {

	owners.insert(p_owner->get_instance_ID());
}

void Resource::unregister_owner(Object *p_owner) {

	owners.erase(p_owner->get_instance_ID());

}


void Resource::notify_change_to_owners() {

	for(Set<ObjectID>::Element *E=owners.front();E;E=E->next()) {

		Object *obj = ObjectDB::get_instance(E->get());
		ERR_EXPLAIN("Object was deleted, while still owning a resource");
		ERR_CONTINUE(!obj); //wtf
		//TODO store string
		obj->call("resource_changed",RES(this));
	}
}

void Resource::set_import_metadata(const Ref<ResourceImportMetadata>& p_metadata) {
#ifdef TOOLS_ENABLED
	import_metadata=p_metadata;
#endif
}

Ref<ResourceImportMetadata> Resource::get_import_metadata() const {

#ifdef TOOLS_ENABLED
	return import_metadata;
#else
	return Ref<ResourceImportMetadata>();
#endif

}

#ifdef TOOLS_ENABLED

uint32_t Resource::hash_edited_version() const {

	uint32_t hash = hash_djb2_one_32(get_edited_version());

	List<PropertyInfo> plist;
	get_property_list(&plist);

	for (List<PropertyInfo>::Element *E=plist.front();E;E=E->next()) {

		if (E->get().type==Variant::OBJECT && E->get().hint==PROPERTY_HINT_RESOURCE_TYPE) {
			RES res = get(E->get().name);
			if (res.is_valid()) {
				hash = hash_djb2_one_32(res->hash_edited_version(),hash);
			}
		}
	}

	return hash;

}

#endif

void Resource::set_local_to_scene(bool p_enable) {

	local_to_scene=p_enable;
}

bool Resource::is_local_to_scene() const {

	return local_to_scene;
}

Node* Resource::get_local_scene() const {

	if (local_scene)
		return local_scene;

	if (_get_local_scene_func) {
		return _get_local_scene_func();
	}

	return NULL;
}

void Resource::setup_local_to_scene() {

	if (get_script_instance())
		get_script_instance()->call("_setup_local_to_scene");
}

Node* (*Resource::_get_local_scene_func)()=NULL;


void Resource::_bind_methods() {

	ClassDB::bind_method(_MD("set_path","path"),&Resource::_set_path);
	ClassDB::bind_method(_MD("take_over_path","path"),&Resource::_take_over_path);
	ClassDB::bind_method(_MD("get_path"),&Resource::get_path);
	ClassDB::bind_method(_MD("set_name","name"),&Resource::set_name);
	ClassDB::bind_method(_MD("get_name"),&Resource::get_name);
	ClassDB::bind_method(_MD("get_rid"),&Resource::get_rid);
	ClassDB::bind_method(_MD("set_import_metadata","metadata"),&Resource::set_import_metadata);
	ClassDB::bind_method(_MD("get_import_metadata"),&Resource::get_import_metadata);
	ClassDB::bind_method(_MD("set_local_to_scene","enable"),&Resource::set_local_to_scene);
	ClassDB::bind_method(_MD("is_local_to_scene"),&Resource::is_local_to_scene);
	ClassDB::bind_method(_MD("get_local_scene:Node"),&Resource::get_local_scene);
	ClassDB::bind_method(_MD("setup_local_to_scene"),&Resource::setup_local_to_scene);

	ClassDB::bind_method(_MD("duplicate","subresources"),&Resource::duplicate,DEFVAL(false));
	ADD_SIGNAL( MethodInfo("changed") );
	ADD_GROUP("Resource","resource_");
	ADD_PROPERTYNZ( PropertyInfo(Variant::BOOL,"resource_local_to_scene" ), _SCS("set_local_to_scene"),_SCS("is_local_to_scene"));
	ADD_PROPERTY( PropertyInfo(Variant::STRING,"resource_path",PROPERTY_HINT_NONE,"",PROPERTY_USAGE_EDITOR ), _SCS("set_path"),_SCS("get_path"));
	ADD_PROPERTYNZ( PropertyInfo(Variant::STRING,"resource_name"), _SCS("set_name"),_SCS("get_name"));

	BIND_VMETHOD( MethodInfo("_setup_local_to_scene") );

}

Resource::Resource() {

#ifdef TOOLS_ENABLED
	last_modified_time=0;
#endif

	subindex=0;
	local_scene=NULL;
}




Resource::~Resource() {

	if (path_cache!="") {
		ResourceCache::lock->write_lock();
		ResourceCache::resources.erase(path_cache);
		ResourceCache::lock->write_unlock();
	}
	if (owners.size()) {
		WARN_PRINT("Resource is still owned");
	}
}

HashMap<String,Resource*> ResourceCache::resources;

RWLock *ResourceCache::lock=NULL;

void ResourceCache::setup() {

	lock = RWLock::create();
}

void ResourceCache::clear() {
	if (resources.size())
		ERR_PRINT("Resources Still in use at Exit!");

	resources.clear();
	memdelete(lock);
}


void ResourceCache::reload_externals() {

	/*
	const String *K=NULL;
	while ((K=resources.next(K))) {
		resources[*K]->reload_external_data();
	}
	*/
}

bool ResourceCache::has(const String& p_path) {

	lock->read_lock();
	bool b = resources.has(p_path);
	lock->read_unlock();


	return b;
}
Resource *ResourceCache::get(const String& p_path) {

	lock->read_lock();

	Resource **res = resources.getptr(p_path);

	lock->read_unlock();

	if (!res) {
		return NULL;
	}

	return *res;
}


void ResourceCache::get_cached_resources(List<Ref<Resource> > *p_resources) {


	lock->read_lock();
	const String* K=NULL;
	while((K=resources.next(K))) {

		Resource *r = resources[*K];
		p_resources->push_back( Ref<Resource>( r ));

	}
	lock->read_unlock();

}

int ResourceCache::get_cached_resource_count() {

	lock->read_lock();
	int rc = resources.size();
	lock->read_unlock();

	return rc;
}

void ResourceCache::dump(const char* p_file,bool p_short) {
#ifdef DEBUG_ENABLED
	lock->read_lock();

	Map<String,int> type_count;


	FileAccess *f=NULL;
	if (p_file) {
		f = FileAccess::open(p_file,FileAccess::WRITE);
		ERR_FAIL_COND(!f);

	}

	const String* K=NULL;
	while((K=resources.next(K))) {

		Resource *r = resources[*K];

		if (!type_count.has(r->get_class())) {
			type_count[r->get_class()]=0;
		}


		type_count[r->get_class()]++;

		if (!p_short) {
			if (f)
				f->store_line(r->get_class()+": "+r->get_path());
		}
	}

	for(Map<String,int>::Element *E=type_count.front();E;E=E->next()) {

		if (f)
			f->store_line(E->key()+" count: "+itos(E->get()));
	}
	if (f) {
		f->close();
		memdelete(f);
	}

	lock->read_unlock();

#endif
}
