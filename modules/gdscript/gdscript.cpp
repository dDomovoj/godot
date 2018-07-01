/*************************************************************************/
/*  gdscript.cpp                                                         */
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

#include "gdscript.h"

#include "engine.h"
#include "gdscript_compiler.h"
#include "global_constants.h"
#include "io/file_access_encrypted.h"
#include "os/file_access.h"
#include "os/os.h"
#include "project_settings.h"

///////////////////////////

GDScriptNativeClass::GDScriptNativeClass(const StringName &p_name) {

	name = p_name;
}

bool GDScriptNativeClass::_get(const StringName &p_name, Variant &r_ret) const {

	bool ok;
	int v = ClassDB::get_integer_constant(name, p_name, &ok);

	if (ok) {
		r_ret = v;
		return true;
	} else {
		return false;
	}
}

void GDScriptNativeClass::_bind_methods() {

	ClassDB::bind_method(D_METHOD("new"), &GDScriptNativeClass::_new);
}

Variant GDScriptNativeClass::_new() {

	Object *o = instance();
	if (!o) {
		ERR_EXPLAIN("Class type: '" + String(name) + "' is not instantiable.");
		ERR_FAIL_COND_V(!o, Variant());
	}

	Reference *ref = Object::cast_to<Reference>(o);
	if (ref) {
		return REF(ref);
	} else {
		return o;
	}
}

Object *GDScriptNativeClass::instance() {

	return ClassDB::instance(name);
}

GDScriptInstance *GDScript::_create_instance(const Variant **p_args, int p_argcount, Object *p_owner, bool p_isref, Variant::CallError &r_error) {

	/* STEP 1, CREATE */

	GDScriptInstance *instance = memnew(GDScriptInstance);
	instance->base_ref = p_isref;
	instance->members.resize(member_indices.size());
	instance->script = Ref<GDScript>(this);
	instance->owner = p_owner;
#ifdef DEBUG_ENABLED
	//needed for hot reloading
	for (Map<StringName, MemberInfo>::Element *E = member_indices.front(); E; E = E->next()) {
		instance->member_indices_cache[E->key()] = E->get().index;
	}
#endif
	instance->owner->set_script_instance(instance);

	/* STEP 2, INITIALIZE AND CONSRTUCT */

#ifndef NO_THREADS
	GDScriptLanguage::singleton->lock->lock();
#endif

	instances.insert(instance->owner);

#ifndef NO_THREADS
	GDScriptLanguage::singleton->lock->unlock();
#endif

	initializer->call(instance, p_args, p_argcount, r_error);

	if (r_error.error != Variant::CallError::CALL_OK) {
		instance->script = Ref<GDScript>();
		instance->owner->set_script_instance(NULL);
#ifndef NO_THREADS
		GDScriptLanguage::singleton->lock->lock();
#endif
		instances.erase(p_owner);
#ifndef NO_THREADS
		GDScriptLanguage::singleton->lock->unlock();
#endif

		ERR_FAIL_COND_V(r_error.error != Variant::CallError::CALL_OK, NULL); //error constructing
	}

	//@TODO make thread safe
	return instance;
}

Variant GDScript::_new(const Variant **p_args, int p_argcount, Variant::CallError &r_error) {

	/* STEP 1, CREATE */

	if (!valid) {
		r_error.error = Variant::CallError::CALL_ERROR_INVALID_METHOD;
		return Variant();
	}

	r_error.error = Variant::CallError::CALL_OK;
	REF ref;
	Object *owner = NULL;

	GDScript *_baseptr = this;
	while (_baseptr->_base) {
		_baseptr = _baseptr->_base;
	}

	ERR_FAIL_COND_V(_baseptr->native.is_null(), Variant());

	if (_baseptr->native.ptr()) {
		owner = _baseptr->native->instance();
	} else {
		owner = memnew(Reference); //by default, no base means use reference
	}

	Reference *r = Object::cast_to<Reference>(owner);
	if (r) {
		ref = REF(r);
	}

	GDScriptInstance *instance = _create_instance(p_args, p_argcount, owner, r != NULL, r_error);
	if (!instance) {
		if (ref.is_null()) {
			memdelete(owner); //no owner, sorry
		}
		return Variant();
	}

	if (ref.is_valid()) {
		return ref;
	} else {
		return owner;
	}
}

bool GDScript::can_instance() const {

	return valid || (!tool && !ScriptServer::is_scripting_enabled());
}

Ref<Script> GDScript::get_base_script() const {

	if (_base) {
		return Ref<GDScript>(_base);
	} else {
		return Ref<Script>();
	}
}

StringName GDScript::get_instance_base_type() const {

	if (native.is_valid())
		return native->get_name();
	if (base.is_valid())
		return base->get_instance_base_type();
	return StringName();
}

struct _GDScriptMemberSort {

	int index;
	StringName name;
	_FORCE_INLINE_ bool operator<(const _GDScriptMemberSort &p_member) const { return index < p_member.index; }
};

#ifdef TOOLS_ENABLED

void GDScript::_placeholder_erased(PlaceHolderScriptInstance *p_placeholder) {

	placeholders.erase(p_placeholder);
}
#endif

void GDScript::get_script_method_list(List<MethodInfo> *p_list) const {

	for (const Map<StringName, GDScriptFunction *>::Element *E = member_functions.front(); E; E = E->next()) {
		GDScriptFunction *func = E->get();
		MethodInfo mi;
		mi.name = E->key();
		for (int i = 0; i < func->get_argument_count(); i++) {
			mi.arguments.push_back(func->get_argument_type(i));
		}

		mi.return_val = func->get_return_type();
		p_list->push_back(mi);
	}
}

void GDScript::get_script_property_list(List<PropertyInfo> *p_list) const {

	const GDScript *sptr = this;
	List<PropertyInfo> props;

	while (sptr) {

		Vector<_GDScriptMemberSort> msort;
		for (Map<StringName, PropertyInfo>::Element *E = sptr->member_info.front(); E; E = E->next()) {

			_GDScriptMemberSort ms;
			ERR_CONTINUE(!sptr->member_indices.has(E->key()));
			ms.index = sptr->member_indices[E->key()].index;
			ms.name = E->key();
			msort.push_back(ms);
		}

		msort.sort();
		msort.invert();
		for (int i = 0; i < msort.size(); i++) {

			props.push_front(sptr->member_info[msort[i].name]);
		}

		sptr = sptr->_base;
	}

	for (List<PropertyInfo>::Element *E = props.front(); E; E = E->next()) {
		p_list->push_back(E->get());
	}
}

bool GDScript::has_method(const StringName &p_method) const {

	return member_functions.has(p_method);
}

MethodInfo GDScript::get_method_info(const StringName &p_method) const {

	const Map<StringName, GDScriptFunction *>::Element *E = member_functions.find(p_method);
	if (!E)
		return MethodInfo();

	GDScriptFunction *func = E->get();
	MethodInfo mi;
	mi.name = E->key();
	for (int i = 0; i < func->get_argument_count(); i++) {
		mi.arguments.push_back(func->get_argument_type(i));
	}

	mi.return_val = func->get_return_type();
	return mi;
}

bool GDScript::get_property_default_value(const StringName &p_property, Variant &r_value) const {

#ifdef TOOLS_ENABLED

	/*
	for (const Map<StringName,Variant>::Element *I=member_default_values.front();I;I=I->next()) {
		print_line("\t"+String(String(I->key())+":"+String(I->get())));
	}
	*/
	const Map<StringName, Variant>::Element *E = member_default_values_cache.find(p_property);
	if (E) {
		r_value = E->get();
		return true;
	}

	if (base_cache.is_valid()) {
		return base_cache->get_property_default_value(p_property, r_value);
	}
#endif
	return false;
}

ScriptInstance *GDScript::instance_create(Object *p_this) {

	if (!tool && !ScriptServer::is_scripting_enabled()) {

#ifdef TOOLS_ENABLED

		//instance a fake script for editing the values
		//plist.invert();

		/*print_line("CREATING PLACEHOLDER");
		for(List<PropertyInfo>::Element *E=plist.front();E;E=E->next()) {
			print_line(E->get().name);
		}*/
		PlaceHolderScriptInstance *si = memnew(PlaceHolderScriptInstance(GDScriptLanguage::get_singleton(), Ref<Script>(this), p_this));
		placeholders.insert(si);
		//_update_placeholder(si);
		_update_exports();
		return si;
#else
		return NULL;
#endif
	}

	GDScript *top = this;
	while (top->_base)
		top = top->_base;

	if (top->native.is_valid()) {
		if (!ClassDB::is_parent_class(p_this->get_class_name(), top->native->get_name())) {

			if (ScriptDebugger::get_singleton()) {
				GDScriptLanguage::get_singleton()->debug_break_parse(get_path(), 0, "Script inherits from native type '" + String(top->native->get_name()) + "', so it can't be instanced in object of type: '" + p_this->get_class() + "'");
			}
			ERR_EXPLAIN("Script inherits from native type '" + String(top->native->get_name()) + "', so it can't be instanced in object of type: '" + p_this->get_class() + "'");
			ERR_FAIL_V(NULL);
		}
	}

	Variant::CallError unchecked_error;
	return _create_instance(NULL, 0, p_this, Object::cast_to<Reference>(p_this), unchecked_error);
}
bool GDScript::instance_has(const Object *p_this) const {

#ifndef NO_THREADS
	GDScriptLanguage::singleton->lock->lock();
#endif
	bool hasit = instances.has((Object *)p_this);

#ifndef NO_THREADS
	GDScriptLanguage::singleton->lock->unlock();
#endif

	return hasit;
}

bool GDScript::has_source_code() const {

	return source != "";
}
String GDScript::get_source_code() const {

	return source;
}
void GDScript::set_source_code(const String &p_code) {

	if (source == p_code)
		return;
	source = p_code;
#ifdef TOOLS_ENABLED
	source_changed_cache = true;
#endif
}

#ifdef TOOLS_ENABLED
void GDScript::_update_exports_values(Map<StringName, Variant> &values, List<PropertyInfo> &propnames) {

	if (base_cache.is_valid()) {
		base_cache->_update_exports_values(values, propnames);
	}

	for (Map<StringName, Variant>::Element *E = member_default_values_cache.front(); E; E = E->next()) {
		values[E->key()] = E->get();
	}

	for (List<PropertyInfo>::Element *E = members_cache.front(); E; E = E->next()) {
		propnames.push_back(E->get());
	}
}
#endif

bool GDScript::_update_exports() {

#ifdef TOOLS_ENABLED

	bool changed = false;

	if (source_changed_cache) {
		//print_line("updating source for "+get_path());
		source_changed_cache = false;
		changed = true;

		String basedir = path;

		if (basedir == "")
			basedir = get_path();

		if (basedir != "")
			basedir = basedir.get_base_dir();

		GDScriptParser parser;
		Error err = parser.parse(source, basedir, true, path);

		if (err == OK) {

			const GDScriptParser::Node *root = parser.get_parse_tree();
			ERR_FAIL_COND_V(root->type != GDScriptParser::Node::TYPE_CLASS, false);

			const GDScriptParser::ClassNode *c = static_cast<const GDScriptParser::ClassNode *>(root);

			if (base_cache.is_valid()) {
				base_cache->inheriters_cache.erase(get_instance_id());
				base_cache = Ref<GDScript>();
			}

			if (c->extends_used && String(c->extends_file) != "" && String(c->extends_file) != get_path()) {

				String path = c->extends_file;
				if (path.is_rel_path()) {

					String base = get_path();
					if (base == "" || base.is_rel_path()) {

						ERR_PRINT(("Could not resolve relative path for parent class: " + path).utf8().get_data());
					} else {
						path = base.get_base_dir().plus_file(path);
					}
				}

				if (path != get_path()) {

					Ref<GDScript> bf = ResourceLoader::load(path);

					if (bf.is_valid()) {

						//print_line("parent is: "+bf->get_path());
						base_cache = bf;
						bf->inheriters_cache.insert(get_instance_id());

						//bf->_update_exports(p_instances,true,false);
					}
				} else {
					ERR_PRINT(("Path extending itself in  " + path).utf8().get_data());
				}
			}

			members_cache.clear();
			member_default_values_cache.clear();

			for (int i = 0; i < c->variables.size(); i++) {
				if (c->variables[i]._export.type == Variant::NIL)
					continue;

				members_cache.push_back(c->variables[i]._export);
				//print_line("found "+c->variables[i]._export.name);
				member_default_values_cache[c->variables[i].identifier] = c->variables[i].default_value;
			}

			_signals.clear();

			for (int i = 0; i < c->_signals.size(); i++) {
				_signals[c->_signals[i].name] = c->_signals[i].arguments;
			}
		}
	} else {
		//print_line("unchanged is "+get_path());
	}

	if (base_cache.is_valid()) {
		if (base_cache->_update_exports()) {
			changed = true;
		}
	}

	if (/*changed &&*/ placeholders.size()) { //hm :(

		//print_line("updating placeholders for "+get_path());

		//update placeholders if any
		Map<StringName, Variant> values;
		List<PropertyInfo> propnames;
		_update_exports_values(values, propnames);

		for (Set<PlaceHolderScriptInstance *>::Element *E = placeholders.front(); E; E = E->next()) {

			E->get()->update(propnames, values);
		}
	}

	return changed;

#else
	return false;
#endif
}

void GDScript::update_exports() {

#ifdef TOOLS_ENABLED

	_update_exports();

	Set<ObjectID> copy = inheriters_cache; //might get modified

	//print_line("update exports for "+get_path()+" ic: "+itos(copy.size()));
	for (Set<ObjectID>::Element *E = copy.front(); E; E = E->next()) {
		Object *id = ObjectDB::get_instance(E->get());
		GDScript *s = Object::cast_to<GDScript>(id);
		if (!s)
			continue;
		s->update_exports();
	}

#endif
}

void GDScript::_set_subclass_path(Ref<GDScript> &p_sc, const String &p_path) {

	p_sc->path = p_path;
	for (Map<StringName, Ref<GDScript> >::Element *E = p_sc->subclasses.front(); E; E = E->next()) {

		_set_subclass_path(E->get(), p_path);
	}
}

Error GDScript::reload(bool p_keep_state) {

#ifndef NO_THREADS
	GDScriptLanguage::singleton->lock->lock();
#endif
	bool has_instances = instances.size();

#ifndef NO_THREADS
	GDScriptLanguage::singleton->lock->unlock();
#endif

	ERR_FAIL_COND_V(!p_keep_state && has_instances, ERR_ALREADY_IN_USE);

	String basedir = path;

	if (basedir == "")
		basedir = get_path();

	if (basedir != "")
		basedir = basedir.get_base_dir();

	if (source.find("%BASE%") != -1) {
		//loading a template, don't parse
		return OK;
	}

	valid = false;
	GDScriptParser parser;
	Error err = parser.parse(source, basedir, false, path);
	if (err) {
		if (ScriptDebugger::get_singleton()) {
			GDScriptLanguage::get_singleton()->debug_break_parse(get_path(), parser.get_error_line(), "Parser Error: " + parser.get_error());
		}
		_err_print_error("GDScript::reload", path.empty() ? "built-in" : (const char *)path.utf8().get_data(), parser.get_error_line(), ("Parse Error: " + parser.get_error()).utf8().get_data(), ERR_HANDLER_SCRIPT);
		ERR_FAIL_V(ERR_PARSE_ERROR);
	}

	bool can_run = ScriptServer::is_scripting_enabled() || parser.is_tool_script();

	GDScriptCompiler compiler;
	err = compiler.compile(&parser, this, p_keep_state);

	if (err) {

		if (can_run) {
			if (ScriptDebugger::get_singleton()) {
				GDScriptLanguage::get_singleton()->debug_break_parse(get_path(), compiler.get_error_line(), "Parser Error: " + compiler.get_error());
			}
			_err_print_error("GDScript::reload", path.empty() ? "built-in" : (const char *)path.utf8().get_data(), compiler.get_error_line(), ("Compile Error: " + compiler.get_error()).utf8().get_data(), ERR_HANDLER_SCRIPT);
			ERR_FAIL_V(ERR_COMPILATION_FAILED);
		} else {
			return err;
		}
	}
#if DEBUG_ENABLED
	for (const List<GDScriptWarning>::Element *E = parser.get_warnings().front(); E; E = E->next()) {
		String msg = "Script warning: " + E->get().get_name() + " (" + path + ") line " + itos(E->get().line) + ": ";
		msg += E->get().get_message();
		WARN_PRINTS(msg);
	}
#endif

	valid = true;

	for (Map<StringName, Ref<GDScript> >::Element *E = subclasses.front(); E; E = E->next()) {

		_set_subclass_path(E->get(), path);
	}

	return OK;
}

ScriptLanguage *GDScript::get_language() const {

	return GDScriptLanguage::get_singleton();
}

void GDScript::get_constants(Map<StringName, Variant> *p_constants) {

	if (p_constants) {
		for (Map<StringName, Variant>::Element *E = constants.front(); E; E = E->next()) {
			(*p_constants)[E->key()] = E->value();
		}
	}
}

void GDScript::get_members(Set<StringName> *p_members) {
	if (p_members) {
		for (Set<StringName>::Element *E = members.front(); E; E = E->next()) {
			p_members->insert(E->get());
		}
	}
}

Variant GDScript::call(const StringName &p_method, const Variant **p_args, int p_argcount, Variant::CallError &r_error) {

	GDScript *top = this;
	while (top) {

		Map<StringName, GDScriptFunction *>::Element *E = top->member_functions.find(p_method);
		if (E) {

			if (!E->get()->is_static()) {
				WARN_PRINT(String("Can't call non-static function: '" + String(p_method) + "' in script.").utf8().get_data());
			}

			return E->get()->call(NULL, p_args, p_argcount, r_error);
		}
		top = top->_base;
	}

	//none found, regular

	return Script::call(p_method, p_args, p_argcount, r_error);
}

bool GDScript::_get(const StringName &p_name, Variant &r_ret) const {

	{

		const GDScript *top = this;
		while (top) {

			{
				const Map<StringName, Variant>::Element *E = top->constants.find(p_name);
				if (E) {

					r_ret = E->get();
					return true;
				}
			}

			{
				const Map<StringName, Ref<GDScript> >::Element *E = subclasses.find(p_name);
				if (E) {

					r_ret = E->get();
					return true;
				}
			}
			top = top->_base;
		}

		if (p_name == GDScriptLanguage::get_singleton()->strings._script_source) {

			r_ret = get_source_code();
			return true;
		}
	}

	return false;
}
bool GDScript::_set(const StringName &p_name, const Variant &p_value) {

	if (p_name == GDScriptLanguage::get_singleton()->strings._script_source) {

		set_source_code(p_value);
		reload();
	} else
		return false;

	return true;
}

void GDScript::_get_property_list(List<PropertyInfo> *p_properties) const {

	p_properties->push_back(PropertyInfo(Variant::STRING, "script/source", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NOEDITOR | PROPERTY_USAGE_INTERNAL));
}

void GDScript::_bind_methods() {

	ClassDB::bind_vararg_method(METHOD_FLAGS_DEFAULT, "new", &GDScript::_new, MethodInfo(Variant::OBJECT, "new"));

	ClassDB::bind_method(D_METHOD("get_as_byte_code"), &GDScript::get_as_byte_code);
}

Vector<uint8_t> GDScript::get_as_byte_code() const {

	GDScriptTokenizerBuffer tokenizer;
	return tokenizer.parse_code_string(source);
};

Error GDScript::load_byte_code(const String &p_path) {

	Vector<uint8_t> bytecode;

	if (p_path.ends_with("gde")) {

		FileAccess *fa = FileAccess::open(p_path, FileAccess::READ);
		ERR_FAIL_COND_V(!fa, ERR_CANT_OPEN);
		FileAccessEncrypted *fae = memnew(FileAccessEncrypted);
		ERR_FAIL_COND_V(!fae, ERR_CANT_OPEN);
		Vector<uint8_t> key;
		key.resize(32);
		for (int i = 0; i < key.size(); i++) {
			key.write[i] = script_encryption_key[i];
		}
		Error err = fae->open_and_parse(fa, key, FileAccessEncrypted::MODE_READ);
		ERR_FAIL_COND_V(err, err);
		bytecode.resize(fae->get_len());
		fae->get_buffer(bytecode.ptrw(), bytecode.size());
		memdelete(fae);
	} else {

		bytecode = FileAccess::get_file_as_array(p_path);
	}
	ERR_FAIL_COND_V(bytecode.size() == 0, ERR_PARSE_ERROR);
	path = p_path;

	String basedir = path;

	if (basedir == "")
		basedir = get_path();

	if (basedir != "")
		basedir = basedir.get_base_dir();

	valid = false;
	GDScriptParser parser;
	Error err = parser.parse_bytecode(bytecode, basedir, get_path());
	if (err) {
		_err_print_error("GDScript::load_byte_code", path.empty() ? "built-in" : (const char *)path.utf8().get_data(), parser.get_error_line(), ("Parse Error: " + parser.get_error()).utf8().get_data(), ERR_HANDLER_SCRIPT);
		ERR_FAIL_V(ERR_PARSE_ERROR);
	}

	GDScriptCompiler compiler;
	err = compiler.compile(&parser, this);

	if (err) {
		_err_print_error("GDScript::load_byte_code", path.empty() ? "built-in" : (const char *)path.utf8().get_data(), compiler.get_error_line(), ("Compile Error: " + compiler.get_error()).utf8().get_data(), ERR_HANDLER_SCRIPT);
		ERR_FAIL_V(ERR_COMPILATION_FAILED);
	}

	valid = true;

	for (Map<StringName, Ref<GDScript> >::Element *E = subclasses.front(); E; E = E->next()) {

		_set_subclass_path(E->get(), path);
	}

	return OK;
}

Error GDScript::load_source_code(const String &p_path) {

	PoolVector<uint8_t> sourcef;
	Error err;
	FileAccess *f = FileAccess::open(p_path, FileAccess::READ, &err);
	if (err) {

		ERR_FAIL_COND_V(err, err);
	}

	int len = f->get_len();
	sourcef.resize(len + 1);
	PoolVector<uint8_t>::Write w = sourcef.write();
	int r = f->get_buffer(w.ptr(), len);
	f->close();
	memdelete(f);
	ERR_FAIL_COND_V(r != len, ERR_CANT_OPEN);
	w[len] = 0;

	String s;
	if (s.parse_utf8((const char *)w.ptr())) {

		ERR_EXPLAIN("Script '" + p_path + "' contains invalid unicode (utf-8), so it was not loaded. Please ensure that scripts are saved in valid utf-8 unicode.");
		ERR_FAIL_V(ERR_INVALID_DATA);
	}

	source = s;
#ifdef TOOLS_ENABLED
	source_changed_cache = true;
#endif
	//print_line("LSC :"+get_path());
	path = p_path;
	return OK;
}

const Map<StringName, GDScriptFunction *> &GDScript::debug_get_member_functions() const {

	return member_functions;
}

StringName GDScript::debug_get_member_by_index(int p_idx) const {

	for (const Map<StringName, MemberInfo>::Element *E = member_indices.front(); E; E = E->next()) {

		if (E->get().index == p_idx)
			return E->key();
	}

	return "<error>";
}

Ref<GDScript> GDScript::get_base() const {

	return base;
}

bool GDScript::has_script_signal(const StringName &p_signal) const {
	if (_signals.has(p_signal))
		return true;
	if (base.is_valid()) {
		return base->has_script_signal(p_signal);
	}
#ifdef TOOLS_ENABLED
	else if (base_cache.is_valid()) {
		return base_cache->has_script_signal(p_signal);
	}

#endif
	return false;
}
void GDScript::get_script_signal_list(List<MethodInfo> *r_signals) const {

	for (const Map<StringName, Vector<StringName> >::Element *E = _signals.front(); E; E = E->next()) {

		MethodInfo mi;
		mi.name = E->key();
		for (int i = 0; i < E->get().size(); i++) {
			PropertyInfo arg;
			arg.name = E->get()[i];
			mi.arguments.push_back(arg);
		}
		r_signals->push_back(mi);
	}

	if (base.is_valid()) {
		base->get_script_signal_list(r_signals);
	}
#ifdef TOOLS_ENABLED
	else if (base_cache.is_valid()) {
		base_cache->get_script_signal_list(r_signals);
	}

#endif
}

GDScript::GDScript() :
		script_list(this) {

	_static_ref = this;
	valid = false;
	subclass_count = 0;
	initializer = NULL;
	_base = NULL;
	_owner = NULL;
	tool = false;
#ifdef TOOLS_ENABLED
	source_changed_cache = false;
#endif

#ifdef DEBUG_ENABLED
	if (GDScriptLanguage::get_singleton()->lock) {
		GDScriptLanguage::get_singleton()->lock->lock();
	}
	GDScriptLanguage::get_singleton()->script_list.add(&script_list);

	if (GDScriptLanguage::get_singleton()->lock) {
		GDScriptLanguage::get_singleton()->lock->unlock();
	}
#endif
}

GDScript::~GDScript() {
	for (Map<StringName, GDScriptFunction *>::Element *E = member_functions.front(); E; E = E->next()) {
		memdelete(E->get());
	}

	for (Map<StringName, Ref<GDScript> >::Element *E = subclasses.front(); E; E = E->next()) {
		E->get()->_owner = NULL; //bye, you are no longer owned cause I died
	}

#ifdef DEBUG_ENABLED
	if (GDScriptLanguage::get_singleton()->lock) {
		GDScriptLanguage::get_singleton()->lock->lock();
	}
	GDScriptLanguage::get_singleton()->script_list.remove(&script_list);

	if (GDScriptLanguage::get_singleton()->lock) {
		GDScriptLanguage::get_singleton()->lock->unlock();
	}
#endif
}

//////////////////////////////
//         INSTANCE         //
//////////////////////////////

bool GDScriptInstance::set(const StringName &p_name, const Variant &p_value) {

	//member
	{
		const Map<StringName, GDScript::MemberInfo>::Element *E = script->member_indices.find(p_name);
		if (E) {
			if (E->get().setter) {
				const Variant *val = &p_value;
				Variant::CallError err;
				call(E->get().setter, &val, 1, err);
				if (err.error == Variant::CallError::CALL_OK) {
					return true; //function exists, call was successful
				}
			} else {
				if (!E->get().data_type.is_type(p_value)) {
					return false; // Type mismatch
				}
				members.write[E->get().index] = p_value;
			}
			return true;
		}
	}

	GDScript *sptr = script.ptr();
	while (sptr) {

		Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(GDScriptLanguage::get_singleton()->strings._set);
		if (E) {

			Variant name = p_name;
			const Variant *args[2] = { &name, &p_value };

			Variant::CallError err;
			Variant ret = E->get()->call(this, (const Variant **)args, 2, err);
			if (err.error == Variant::CallError::CALL_OK && ret.get_type() == Variant::BOOL && ret.operator bool())
				return true;
		}
		sptr = sptr->_base;
	}

	return false;
}

bool GDScriptInstance::get(const StringName &p_name, Variant &r_ret) const {

	const GDScript *sptr = script.ptr();
	while (sptr) {

		{
			const Map<StringName, GDScript::MemberInfo>::Element *E = script->member_indices.find(p_name);
			if (E) {
				if (E->get().getter) {
					Variant::CallError err;
					r_ret = const_cast<GDScriptInstance *>(this)->call(E->get().getter, NULL, 0, err);
					if (err.error == Variant::CallError::CALL_OK) {
						return true;
					}
				}
				r_ret = members[E->get().index];
				return true; //index found
			}
		}

		{

			const GDScript *sl = sptr;
			while (sl) {
				const Map<StringName, Variant>::Element *E = sl->constants.find(p_name);
				if (E) {
					r_ret = E->get();
					return true; //index found
				}
				sl = sl->_base;
			}
		}

		{
			const Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(GDScriptLanguage::get_singleton()->strings._get);
			if (E) {

				Variant name = p_name;
				const Variant *args[1] = { &name };

				Variant::CallError err;
				Variant ret = const_cast<GDScriptFunction *>(E->get())->call(const_cast<GDScriptInstance *>(this), (const Variant **)args, 1, err);
				if (err.error == Variant::CallError::CALL_OK && ret.get_type() != Variant::NIL) {
					r_ret = ret;
					return true;
				}
			}
		}
		sptr = sptr->_base;
	}

	return false;
}

Variant::Type GDScriptInstance::get_property_type(const StringName &p_name, bool *r_is_valid) const {

	const GDScript *sptr = script.ptr();
	while (sptr) {

		if (sptr->member_info.has(p_name)) {
			if (r_is_valid)
				*r_is_valid = true;
			return sptr->member_info[p_name].type;
		}
		sptr = sptr->_base;
	}

	if (r_is_valid)
		*r_is_valid = false;
	return Variant::NIL;
}

void GDScriptInstance::get_property_list(List<PropertyInfo> *p_properties) const {
	// exported members, not doen yet!

	const GDScript *sptr = script.ptr();
	List<PropertyInfo> props;

	while (sptr) {

		const Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(GDScriptLanguage::get_singleton()->strings._get_property_list);
		if (E) {

			Variant::CallError err;
			Variant ret = const_cast<GDScriptFunction *>(E->get())->call(const_cast<GDScriptInstance *>(this), NULL, 0, err);
			if (err.error == Variant::CallError::CALL_OK) {

				if (ret.get_type() != Variant::ARRAY) {

					ERR_EXPLAIN("Wrong type for _get_property list, must be an array of dictionaries.");
					ERR_FAIL();
				}
				Array arr = ret;
				for (int i = 0; i < arr.size(); i++) {

					Dictionary d = arr[i];
					ERR_CONTINUE(!d.has("name"));
					ERR_CONTINUE(!d.has("type"));
					PropertyInfo pinfo;
					pinfo.type = Variant::Type(d["type"].operator int());
					ERR_CONTINUE(pinfo.type < 0 || pinfo.type >= Variant::VARIANT_MAX);
					pinfo.name = d["name"];
					ERR_CONTINUE(pinfo.name == "");
					if (d.has("hint"))
						pinfo.hint = PropertyHint(d["hint"].operator int());
					if (d.has("hint_string"))
						pinfo.hint_string = d["hint_string"];
					if (d.has("usage"))
						pinfo.usage = d["usage"];

					props.push_back(pinfo);
				}
			}
		}

		//instance a fake script for editing the values

		Vector<_GDScriptMemberSort> msort;
		for (Map<StringName, PropertyInfo>::Element *E = sptr->member_info.front(); E; E = E->next()) {

			_GDScriptMemberSort ms;
			ERR_CONTINUE(!sptr->member_indices.has(E->key()));
			ms.index = sptr->member_indices[E->key()].index;
			ms.name = E->key();
			msort.push_back(ms);
		}

		msort.sort();
		msort.invert();
		for (int i = 0; i < msort.size(); i++) {

			props.push_front(sptr->member_info[msort[i].name]);
		}

		sptr = sptr->_base;
	}

	for (List<PropertyInfo>::Element *E = props.front(); E; E = E->next()) {

		p_properties->push_back(E->get());
	}
}

void GDScriptInstance::get_method_list(List<MethodInfo> *p_list) const {

	const GDScript *sptr = script.ptr();
	while (sptr) {

		for (Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.front(); E; E = E->next()) {

			MethodInfo mi;
			mi.name = E->key();
			mi.flags |= METHOD_FLAG_FROM_SCRIPT;
			for (int i = 0; i < E->get()->get_argument_count(); i++)
				mi.arguments.push_back(PropertyInfo(Variant::NIL, "arg" + itos(i)));
			p_list->push_back(mi);
		}
		sptr = sptr->_base;
	}
}

bool GDScriptInstance::has_method(const StringName &p_method) const {

	const GDScript *sptr = script.ptr();
	while (sptr) {
		const Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(p_method);
		if (E)
			return true;
		sptr = sptr->_base;
	}

	return false;
}
Variant GDScriptInstance::call(const StringName &p_method, const Variant **p_args, int p_argcount, Variant::CallError &r_error) {

	//printf("calling %ls:%i method %ls\n", script->get_path().c_str(), -1, String(p_method).c_str());

	GDScript *sptr = script.ptr();
	while (sptr) {
		Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(p_method);
		if (E) {
			return E->get()->call(this, p_args, p_argcount, r_error);
		}
		sptr = sptr->_base;
	}
	r_error.error = Variant::CallError::CALL_ERROR_INVALID_METHOD;
	return Variant();
}

void GDScriptInstance::call_multilevel(const StringName &p_method, const Variant **p_args, int p_argcount) {

	GDScript *sptr = script.ptr();
	Variant::CallError ce;

	while (sptr) {
		Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(p_method);
		if (E) {
			E->get()->call(this, p_args, p_argcount, ce);
		}
		sptr = sptr->_base;
	}
}

void GDScriptInstance::_ml_call_reversed(GDScript *sptr, const StringName &p_method, const Variant **p_args, int p_argcount) {

	if (sptr->_base)
		_ml_call_reversed(sptr->_base, p_method, p_args, p_argcount);

	Variant::CallError ce;

	Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(p_method);
	if (E) {
		E->get()->call(this, p_args, p_argcount, ce);
	}
}

void GDScriptInstance::call_multilevel_reversed(const StringName &p_method, const Variant **p_args, int p_argcount) {

	if (script.ptr()) {
		_ml_call_reversed(script.ptr(), p_method, p_args, p_argcount);
	}
}

void GDScriptInstance::notification(int p_notification) {

	//notification is not virtual, it gets called at ALL levels just like in C.
	Variant value = p_notification;
	const Variant *args[1] = { &value };

	GDScript *sptr = script.ptr();
	while (sptr) {
		Map<StringName, GDScriptFunction *>::Element *E = sptr->member_functions.find(GDScriptLanguage::get_singleton()->strings._notification);
		if (E) {
			Variant::CallError err;
			E->get()->call(this, args, 1, err);
			if (err.error != Variant::CallError::CALL_OK) {
				//print error about notification call
			}
		}
		sptr = sptr->_base;
	}
}

Ref<Script> GDScriptInstance::get_script() const {

	return script;
}

ScriptLanguage *GDScriptInstance::get_language() {

	return GDScriptLanguage::get_singleton();
}

MultiplayerAPI::RPCMode GDScriptInstance::get_rpc_mode(const StringName &p_method) const {

	const GDScript *cscript = script.ptr();

	while (cscript) {
		const Map<StringName, GDScriptFunction *>::Element *E = cscript->member_functions.find(p_method);
		if (E) {

			if (E->get()->get_rpc_mode() != MultiplayerAPI::RPC_MODE_DISABLED) {
				return E->get()->get_rpc_mode();
			}
		}
		cscript = cscript->_base;
	}

	return MultiplayerAPI::RPC_MODE_DISABLED;
}

MultiplayerAPI::RPCMode GDScriptInstance::get_rset_mode(const StringName &p_variable) const {

	const GDScript *cscript = script.ptr();

	while (cscript) {
		const Map<StringName, GDScript::MemberInfo>::Element *E = cscript->member_indices.find(p_variable);
		if (E) {

			if (E->get().rpc_mode) {
				return E->get().rpc_mode;
			}
		}
		cscript = cscript->_base;
	}

	return MultiplayerAPI::RPC_MODE_DISABLED;
}

void GDScriptInstance::reload_members() {

#ifdef DEBUG_ENABLED

	members.resize(script->member_indices.size()); //resize

	Vector<Variant> new_members;
	new_members.resize(script->member_indices.size());

	//pass the values to the new indices
	for (Map<StringName, GDScript::MemberInfo>::Element *E = script->member_indices.front(); E; E = E->next()) {

		if (member_indices_cache.has(E->key())) {
			Variant value = members[member_indices_cache[E->key()]];
			new_members.write[E->get().index] = value;
		}
	}

	//apply
	members = new_members;

	//pass the values to the new indices
	member_indices_cache.clear();
	for (Map<StringName, GDScript::MemberInfo>::Element *E = script->member_indices.front(); E; E = E->next()) {

		member_indices_cache[E->key()] = E->get().index;
	}

#endif
}

GDScriptInstance::GDScriptInstance() {
	owner = NULL;
	base_ref = false;
}

GDScriptInstance::~GDScriptInstance() {
	if (script.is_valid() && owner) {
#ifndef NO_THREADS
		GDScriptLanguage::singleton->lock->lock();
#endif

		script->instances.erase(owner);
#ifndef NO_THREADS
		GDScriptLanguage::singleton->lock->unlock();
#endif
	}
}

/************* SCRIPT LANGUAGE **************/

GDScriptLanguage *GDScriptLanguage::singleton = NULL;

String GDScriptLanguage::get_name() const {

	return "GDScript";
}

/* LANGUAGE FUNCTIONS */

void GDScriptLanguage::_add_global(const StringName &p_name, const Variant &p_value) {

	if (globals.has(p_name)) {
		//overwrite existing
		global_array.write[globals[p_name]] = p_value;
		return;
	}
	globals[p_name] = global_array.size();
	global_array.push_back(p_value);
	_global_array = global_array.ptrw();
}

void GDScriptLanguage::add_global_constant(const StringName &p_variable, const Variant &p_value) {

	_add_global(p_variable, p_value);
}

void GDScriptLanguage::add_named_global_constant(const StringName &p_name, const Variant &p_value) {
	named_globals[p_name] = p_value;
}

void GDScriptLanguage::remove_named_global_constant(const StringName &p_name) {
	ERR_FAIL_COND(!named_globals.has(p_name));
	named_globals.erase(p_name);
}

void GDScriptLanguage::init() {

	//populate global constants
	int gcc = GlobalConstants::get_global_constant_count();
	for (int i = 0; i < gcc; i++) {

		_add_global(StaticCString::create(GlobalConstants::get_global_constant_name(i)), GlobalConstants::get_global_constant_value(i));
	}

	_add_global(StaticCString::create("PI"), Math_PI);
	_add_global(StaticCString::create("TAU"), Math_TAU);
	_add_global(StaticCString::create("INF"), Math_INF);
	_add_global(StaticCString::create("NAN"), Math_NAN);

	//populate native classes

	List<StringName> class_list;
	ClassDB::get_class_list(&class_list);
	for (List<StringName>::Element *E = class_list.front(); E; E = E->next()) {

		StringName n = E->get();
		String s = String(n);
		if (s.begins_with("_"))
			n = s.substr(1, s.length());

		if (globals.has(n))
			continue;
		Ref<GDScriptNativeClass> nc = memnew(GDScriptNativeClass(E->get()));
		_add_global(n, nc);
	}

	//populate singletons

	List<Engine::Singleton> singletons;
	Engine::get_singleton()->get_singletons(&singletons);
	for (List<Engine::Singleton>::Element *E = singletons.front(); E; E = E->next()) {

		_add_global(E->get().name, E->get().ptr);
	}
}

String GDScriptLanguage::get_type() const {

	return "GDScript";
}
String GDScriptLanguage::get_extension() const {

	return "gd";
}
Error GDScriptLanguage::execute_file(const String &p_path) {

	// ??
	return OK;
}
void GDScriptLanguage::finish() {
}

void GDScriptLanguage::profiling_start() {

#ifdef DEBUG_ENABLED
	if (lock) {
		lock->lock();
	}

	SelfList<GDScriptFunction> *elem = function_list.first();
	while (elem) {
		elem->self()->profile.call_count = 0;
		elem->self()->profile.self_time = 0;
		elem->self()->profile.total_time = 0;
		elem->self()->profile.frame_call_count = 0;
		elem->self()->profile.frame_self_time = 0;
		elem->self()->profile.frame_total_time = 0;
		elem->self()->profile.last_frame_call_count = 0;
		elem->self()->profile.last_frame_self_time = 0;
		elem->self()->profile.last_frame_total_time = 0;
		elem = elem->next();
	}

	profiling = true;
	if (lock) {
		lock->unlock();
	}

#endif
}

void GDScriptLanguage::profiling_stop() {

#ifdef DEBUG_ENABLED
	if (lock) {
		lock->lock();
	}

	profiling = false;
	if (lock) {
		lock->unlock();
	}

#endif
}

int GDScriptLanguage::profiling_get_accumulated_data(ProfilingInfo *p_info_arr, int p_info_max) {

	int current = 0;
#ifdef DEBUG_ENABLED
	if (lock) {
		lock->lock();
	}

	SelfList<GDScriptFunction> *elem = function_list.first();
	while (elem) {
		if (current >= p_info_max)
			break;
		p_info_arr[current].call_count = elem->self()->profile.call_count;
		p_info_arr[current].self_time = elem->self()->profile.self_time;
		p_info_arr[current].total_time = elem->self()->profile.total_time;
		p_info_arr[current].signature = elem->self()->profile.signature;
		elem = elem->next();
		current++;
	}

	if (lock) {
		lock->unlock();
	}

#endif

	return current;
}

int GDScriptLanguage::profiling_get_frame_data(ProfilingInfo *p_info_arr, int p_info_max) {

	int current = 0;

#ifdef DEBUG_ENABLED
	if (lock) {
		lock->lock();
	}

	SelfList<GDScriptFunction> *elem = function_list.first();
	while (elem) {
		if (current >= p_info_max)
			break;
		if (elem->self()->profile.last_frame_call_count > 0) {
			p_info_arr[current].call_count = elem->self()->profile.last_frame_call_count;
			p_info_arr[current].self_time = elem->self()->profile.last_frame_self_time;
			p_info_arr[current].total_time = elem->self()->profile.last_frame_total_time;
			p_info_arr[current].signature = elem->self()->profile.signature;
			//print_line(String(elem->self()->profile.signature)+": "+itos(elem->self()->profile.last_frame_call_count));
			current++;
		}
		elem = elem->next();
	}

	if (lock) {
		lock->unlock();
	}

#endif

	return current;
}

struct GDScriptDepSort {

	//must support sorting so inheritance works properly (parent must be reloaded first)
	bool operator()(const Ref<GDScript> &A, const Ref<GDScript> &B) const {

		if (A == B)
			return false; //shouldn't happen but..
		const GDScript *I = B->get_base().ptr();
		while (I) {
			if (I == A.ptr()) {
				// A is a base of B
				return true;
			}

			I = I->get_base().ptr();
		}

		return false; //not a base
	}
};

void GDScriptLanguage::reload_all_scripts() {

#ifdef DEBUG_ENABLED
	print_line("RELOAD ALL SCRIPTS");
	if (lock) {
		lock->lock();
	}

	List<Ref<GDScript> > scripts;

	SelfList<GDScript> *elem = script_list.first();
	while (elem) {
		if (elem->self()->get_path().is_resource_file()) {
			print_line("FOUND: " + elem->self()->get_path());
			scripts.push_back(Ref<GDScript>(elem->self())); //cast to gdscript to avoid being erased by accident
		}
		elem = elem->next();
	}

	if (lock) {
		lock->unlock();
	}

	//as scripts are going to be reloaded, must proceed without locking here

	scripts.sort_custom<GDScriptDepSort>(); //update in inheritance dependency order

	for (List<Ref<GDScript> >::Element *E = scripts.front(); E; E = E->next()) {

		print_line("RELOADING: " + E->get()->get_path());
		E->get()->load_source_code(E->get()->get_path());
		E->get()->reload(true);
	}
#endif
}

void GDScriptLanguage::reload_tool_script(const Ref<Script> &p_script, bool p_soft_reload) {

#ifdef DEBUG_ENABLED

	if (lock) {
		lock->lock();
	}

	List<Ref<GDScript> > scripts;

	SelfList<GDScript> *elem = script_list.first();
	while (elem) {
		if (elem->self()->get_path().is_resource_file()) {

			scripts.push_back(Ref<GDScript>(elem->self())); //cast to gdscript to avoid being erased by accident
		}
		elem = elem->next();
	}

	if (lock) {
		lock->unlock();
	}

	//when someone asks you why dynamically typed languages are easier to write....

	Map<Ref<GDScript>, Map<ObjectID, List<Pair<StringName, Variant> > > > to_reload;

	//as scripts are going to be reloaded, must proceed without locking here

	scripts.sort_custom<GDScriptDepSort>(); //update in inheritance dependency order

	for (List<Ref<GDScript> >::Element *E = scripts.front(); E; E = E->next()) {

		bool reload = E->get() == p_script || to_reload.has(E->get()->get_base());

		if (!reload)
			continue;

		to_reload.insert(E->get(), Map<ObjectID, List<Pair<StringName, Variant> > >());

		if (!p_soft_reload) {

			//save state and remove script from instances
			Map<ObjectID, List<Pair<StringName, Variant> > > &map = to_reload[E->get()];

			while (E->get()->instances.front()) {
				Object *obj = E->get()->instances.front()->get();
				//save instance info
				List<Pair<StringName, Variant> > state;
				if (obj->get_script_instance()) {

					obj->get_script_instance()->get_property_state(state);
					map[obj->get_instance_id()] = state;
					obj->set_script(RefPtr());
				}
			}

//same thing for placeholders
#ifdef TOOLS_ENABLED

			while (E->get()->placeholders.size()) {
				Object *obj = E->get()->placeholders.front()->get()->get_owner();

				//save instance info
				if (obj->get_script_instance()) {

					map.insert(obj->get_instance_id(), List<Pair<StringName, Variant> >());
					List<Pair<StringName, Variant> > &state = map[obj->get_instance_id()];
					obj->get_script_instance()->get_property_state(state);
					obj->set_script(RefPtr());
				} else {
					// no instance found. Let's remove it so we don't loop forever
					E->get()->placeholders.erase(E->get()->placeholders.front()->get());
				}
			}

#endif

			for (Map<ObjectID, List<Pair<StringName, Variant> > >::Element *F = E->get()->pending_reload_state.front(); F; F = F->next()) {
				map[F->key()] = F->get(); //pending to reload, use this one instead
			}
		}
	}

	for (Map<Ref<GDScript>, Map<ObjectID, List<Pair<StringName, Variant> > > >::Element *E = to_reload.front(); E; E = E->next()) {

		Ref<GDScript> scr = E->key();
		scr->reload(p_soft_reload);

		//restore state if saved
		for (Map<ObjectID, List<Pair<StringName, Variant> > >::Element *F = E->get().front(); F; F = F->next()) {

			Object *obj = ObjectDB::get_instance(F->key());
			if (!obj)
				continue;

			if (!p_soft_reload) {
				//clear it just in case (may be a pending reload state)
				obj->set_script(RefPtr());
			}
			obj->set_script(scr.get_ref_ptr());
			if (!obj->get_script_instance()) {
				//failed, save reload state for next time if not saved
				if (!scr->pending_reload_state.has(obj->get_instance_id())) {
					scr->pending_reload_state[obj->get_instance_id()] = F->get();
				}
				continue;
			}

			for (List<Pair<StringName, Variant> >::Element *G = F->get().front(); G; G = G->next()) {
				obj->get_script_instance()->set(G->get().first, G->get().second);
			}

			scr->pending_reload_state.erase(obj->get_instance_id()); //as it reloaded, remove pending state
		}

		//if instance states were saved, set them!
	}

#endif
}

void GDScriptLanguage::frame() {

	//print_line("calls: "+itos(calls));
	calls = 0;

#ifdef DEBUG_ENABLED
	if (profiling) {
		if (lock) {
			lock->lock();
		}

		SelfList<GDScriptFunction> *elem = function_list.first();
		while (elem) {
			elem->self()->profile.last_frame_call_count = elem->self()->profile.frame_call_count;
			elem->self()->profile.last_frame_self_time = elem->self()->profile.frame_self_time;
			elem->self()->profile.last_frame_total_time = elem->self()->profile.frame_total_time;
			elem->self()->profile.frame_call_count = 0;
			elem->self()->profile.frame_self_time = 0;
			elem->self()->profile.frame_total_time = 0;
			elem = elem->next();
		}

		if (lock) {
			lock->unlock();
		}
	}

#endif
}

/* EDITOR FUNCTIONS */
void GDScriptLanguage::get_reserved_words(List<String> *p_words) const {

	static const char *_reserved_words[] = {
		// operators
		"and",
		"in",
		"not",
		"or",
		// types and values
		"false",
		"float",
		"int",
		"bool",
		"null",
		"PI",
		"TAU",
		"INF",
		"NAN",
		"self",
		"true",
		"void",
		// functions
		"as",
		"assert",
		"breakpoint",
		"class",
		"class_name",
		"extends",
		"is",
		"func",
		"preload",
		"setget",
		"signal",
		"tool",
		"yield",
		// var
		"const",
		"enum",
		"export",
		"onready",
		"static",
		"var",
		// control flow
		"break",
		"continue",
		"if",
		"elif",
		"else",
		"for",
		"pass",
		"return",
		"match",
		"while",
		"remote",
		"sync",
		"master",
		"slave",
		"remotesync",
		"mastersync",
		"slavesync",
		0
	};

	const char **w = _reserved_words;

	while (*w) {

		p_words->push_back(*w);
		w++;
	}

	for (int i = 0; i < GDScriptFunctions::FUNC_MAX; i++) {
		p_words->push_back(GDScriptFunctions::get_func_name(GDScriptFunctions::Function(i)));
	}
}

bool GDScriptLanguage::handles_global_class_type(const String &p_type) const {

	return p_type == "GDScript";
}

String GDScriptLanguage::get_global_class_name(const String &p_path, String *r_base_type) const {

	PoolVector<uint8_t> sourcef;
	Error err;
	FileAccess *f = FileAccess::open(p_path, FileAccess::READ, &err);
	if (err) {
		return String();
	}

	int len = f->get_len();
	sourcef.resize(len + 1);
	PoolVector<uint8_t>::Write w = sourcef.write();
	int r = f->get_buffer(w.ptr(), len);
	f->close();
	memdelete(f);
	ERR_FAIL_COND_V(r != len, String());
	w[len] = 0;

	String s;
	if (s.parse_utf8((const char *)w.ptr())) {
		return String();
	}

	GDScriptParser parser;

	parser.parse(s, p_path.get_base_dir(), true, p_path);

	if (parser.get_parse_tree() && parser.get_parse_tree()->type == GDScriptParser::Node::TYPE_CLASS) {

		const GDScriptParser::ClassNode *c = static_cast<const GDScriptParser::ClassNode *>(parser.get_parse_tree());
		if (r_base_type) {
			GDScriptParser::DataType base_type;
			if (c->base_type.has_type) {
				base_type = c->base_type;
				while (base_type.has_type && base_type.kind != GDScriptParser::DataType::NATIVE) {
					switch (base_type.kind) {
						case GDScriptParser::DataType::CLASS: {
							base_type = base_type.class_type->base_type;
						} break;
						case GDScriptParser::DataType::GDSCRIPT: {
							Ref<GDScript> gds = base_type.script_type;
							if (gds.is_valid()) {
								base_type.kind = GDScriptParser::DataType::NATIVE;
								base_type.native_type = gds->get_instance_base_type();
							} else {
								base_type = GDScriptParser::DataType();
							}
						} break;
						default: {
							base_type = GDScriptParser::DataType();
						} break;
					}
				}
			}
			if (base_type.has_type) {
				*r_base_type = base_type.native_type;
			} else {
				// Fallback
				if (c->extends_used && c->extends_class.size() == 1) {
					*r_base_type = c->extends_class[0];
				} else if (!c->extends_used) {
					*r_base_type = "Reference";
				}
			}
		}
		return c->name;
	}

	return String();
}

#ifdef DEBUG_ENABLED
String GDScriptWarning::get_message() const {

#define CHECK_SYMBOLS(m_amount) ERR_FAIL_COND_V(symbols.size() < m_amount, String());

	switch (code) {
		case UNASSIGNED_VARIABLE_OP_ASSIGN: {
			CHECK_SYMBOLS(1);
			return "Using assignment with operation but the variable '" + symbols[0] + "' was not previously assigned a value.";
		} break;
		case UNASSIGNED_VARIABLE: {
			CHECK_SYMBOLS(1);
			return "The variable '" + symbols[0] + "' was used but never assigned a value.";
		} break;
		case UNUSED_VARIABLE: {
			CHECK_SYMBOLS(1);
			return "The local variable '" + symbols[0] + "' is declared but never used in the block.";
		} break;
		case UNUSED_CLASS_VARIABLE: {
			CHECK_SYMBOLS(1);
			return "The class variable '" + symbols[0] + "' is declared but never used in the script.";
		} break;
		case UNUSED_ARGUMENT: {
			CHECK_SYMBOLS(2);
			return "The argument '" + symbols[1] + "' is never used in the function '" + symbols[0] + "'.";
		} break;
		case UNREACHABLE_CODE: {
			CHECK_SYMBOLS(1);
			return "Unreachable code (statement after return) in function '" + symbols[0] + "()'.";
		} break;
		case STANDALONE_EXPRESSION: {
			return "Standalone expression (the line has no effect).";
		} break;
		case VOID_ASSIGNMENT: {
			CHECK_SYMBOLS(1);
			return "Assignment operation, but the function '" + symbols[0] + "()' returns void.";
		} break;
		case NARROWING_CONVERSION: {
			return "Narrowing coversion (float is converted to int and lose precision).";
		} break;
		case FUNCTION_MAY_YIELD: {
			CHECK_SYMBOLS(1);
			return "Assigned variable is typed but the function '" + symbols[0] + "()' may yield and return a GDScriptFunctionState instead.";
		} break;
		case VARIABLE_CONFLICTS_FUNCTION: {
			CHECK_SYMBOLS(1);
			return "Variable declaration of '" + symbols[0] + "' conflicts with a function of the same name.";
		} break;
		case FUNCTION_CONFLICTS_VARIABLE: {
			CHECK_SYMBOLS(1);
			return "Function declaration of '" + symbols[0] + "()' conflicts with a variable of the same name.";
		} break;
		case FUNCTION_CONFLICTS_CONSTANT: {
			CHECK_SYMBOLS(1);
			return "Function declaration of '" + symbols[0] + "()' conflicts with a constant of the same name.";
		} break;
		case INCOMPATIBLE_TERNARY: {
			return "Values of the ternary conditional are not mutually compatible.";
		} break;
		case UNUSED_SIGNAL: {
			CHECK_SYMBOLS(1);
			return "The signal '" + symbols[0] + "' is declared but never emitted.";
		} break;
		case RETURN_VALUE_DISCARDED: {
			CHECK_SYMBOLS(1);
			return "The function '" + symbols[0] + "()' returns a value, but this value is never used.";
		} break;
		case PROPERTY_USED_AS_FUNCTION: {
			CHECK_SYMBOLS(2);
			return "The method '" + symbols[0] + "()' was not found in base '" + symbols[1] + "' but there's a property with the same name. Did you mean to access it?";
		} break;
		case CONSTANT_USED_AS_FUNCTION: {
			CHECK_SYMBOLS(2);
			return "The method '" + symbols[0] + "()' was not found in base '" + symbols[1] + "' but there's a constant with the same name. Did you mean to access it?";
		} break;
		case FUNCTION_USED_AS_PROPERTY: {
			CHECK_SYMBOLS(2);
			return "The property '" + symbols[0] + "' was not found in base '" + symbols[1] + "' but there's a method with the same name. Did you mean to call it?";
		} break;
		case INTEGER_DIVISION: {
			return "Integer division, decimal part will be discarded.";
		} break;
		case UNSAFE_PROPERTY_ACCESS: {
			CHECK_SYMBOLS(2);
			return "The property '" + symbols[0] + "' is not present on the inferred type '" + symbols[1] + "' (but may be present on a subtype).";
		} break;
		case UNSAFE_METHOD_ACCESS: {
			CHECK_SYMBOLS(2);
			return "The method '" + symbols[0] + "' is not present on the inferred type '" + symbols[1] + "' (but may be present on a subtype).";
		} break;
		case UNSAFE_CAST: {
			CHECK_SYMBOLS(1);
			return "The value is cast to '" + symbols[0] + "' but has an unkown type.";
		} break;
		case UNSAFE_CALL_ARGUMENT: {
			CHECK_SYMBOLS(4);
			return "The argument '" + symbols[0] + "' of the function '" + symbols[1] + "' requires a the subtype '" + symbols[2] + "' but the supertype '" + symbols[3] + "' was provided";
		} break;
	}
	ERR_EXPLAIN("Invalid GDScript waring code: " + get_name_from_code(code));
	ERR_FAIL_V(String());

#undef CHECK_SYMBOLS
}

String GDScriptWarning::get_name() const {
	return get_name_from_code(code);
}

String GDScriptWarning::get_name_from_code(Code p_code) {
	ERR_FAIL_COND_V(p_code < 0 || p_code >= WARNING_MAX, String());

	static const char *names[] = {
		"UNASSIGNED_VARIABLE",
		"UNASSIGNED_VARIABLE_OP_ASSIGN",
		"UNUSED_VARIABLE",
		"UNUSED_CLASS_VARIABLE",
		"UNUSED_ARGUMENT",
		"UNREACHABLE_CODE",
		"STANDALONE_EXPRESSION",
		"VOID_ASSIGNMENT",
		"NARROWING_CONVERSION",
		"FUNCTION_MAY_YIELD",
		"VARIABLE_CONFLICTS_FUNCTION",
		"FUNCTION_CONFLICTS_VARIABLE",
		"FUNCTION_CONFLICTS_CONSTANT",
		"INCOMPATIBLE_TERNARY",
		"UNUSED_SIGNAL",
		"RETURN_VALUE_DISCARDED",
		"PROPERTY_USED_AS_FUNCTION",
		"CONSTANT_USED_AS_FUNCTION",
		"FUNCTION_USED_AS_PROPERTY",
		"INTEGER_DIVISION",
		"UNSAFE_PROPERTY_ACCESS",
		"UNSAFE_METHOD_ACCESS",
		"UNSAFE_CAST",
		"UNSAFE_CALL_ARGUMENT",
		NULL
	};

	return names[(int)p_code];
}

GDScriptWarning::Code GDScriptWarning::get_code_from_name(const String &p_name) {
	for (int i = 0; i < WARNING_MAX; i++) {
		if (get_name_from_code((Code)i) == p_name) {
			return (Code)i;
		}
	}

	ERR_EXPLAIN("Invalid GDScript waring name: " + p_name);
	ERR_FAIL_V(WARNING_MAX);
}

#endif // DEBUG_ENABLED

GDScriptLanguage::GDScriptLanguage() {

	calls = 0;
	ERR_FAIL_COND(singleton);
	singleton = this;
	strings._init = StaticCString::create("_init");
	strings._notification = StaticCString::create("_notification");
	strings._set = StaticCString::create("_set");
	strings._get = StaticCString::create("_get");
	strings._get_property_list = StaticCString::create("_get_property_list");
	strings._script_source = StaticCString::create("script/source");
	_debug_parse_err_line = -1;
	_debug_parse_err_file = "";

#ifdef NO_THREADS
	lock = NULL;
#else
	lock = Mutex::create();
#endif
	profiling = false;
	script_frame_time = 0;

	_debug_call_stack_pos = 0;
	int dmcs = GLOBAL_DEF("debug/settings/gdscript/max_call_stack", 1024);
	if (ScriptDebugger::get_singleton()) {
		//debugging enabled!

		_debug_max_call_stack = dmcs;
		if (_debug_max_call_stack < 1024)
			_debug_max_call_stack = 1024;
		_call_stack = memnew_arr(CallLevel, _debug_max_call_stack + 1);

	} else {
		_debug_max_call_stack = 0;
		_call_stack = NULL;
	}

#ifdef DEBUG_ENABLED
	GLOBAL_DEF("debug/gdscript/warnings/enable", true);
	GLOBAL_DEF("debug/gdscript/warnings/treat_warnings_as_errors", false);
	for (int i = 0; i < (int)GDScriptWarning::WARNING_MAX; i++) {
		String warning = GDScriptWarning::get_name_from_code((GDScriptWarning::Code)i).to_lower();
		GLOBAL_DEF("debug/gdscript/warnings/" + warning, !warning.begins_with("unsafe_"));
	}
#endif // DEBUG_ENABLED
}

GDScriptLanguage::~GDScriptLanguage() {

	if (lock) {
		memdelete(lock);
		lock = NULL;
	}
	if (_call_stack) {
		memdelete_arr(_call_stack);
	}
	singleton = NULL;
}

/*************** RESOURCE ***************/

RES ResourceFormatLoaderGDScript::load(const String &p_path, const String &p_original_path, Error *r_error) {

	if (r_error)
		*r_error = ERR_FILE_CANT_OPEN;

	GDScript *script = memnew(GDScript);

	Ref<GDScript> scriptres(script);

	if (p_path.ends_with(".gde") || p_path.ends_with(".gdc")) {

		script->set_script_path(p_original_path); // script needs this.
		script->set_path(p_original_path);
		Error err = script->load_byte_code(p_path);

		if (err != OK) {

			ERR_FAIL_COND_V(err != OK, RES());
		}

	} else {
		Error err = script->load_source_code(p_path);

		if (err != OK) {

			ERR_FAIL_COND_V(err != OK, RES());
		}

		script->set_script_path(p_original_path); // script needs this.
		script->set_path(p_original_path);
		//script->set_name(p_path.get_file());

		script->reload();
	}
	if (r_error)
		*r_error = OK;

	return scriptres;
}
void ResourceFormatLoaderGDScript::get_recognized_extensions(List<String> *p_extensions) const {

	p_extensions->push_back("gd");
	p_extensions->push_back("gdc");
	p_extensions->push_back("gde");
}

bool ResourceFormatLoaderGDScript::handles_type(const String &p_type) const {

	return (p_type == "Script" || p_type == "GDScript");
}

String ResourceFormatLoaderGDScript::get_resource_type(const String &p_path) const {

	String el = p_path.get_extension().to_lower();
	if (el == "gd" || el == "gdc" || el == "gde")
		return "GDScript";
	return "";
}

Error ResourceFormatSaverGDScript::save(const String &p_path, const RES &p_resource, uint32_t p_flags) {

	Ref<GDScript> sqscr = p_resource;
	ERR_FAIL_COND_V(sqscr.is_null(), ERR_INVALID_PARAMETER);

	String source = sqscr->get_source_code();

	Error err;
	FileAccess *file = FileAccess::open(p_path, FileAccess::WRITE, &err);

	if (err) {

		ERR_FAIL_COND_V(err, err);
	}

	file->store_string(source);
	if (file->get_error() != OK && file->get_error() != ERR_FILE_EOF) {
		memdelete(file);
		return ERR_CANT_CREATE;
	}
	file->close();
	memdelete(file);

	if (ScriptServer::is_reload_scripts_on_save_enabled()) {
		GDScriptLanguage::get_singleton()->reload_tool_script(p_resource, false);
	}

	return OK;
}

void ResourceFormatSaverGDScript::get_recognized_extensions(const RES &p_resource, List<String> *p_extensions) const {

	if (Object::cast_to<GDScript>(*p_resource)) {
		p_extensions->push_back("gd");
	}
}
bool ResourceFormatSaverGDScript::recognize(const RES &p_resource) const {

	return Object::cast_to<GDScript>(*p_resource) != NULL;
}
