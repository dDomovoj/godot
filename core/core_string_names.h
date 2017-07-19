/*************************************************************************/
/*  core_string_names.h                                                  */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                    http://www.godotengine.org                         */
/*************************************************************************/
/* Copyright (c) 2007-2017 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2017 Godot Engine contributors (cf. AUTHORS.md)    */
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
#ifndef CORE_STRING_NAMES_H
#define CORE_STRING_NAMES_H

#include "string_db.h"

class CoreStringNames {

	friend void register_core_types();
	friend void unregister_core_types();

	static CoreStringNames *singleton;

	static void create() { singleton = memnew(CoreStringNames); }
	static void free() {
		memdelete(singleton);
		singleton = NULL;
	}

	CoreStringNames();

public:
	_FORCE_INLINE_ static CoreStringNames *get_singleton() { return singleton; }

	StringName _free;
	StringName changed;
	StringName _meta;
	StringName _script;
	StringName script_changed;
	StringName ___pdcdata;
	StringName __getvar;
	StringName _iter_init;
	StringName _iter_next;
	StringName _iter_get;
	StringName get_rid;
	StringName _custom_features;
#ifdef TOOLS_ENABLED
	StringName _sections_unfolded;
#endif
};

#endif // SCENE_STRING_NAMES_H
