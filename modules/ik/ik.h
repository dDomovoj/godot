/*************************************************************************/
/*  ik.h                                                           */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                    http://www.godotengine.org                         */
/*************************************************************************/
/* Copyright (c) 2007-2016 Juan Linietsky, Ariel Manzur.                 */
/* This file is (c) 2016 Sergey Lapin <slapinid@gmail.com>               */
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
#ifndef IK_H
#define IK_H

#include "scene/3d/skeleton.h"
class InverseKinematics : public Spatial {
	OBJ_TYPE(InverseKinematics, Spatial);
	bool bound;
	String ik_bone;
	int ik_bone_no;
	int tail_bone;
	int chain_size;
	Skeleton *skel;
	List<int> chain;
	void _check_bind();
	void _check_unbind();
	int iterations;
	float precision;
	float speed;

protected:
	bool _set(const StringName& p_name, const Variant& p_value);
	bool _get(const StringName& p_name,Variant &r_ret) const;
	void _get_property_list( List<PropertyInfo> *p_list) const;

	void _notification(int p_what);
	static void _bind_methods();
public:
	Skeleton *get_skeleton();
	void set_bone_name(const String& p_name);
	String get_bone_name() const;
	void set_iterations(int itn);
	int get_iterations() const;
	void set_chain_size(int cs);
	int get_chain_size() const;
	void set_precision(float p);
	float get_precision() const;
	void set_speed(float p);
	float get_speed() const;
	InverseKinematics();
};

#endif

