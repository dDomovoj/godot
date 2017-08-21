/*************************************************************************/
/*  arvr_server.cpp                                                      */
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

#include "arvr_server.h"
#include "arvr/arvr_interface.h"
#include "arvr/arvr_positional_tracker.h"
#include "project_settings.h"

ARVRServer *ARVRServer::singleton = NULL;

ARVRServer *ARVRServer::get_singleton() {
	return singleton;
};

void ARVRServer::_bind_methods() {
	ClassDB::bind_method(D_METHOD("get_world_scale"), &ARVRServer::get_world_scale);
	ClassDB::bind_method(D_METHOD("set_world_scale"), &ARVRServer::set_world_scale);
	ClassDB::bind_method(D_METHOD("get_reference_frame"), &ARVRServer::get_reference_frame);
	ClassDB::bind_method(D_METHOD("request_reference_frame", "ignore_tilt", "keep_height"), &ARVRServer::request_reference_frame);

	ADD_PROPERTY(PropertyInfo(Variant::REAL, "world_scale"), "set_world_scale", "get_world_scale");

	ClassDB::bind_method(D_METHOD("get_interface_count"), &ARVRServer::get_interface_count);
	ClassDB::bind_method(D_METHOD("get_interface", "idx"), &ARVRServer::get_interface);
	ClassDB::bind_method(D_METHOD("find_interface", "name"), &ARVRServer::find_interface);
	ClassDB::bind_method(D_METHOD("get_tracker_count"), &ARVRServer::get_tracker_count);
	ClassDB::bind_method(D_METHOD("get_tracker", "idx"), &ARVRServer::get_tracker);

	ClassDB::bind_method(D_METHOD("set_primary_interface"), &ARVRServer::set_primary_interface);

	ClassDB::bind_method(D_METHOD("add_interface"), &ARVRServer::add_interface);
	ClassDB::bind_method(D_METHOD("remove_interface"), &ARVRServer::remove_interface);

	BIND_ENUM_CONSTANT(TRACKER_CONTROLLER);
	BIND_ENUM_CONSTANT(TRACKER_BASESTATION);
	BIND_ENUM_CONSTANT(TRACKER_ANCHOR);
	BIND_ENUM_CONSTANT(TRACKER_UNKNOWN);
	BIND_ENUM_CONSTANT(TRACKER_ANY_KNOWN);
	BIND_ENUM_CONSTANT(TRACKER_ANY);

	ADD_SIGNAL(MethodInfo("interface_added", PropertyInfo(Variant::STRING, "name")));
	ADD_SIGNAL(MethodInfo("interface_removed", PropertyInfo(Variant::STRING, "name")));

	ADD_SIGNAL(MethodInfo("tracker_added", PropertyInfo(Variant::STRING, "name"), PropertyInfo(Variant::INT, "type")));
	ADD_SIGNAL(MethodInfo("tracker_removed", PropertyInfo(Variant::STRING, "name")));
};

real_t ARVRServer::get_world_scale() const {
	return world_scale;
};

void ARVRServer::set_world_scale(real_t p_world_scale) {
	if (world_scale < 0.01) {
		world_scale = 0.01;
	} else if (world_scale > 1000.0) {
		world_scale = 1000.0;
	};

	world_scale = p_world_scale;
};

Transform ARVRServer::get_world_origin() const {
	return world_origin;
};

void ARVRServer::set_world_origin(const Transform p_world_origin) {
	world_origin = p_world_origin;
};

Transform ARVRServer::get_reference_frame() const {
	return reference_frame;
};

void ARVRServer::request_reference_frame(bool p_ignore_tilt, bool p_keep_height) {
	if (primary_interface != NULL) {
		// clear our current reference frame or we'll end up double adjusting it
		reference_frame = Transform();

		// requesting our EYE_MONO transform should return our current HMD position
		Transform new_reference_frame = primary_interface->get_transform_for_eye(ARVRInterface::EYE_MONO, Transform());

		// remove our tilt
		if (p_ignore_tilt) {
			// take the Y out of our Z
			new_reference_frame.basis.set_axis(2, Vector3(new_reference_frame.basis.elements[0][2], 0.0, new_reference_frame.basis.elements[2][2]).normalized());

			// Y is straight up
			new_reference_frame.basis.set_axis(1, Vector3(0.0, 1.0, 0.0));

			// and X is our cross reference
			new_reference_frame.basis.set_axis(0, new_reference_frame.basis.get_axis(1).cross(new_reference_frame.basis.get_axis(2)).normalized());
		};

		// don't negate our height
		if (p_keep_height) {
			new_reference_frame.origin.y = 0.0;
		};

		reference_frame = new_reference_frame.inverse();
	};
};

void ARVRServer::add_interface(const Ref<ARVRInterface> &p_interface) {
	ERR_FAIL_COND(p_interface.is_null());

	int idx = -1;
	for (int i = 0; i < interfaces.size(); i++) {

		if (interfaces[i] == p_interface) {
			ERR_PRINT("Interface was already added");
			return;
		};
	};

	print_line("Registered interface " + p_interface->get_name());

	interfaces.push_back(p_interface);
	emit_signal("interface_added", p_interface->get_name());
};

void ARVRServer::remove_interface(const Ref<ARVRInterface> &p_interface) {
	ERR_FAIL_COND(p_interface.is_null());

	int idx = -1;
	for (int i = 0; i < interfaces.size(); i++) {

		if (interfaces[i] == p_interface) {

			idx = i;
			break;
		};
	};

	ERR_FAIL_COND(idx == -1);

	print_line("Removed interface" + p_interface->get_name());

	emit_signal("interface_removed", p_interface->get_name());
	interfaces.remove(idx);
};

int ARVRServer::get_interface_count() const {
	return interfaces.size();
};

Ref<ARVRInterface> ARVRServer::get_interface(int p_index) const {
	ERR_FAIL_INDEX_V(p_index, interfaces.size(), NULL);

	return interfaces[p_index];
};

Ref<ARVRInterface> ARVRServer::find_interface(const String &p_name) const {
	int idx = -1;
	for (int i = 0; i < interfaces.size(); i++) {

		if (interfaces[i]->get_name() == p_name) {

			idx = i;
			break;
		};
	};

	ERR_FAIL_COND_V(idx == -1, NULL);

	return interfaces[idx];
};

/*
	A little extra info on the tracker ids, these are unique per tracker type so we get soem consistency in recognising our trackers, specifically controllers.

	The first controller that is turned of will get ID 1, the second will get ID 2, etc.
	The magic happens when one of the controllers is turned off, say controller 1 turns off, controller 2 will remain controller 2, controller 3 will remain controller 3.
	If controller number 1 is turned on again it again gets ID 1 unless another new controller was turned on since.

	The most likely scenario however is a controller that runs out of battery and another controller being used to replace it.
	Because the controllers are often linked to physical objects, say you're holding a shield in controller 1, your left hand, and a gun in controller 2, your right hand, and controller 1 dies:
	- using our tracker index would suddenly make the gun disappear and the shield jump into your right hand because controller 2 becomes controller 1.
	- using this approach the shield disappears or is no longer tracked, but the gun stays firmly in your right hand because that is still controller 2, further more, if controller 1 is replaced the shield will return.
*/

bool ARVRServer::is_tracker_id_in_use_for_type(TrackerType p_tracker_type, int p_tracker_id) const {
	for (int i = 0; i < trackers.size(); i++) {
		if (trackers[i]->get_type() == p_tracker_type && trackers[i]->get_tracker_id() == p_tracker_id) {
			return true;
		};
	};

	// all good
	return false;
};

int ARVRServer::get_free_tracker_id_for_type(TrackerType p_tracker_type) {
	// we start checking at 1, 0 means that it's not a controller..
	int tracker_id = 1;

	while (is_tracker_id_in_use_for_type(p_tracker_type, tracker_id)) {
		// try the next one
		tracker_id++;
	};

	return tracker_id;
};

void ARVRServer::add_tracker(ARVRPositionalTracker *p_tracker) {
	ERR_FAIL_NULL(p_tracker);

	trackers.push_back(p_tracker);
	emit_signal("tracker_added", p_tracker->get_name(), p_tracker->get_type());
};

void ARVRServer::remove_tracker(ARVRPositionalTracker *p_tracker) {
	ERR_FAIL_NULL(p_tracker);

	int idx = -1;
	for (int i = 0; i < trackers.size(); i++) {

		if (trackers[i] == p_tracker) {

			idx = i;
			break;
		};
	};

	ERR_FAIL_COND(idx == -1);

	emit_signal("tracker_removed", p_tracker->get_name());
	trackers.remove(idx);
};

int ARVRServer::get_tracker_count() const {
	return trackers.size();
};

ARVRPositionalTracker *ARVRServer::get_tracker(int p_index) const {
	ERR_FAIL_INDEX_V(p_index, trackers.size(), NULL);

	return trackers[p_index];
};

ARVRPositionalTracker *ARVRServer::find_by_type_and_id(TrackerType p_tracker_type, int p_tracker_id) const {
	ERR_FAIL_COND_V(p_tracker_id == 0, NULL);

	for (int i = 0; i < trackers.size(); i++) {
		if (trackers[i]->get_type() == p_tracker_type && trackers[i]->get_tracker_id() == p_tracker_id) {
			return trackers[i];
		};
	};

	return NULL;
};

Ref<ARVRInterface> ARVRServer::get_primary_interface() const {
	return primary_interface;
};

void ARVRServer::set_primary_interface(const Ref<ARVRInterface> &p_primary_interface) {
	primary_interface = p_primary_interface;

	print_line("Primary interface set to: " + primary_interface->get_name());
};

void ARVRServer::clear_primary_interface_if(const Ref<ARVRInterface> &p_primary_interface) {
	if (primary_interface == p_primary_interface) {
		print_line("Clearing primary interface");
		primary_interface.unref();
	};
};

ARVRServer::ARVRServer() {
	singleton = this;
	world_scale = 1.0;
};

ARVRServer::~ARVRServer() {
	primary_interface.unref();

	while (interfaces.size() > 0) {
		interfaces.remove(0);
	}

	while (trackers.size() > 0) {
		trackers.remove(0);
	}

	singleton = NULL;
};
