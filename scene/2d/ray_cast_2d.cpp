/*************************************************************************/
/*  ray_cast_2d.cpp                                                      */
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
#include "ray_cast_2d.h"
#include "collision_object_2d.h"
#include "physics_body_2d.h"
#include "servers/physics_2d_server.h"

void RayCast2D::set_cast_to(const Vector2 &p_point) {

	cast_to = p_point;
	if (is_inside_tree() && (get_tree()->is_editor_hint() || get_tree()->is_debugging_collisions_hint()))
		update();
}

Vector2 RayCast2D::get_cast_to() const {

	return cast_to;
}

void RayCast2D::set_layer_mask(uint32_t p_mask) {

	layer_mask = p_mask;
}

uint32_t RayCast2D::get_layer_mask() const {

	return layer_mask;
}

void RayCast2D::set_type_mask(uint32_t p_mask) {

	type_mask = p_mask;
}

uint32_t RayCast2D::get_type_mask() const {

	return type_mask;
}

bool RayCast2D::is_colliding() const {

	return collided;
}
Object *RayCast2D::get_collider() const {

	if (against == 0)
		return NULL;

	return ObjectDB::get_instance(against);
}

int RayCast2D::get_collider_shape() const {

	return against_shape;
}
Vector2 RayCast2D::get_collision_point() const {

	return collision_point;
}
Vector2 RayCast2D::get_collision_normal() const {

	return collision_normal;
}

void RayCast2D::set_enabled(bool p_enabled) {

	enabled = p_enabled;
	if (is_inside_tree() && !get_tree()->is_editor_hint())
		set_fixed_process(p_enabled);
	if (!p_enabled)
		collided = false;
}

bool RayCast2D::is_enabled() const {

	return enabled;
}

void RayCast2D::set_exclude_parent_body(bool p_exclude_parent_body) {

	if (exclude_parent_body == p_exclude_parent_body)
		return;

	exclude_parent_body = p_exclude_parent_body;

	if (!is_inside_tree())
		return;

	if (get_parent()->cast_to<PhysicsBody2D>()) {
		if (exclude_parent_body)
			exclude.insert(get_parent()->cast_to<PhysicsBody2D>()->get_rid());
		else
			exclude.erase(get_parent()->cast_to<PhysicsBody2D>()->get_rid());
	}
}

bool RayCast2D::get_exclude_parent_body() const {

	return exclude_parent_body;
}

void RayCast2D::_notification(int p_what) {

	switch (p_what) {

		case NOTIFICATION_ENTER_TREE: {

			if (enabled && !get_tree()->is_editor_hint())
				set_fixed_process(true);
			else
				set_fixed_process(false);

			if (get_parent()->cast_to<PhysicsBody2D>()) {
				if (exclude_parent_body)
					exclude.insert(get_parent()->cast_to<PhysicsBody2D>()->get_rid());
				else
					exclude.erase(get_parent()->cast_to<PhysicsBody2D>()->get_rid());
			}
		} break;
		case NOTIFICATION_EXIT_TREE: {

			if (enabled)
				set_fixed_process(false);

		} break;

		case NOTIFICATION_DRAW: {

			if (!get_tree()->is_editor_hint() && !get_tree()->is_debugging_collisions_hint())
				break;
			Transform2D xf;
			xf.rotate(cast_to.angle());
			xf.translate(Vector2(0, cast_to.length()));

			//Vector2 tip = Vector2(0,s->get_length());
			Color dcol = get_tree()->get_debug_collisions_color(); //0.9,0.2,0.2,0.4);
			draw_line(Vector2(), cast_to, dcol, 3);
			Vector<Vector2> pts;
			float tsize = 4;
			pts.push_back(xf.xform(Vector2(0, tsize)));
			pts.push_back(xf.xform(Vector2(0.707 * tsize, 0)));
			pts.push_back(xf.xform(Vector2(-0.707 * tsize, 0)));
			Vector<Color> cols;
			for (int i = 0; i < 3; i++)
				cols.push_back(dcol);

			draw_primitive(pts, cols, Vector<Vector2>()); //small arrow

		} break;

		case NOTIFICATION_FIXED_PROCESS: {

			if (!enabled)
				break;

			_update_raycast_state();

		} break;
	}
}

void RayCast2D::_update_raycast_state() {
	Ref<World2D> w2d = get_world_2d();
	ERR_FAIL_COND(w2d.is_null());

	Physics2DDirectSpaceState *dss = Physics2DServer::get_singleton()->space_get_direct_state(w2d->get_space());
	ERR_FAIL_COND(!dss);

	Transform2D gt = get_global_transform();

	Vector2 to = cast_to;
	if (to == Vector2())
		to = Vector2(0, 0.01);

	Physics2DDirectSpaceState::RayResult rr;

	if (dss->intersect_ray(gt.get_origin(), gt.xform(to), rr, exclude, layer_mask, type_mask)) {

		collided = true;
		against = rr.collider_id;
		collision_point = rr.position;
		collision_normal = rr.normal;
		against_shape = rr.shape;
	} else {
		collided = false;
	}
}

void RayCast2D::force_raycast_update() {
	_update_raycast_state();
}

void RayCast2D::add_exception_rid(const RID &p_rid) {

	exclude.insert(p_rid);
}

void RayCast2D::add_exception(const Object *p_object) {

	ERR_FAIL_NULL(p_object);
	CollisionObject2D *co = ((Object *)p_object)->cast_to<CollisionObject2D>();
	if (!co)
		return;
	add_exception_rid(co->get_rid());
}

void RayCast2D::remove_exception_rid(const RID &p_rid) {

	exclude.erase(p_rid);
}

void RayCast2D::remove_exception(const Object *p_object) {

	ERR_FAIL_NULL(p_object);
	CollisionObject2D *co = ((Object *)p_object)->cast_to<CollisionObject2D>();
	if (!co)
		return;
	remove_exception_rid(co->get_rid());
}

void RayCast2D::clear_exceptions() {

	exclude.clear();
}

void RayCast2D::_bind_methods() {

	ClassDB::bind_method(D_METHOD("set_enabled", "enabled"), &RayCast2D::set_enabled);
	ClassDB::bind_method(D_METHOD("is_enabled"), &RayCast2D::is_enabled);

	ClassDB::bind_method(D_METHOD("set_cast_to", "local_point"), &RayCast2D::set_cast_to);
	ClassDB::bind_method(D_METHOD("get_cast_to"), &RayCast2D::get_cast_to);

	ClassDB::bind_method(D_METHOD("is_colliding"), &RayCast2D::is_colliding);
	ClassDB::bind_method(D_METHOD("force_raycast_update"), &RayCast2D::force_raycast_update);

	ClassDB::bind_method(D_METHOD("get_collider"), &RayCast2D::get_collider);
	ClassDB::bind_method(D_METHOD("get_collider_shape"), &RayCast2D::get_collider_shape);
	ClassDB::bind_method(D_METHOD("get_collision_point"), &RayCast2D::get_collision_point);
	ClassDB::bind_method(D_METHOD("get_collision_normal"), &RayCast2D::get_collision_normal);

	ClassDB::bind_method(D_METHOD("add_exception_rid", "rid"), &RayCast2D::add_exception_rid);
	ClassDB::bind_method(D_METHOD("add_exception", "node"), &RayCast2D::add_exception);

	ClassDB::bind_method(D_METHOD("remove_exception_rid", "rid"), &RayCast2D::remove_exception_rid);
	ClassDB::bind_method(D_METHOD("remove_exception", "node"), &RayCast2D::remove_exception);

	ClassDB::bind_method(D_METHOD("clear_exceptions"), &RayCast2D::clear_exceptions);

	ClassDB::bind_method(D_METHOD("set_layer_mask", "mask"), &RayCast2D::set_layer_mask);
	ClassDB::bind_method(D_METHOD("get_layer_mask"), &RayCast2D::get_layer_mask);

	ClassDB::bind_method(D_METHOD("set_type_mask", "mask"), &RayCast2D::set_type_mask);
	ClassDB::bind_method(D_METHOD("get_type_mask"), &RayCast2D::get_type_mask);

	ClassDB::bind_method(D_METHOD("set_exclude_parent_body", "mask"), &RayCast2D::set_exclude_parent_body);
	ClassDB::bind_method(D_METHOD("get_exclude_parent_body"), &RayCast2D::get_exclude_parent_body);

	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "enabled"), "set_enabled", "is_enabled");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_parent"), "set_exclude_parent_body", "get_exclude_parent_body");
	ADD_PROPERTY(PropertyInfo(Variant::VECTOR2, "cast_to"), "set_cast_to", "get_cast_to");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "layer_mask", PROPERTY_HINT_LAYERS_2D_PHYSICS), "set_layer_mask", "get_layer_mask");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "type_mask", PROPERTY_HINT_FLAGS, "Static,Kinematic,Rigid,Character,Area"), "set_type_mask", "get_type_mask");
}

RayCast2D::RayCast2D() {

	enabled = false;
	against = 0;
	collided = false;
	against_shape = 0;
	layer_mask = 1;
	type_mask = Physics2DDirectSpaceState::TYPE_MASK_COLLISION;
	cast_to = Vector2(0, 50);
	exclude_parent_body = true;
}
