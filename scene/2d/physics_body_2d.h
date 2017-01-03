/*************************************************************************/
/*  physics_body_2d.h                                                    */
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
#ifndef PHYSICS_BODY_2D_H
#define PHYSICS_BODY_2D_H

#include "scene/2d/collision_object_2d.h"
#include "servers/physics_2d_server.h"
#include "vset.h"


class PhysicsBody2D : public CollisionObject2D {

	GDCLASS(PhysicsBody2D,CollisionObject2D);

	uint32_t mask;
	uint32_t collision_mask;
	Vector2 one_way_collision_direction;
	float one_way_collision_max_depth;


	void _set_layers(uint32_t p_mask);
	uint32_t _get_layers() const;

protected:

	void _notification(int p_what);
	PhysicsBody2D(Physics2DServer::BodyMode p_mode);

	static void _bind_methods();
public:

	void set_layer_mask(uint32_t p_mask);
	uint32_t get_layer_mask() const;

	void set_collision_mask(uint32_t p_mask);
	uint32_t get_collision_mask() const;


	void set_collision_mask_bit(int p_bit, bool p_value);
	bool get_collision_mask_bit(int p_bit) const;

	void set_layer_mask_bit(int p_bit, bool p_value);
	bool get_layer_mask_bit(int p_bit) const;

	void add_collision_exception_with(Node* p_node); //must be physicsbody
	void remove_collision_exception_with(Node* p_node);

	void set_one_way_collision_direction(const Vector2& p_dir);
	Vector2 get_one_way_collision_direction() const;

	void set_one_way_collision_max_depth(float p_dir);
	float get_one_way_collision_max_depth() const;

	PhysicsBody2D();

};

class StaticBody2D : public PhysicsBody2D {

	GDCLASS(StaticBody2D,PhysicsBody2D);

	Vector2 constant_linear_velocity;
	real_t constant_angular_velocity;

	real_t bounce;
	real_t friction;


protected:

	static void _bind_methods();

public:

	void set_friction(real_t p_friction);
	real_t get_friction() const;

	void set_bounce(real_t p_bounce);
	real_t get_bounce() const;


	void set_constant_linear_velocity(const Vector2& p_vel);
	void set_constant_angular_velocity(real_t p_vel);

	Vector2 get_constant_linear_velocity() const;
	real_t get_constant_angular_velocity() const;

	StaticBody2D();
	~StaticBody2D();

};

class RigidBody2D : public PhysicsBody2D {

	GDCLASS(RigidBody2D,PhysicsBody2D);
public:

	enum Mode {
		MODE_RIGID,
		MODE_STATIC,
		MODE_CHARACTER,
		MODE_KINEMATIC,
	};

	enum CCDMode {
		CCD_MODE_DISABLED,
		CCD_MODE_CAST_RAY,
		CCD_MODE_CAST_SHAPE,
	};

private:

	bool can_sleep;
	Physics2DDirectBodyState *state;
	Mode mode;

	real_t bounce;
	real_t mass;
	real_t friction;
	real_t gravity_scale;
	real_t linear_damp;
	real_t angular_damp;

	Vector2 linear_velocity;
	real_t angular_velocity;
	bool sleeping;


	int max_contacts_reported;

	bool custom_integrator;

	CCDMode ccd_mode;


	struct ShapePair {

		int body_shape;
		int local_shape;
		bool tagged;
		bool operator<(const ShapePair& p_sp) const {
			if (body_shape==p_sp.body_shape)
				return local_shape < p_sp.local_shape;
			else
				return body_shape < p_sp.body_shape;
		}

		ShapePair() {}
		ShapePair(int p_bs, int p_ls) { body_shape=p_bs; local_shape=p_ls; }
	};
	struct RigidBody2D_RemoveAction {


		ObjectID body_id;
		ShapePair pair;

	};
	struct BodyState {

		//int rc;
		bool in_scene;
		VSet<ShapePair> shapes;
	};

	struct ContactMonitor {

		bool locked;
		Map<ObjectID,BodyState> body_map;

	};


	ContactMonitor *contact_monitor;
	void _body_enter_tree(ObjectID p_id);
	void _body_exit_tree(ObjectID p_id);


	void _body_inout(int p_status, ObjectID p_instance, int p_body_shape,int p_local_shape);
	void _direct_state_changed(Object *p_state);

	bool _test_motion(const Vector2& p_motion,float p_margin=0.08,const Ref<Physics2DTestMotionResult>& p_result=Ref<Physics2DTestMotionResult>());

protected:

	static void _bind_methods();
public:

	void set_mode(Mode p_mode);
	Mode get_mode() const;

	void set_mass(real_t p_mass);
	real_t get_mass() const;

	void set_inertia(real_t p_inertia);
	real_t get_inertia() const;

	void set_weight(real_t p_weight);
	real_t get_weight() const;

	void set_friction(real_t p_friction);
	real_t get_friction() const;

	void set_bounce(real_t p_bounce);
	real_t get_bounce() const;

	void set_gravity_scale(real_t p_gravity_scale);
	real_t get_gravity_scale() const;

	void set_linear_damp(real_t p_linear_damp);
	real_t get_linear_damp() const;

	void set_angular_damp(real_t p_angular_damp);
	real_t get_angular_damp() const;

	void set_linear_velocity(const Vector2& p_velocity);
	Vector2 get_linear_velocity() const;

	void set_axis_velocity(const Vector2& p_axis);

	void set_angular_velocity(real_t p_velocity);
	real_t get_angular_velocity() const;

	void set_use_custom_integrator(bool p_enable);
	bool is_using_custom_integrator();

	void set_sleeping(bool p_sleeping);
	bool is_sleeping() const;

	void set_can_sleep(bool p_active);
	bool is_able_to_sleep() const;

	void set_contact_monitor(bool p_enabled);
	bool is_contact_monitor_enabled() const;

	void set_max_contacts_reported(int p_amount);
	int get_max_contacts_reported() const;

	void set_continuous_collision_detection_mode(CCDMode p_mode);
	CCDMode get_continuous_collision_detection_mode() const;

	void apply_impulse(const Vector2& p_offset, const Vector2& p_impulse);

	void set_applied_force(const Vector2& p_force);
	Vector2 get_applied_force() const;

	void set_applied_torque(const float p_torque);
	float get_applied_torque() const;

	void add_force(const Vector2& p_offset, const Vector2& p_force);



	Array get_colliding_bodies() const; //function for script

	RigidBody2D();
	~RigidBody2D();

};

VARIANT_ENUM_CAST(RigidBody2D::Mode);
VARIANT_ENUM_CAST(RigidBody2D::CCDMode);



class KinematicBody2D : public PhysicsBody2D {

	GDCLASS(KinematicBody2D,PhysicsBody2D);

	float margin;
	bool colliding;
	Vector2 collision;
	Vector2 normal;
	Vector2 collider_vel;
	ObjectID collider;
	int collider_shape;
	Variant collider_metadata;
	Vector2 travel;

	Vector2 move_and_slide_floor_velocity;
	bool move_and_slide_on_floor;
	bool move_and_slide_on_ceiling;
	bool move_and_slide_on_wall;
	Array move_and_slide_colliders;

	Variant _get_collider() const;

	_FORCE_INLINE_ bool _ignores_mode(Physics2DServer::BodyMode) const;
protected:

	static void _bind_methods();
public:

	Vector2 move(const Vector2& p_motion);
	Vector2 move_to(const Vector2& p_position);

	bool test_move(const Matrix32 &p_from, const Vector2& p_motion);
	bool is_colliding() const;

	Vector2 get_travel() const;
	void revert_motion();

	Vector2 get_collision_pos() const;
	Vector2 get_collision_normal() const;
	Vector2 get_collider_velocity() const;
	ObjectID get_collider() const;
	int get_collider_shape() const;
	Variant get_collider_metadata() const;

	void set_collision_margin(float p_margin);
	float get_collision_margin() const;

	Vector2 move_and_slide(const Vector2& p_linear_velocity, const Vector2& p_floor_direction=Vector2(0,0), float p_slope_stop_min_velocity=5, int p_max_bounces=4);
	bool is_move_and_slide_on_floor() const;
	bool is_move_and_slide_on_wall() const;
	bool is_move_and_slide_on_ceiling() const;
	Array get_move_and_slide_colliders() const;


	KinematicBody2D();
	~KinematicBody2D();

};


#endif // PHYSICS_BODY_2D_H
