/*************************************************************************/
/*  basis.cpp                                                            */
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
#include <godot/basis.h>

#include "core/math/matrix3.h"
#include "core/variant.h"

#ifdef __cplusplus
extern "C" {
#endif

void _basis_api_anchor() {}

void GDAPI godot_basis_new_with_rows(godot_basis *r_dest, const godot_vector3 *p_x_axis, const godot_vector3 *p_y_axis, const godot_vector3 *p_z_axis) {
	const Vector3 *x_axis = (const Vector3 *)p_x_axis;
	const Vector3 *y_axis = (const Vector3 *)p_y_axis;
	const Vector3 *z_axis = (const Vector3 *)p_z_axis;
	Basis *dest = (Basis *)r_dest;
	*dest = Basis(*x_axis, *y_axis, *z_axis);
}

void GDAPI godot_basis_new_with_axis_and_angle(godot_basis *r_dest, const godot_vector3 *p_axis, const godot_real p_phi) {
	const Vector3 *axis = (const Vector3 *)p_axis;
	Basis *dest = (Basis *)r_dest;
	*dest = Basis(*axis, p_phi);
}

void GDAPI godot_basis_new_with_euler(godot_basis *r_dest, const godot_vector3 *p_euler) {
	const Vector3 *euler = (const Vector3 *)p_euler;
	Basis *dest = (Basis *)r_dest;
	*dest = Basis(*euler);
}

godot_string GDAPI godot_basis_as_string(const godot_basis *p_self) {
	godot_string ret;
	const Basis *self = (const Basis *)p_self;
	memnew_placement(&ret, String(*self));
	return ret;
}

godot_basis GDAPI godot_basis_inverse(const godot_basis *p_self) {
	godot_basis dest;
	const Basis *self = (const Basis *)p_self;
	*((Basis *)&dest) = self->inverse();
	return dest;
}

godot_basis GDAPI godot_basis_transposed(const godot_basis *p_self) {
	godot_basis dest;
	const Basis *self = (const Basis *)p_self;
	*((Basis *)&dest) = self->transposed();
	return dest;
}

godot_basis GDAPI godot_basis_orthonormalized(const godot_basis *p_self) {
	godot_basis dest;
	const Basis *self = (const Basis *)p_self;
	*((Basis *)&dest) = self->orthonormalized();
	return dest;
}

godot_real GDAPI godot_basis_determinant(const godot_basis *p_self) {
	const Basis *self = (const Basis *)p_self;
	return self->determinant();
}

godot_basis GDAPI godot_basis_rotated(const godot_basis *p_self, const godot_vector3 *p_axis, const godot_real p_phi) {
	godot_basis dest;
	const Basis *self = (const Basis *)p_self;
	const Vector3 *axis = (const Vector3 *)p_axis;
	*((Basis *)&dest) = self->rotated(*axis, p_phi);
	return dest;
}

godot_basis GDAPI godot_basis_scaled(const godot_basis *p_self, const godot_vector3 *p_scale) {
	godot_basis dest;
	const Basis *self = (const Basis *)p_self;
	const Vector3 *scale = (const Vector3 *)p_scale;
	*((Basis *)&dest) = self->scaled(*scale);
	return dest;
}

void GDAPI godot_basis_set_scale(godot_basis *p_self, const godot_vector3 *p_scale) {
	Basis *self = (Basis *)p_self;
	const Vector3 *scale = (const Vector3 *)p_scale;
	self->set_scale(*scale);
}

void GDAPI godot_basis_set_rotation_euler(godot_basis *p_self, const godot_vector3 *p_euler) {
	Basis *self = (Basis *)p_self;
	const Vector3 *euler = (const Vector3 *)p_euler;
	self->set_rotation_euler(*euler);
}

void GDAPI godot_basis_set_rotation_axis_angle(godot_basis *p_self, const godot_vector3 *p_axis, const godot_real p_angle) {
	Basis *self = (Basis *)p_self;
	const Vector3 *axis = (const Vector3 *)p_axis;
	self->set_rotation_axis_angle(*axis, p_angle);
}

godot_vector3 GDAPI godot_basis_get_scale(const godot_basis *p_self) {
	godot_vector3 dest;
	const Basis *self = (const Basis *)p_self;
	*((Vector3 *)&dest) = self->get_scale();
	return dest;
}

godot_vector3 GDAPI godot_basis_get_euler(const godot_basis *p_self) {
	godot_vector3 dest;
	const Basis *self = (const Basis *)p_self;
	*((Vector3 *)&dest) = self->get_euler();
	return dest;
}

godot_real GDAPI godot_basis_tdotx(const godot_basis *p_self, const godot_vector3 *p_with) {
	const Basis *self = (const Basis *)p_self;
	const Vector3 *with = (const Vector3 *)p_with;
	return self->tdotx(*with);
}

godot_real GDAPI godot_basis_tdoty(const godot_basis *p_self, const godot_vector3 *p_with) {
	const Basis *self = (const Basis *)p_self;
	const Vector3 *with = (const Vector3 *)p_with;
	return self->tdoty(*with);
}

godot_real GDAPI godot_basis_tdotz(const godot_basis *p_self, const godot_vector3 *p_with) {
	const Basis *self = (const Basis *)p_self;
	const Vector3 *with = (const Vector3 *)p_with;
	return self->tdotz(*with);
}

godot_vector3 GDAPI godot_basis_xform(const godot_basis *p_self, const godot_vector3 *p_v) {
	godot_vector3 dest;
	const Basis *self = (const Basis *)p_self;
	const Vector3 *v = (const Vector3 *)p_v;
	*((Vector3 *)&dest) = self->xform(*v);
	return dest;
}

godot_vector3 GDAPI godot_basis_xform_inv(const godot_basis *p_self, const godot_vector3 *p_v) {
	godot_vector3 dest;
	const Basis *self = (const Basis *)p_self;
	const Vector3 *v = (const Vector3 *)p_v;
	*((Vector3 *)&dest) = self->xform_inv(*v);
	return dest;
}

godot_int GDAPI godot_basis_get_orthogonal_index(const godot_basis *p_self) {
	const Basis *self = (const Basis *)p_self;
	return self->get_orthogonal_index();
}

void GDAPI godot_basis_new(godot_basis *r_dest) {
	Basis *dest = (Basis *)r_dest;
	*dest = Basis();
}

void GDAPI godot_basis_new_with_euler_quat(godot_basis *r_dest, const godot_quat *p_euler) {
	Basis *dest = (Basis *)r_dest;
	const Quat *euler = (const Quat *)p_euler;
	*dest = Basis(*euler);
}

// p_elements is a pointer to an array of 3 (!!) vector3
void GDAPI godot_basis_get_elements(godot_basis *p_self, godot_vector3 *p_elements) {
	const Basis *self = (const Basis *)p_self;
	Vector3 *elements = (Vector3 *)p_elements;
	elements[0] = self->elements[0];
	elements[1] = self->elements[1];
	elements[2] = self->elements[2];
}

godot_vector3 GDAPI godot_basis_get_axis(const godot_basis *p_self, const godot_int p_axis) {
	godot_vector3 dest;
	Vector3 *d = (Vector3 *)&dest;
	const Basis *self = (const Basis *)p_self;
	*d = self->get_axis(p_axis);
	return dest;
}

void GDAPI godot_basis_set_axis(godot_basis *p_self, const godot_int p_axis, const godot_vector3 *p_value) {
	Basis *self = (Basis *)p_self;
	const Vector3 *value = (const Vector3 *)p_value;
	self->set_axis(p_axis, *value);
}

godot_vector3 GDAPI godot_basis_get_row(const godot_basis *p_self, const godot_int p_row) {
	godot_vector3 dest;
	Vector3 *d = (Vector3 *)&dest;
	const Basis *self = (const Basis *)p_self;
	*d = self->get_row(p_row);
	return dest;
}

void GDAPI godot_basis_set_row(godot_basis *p_self, const godot_int p_row, const godot_vector3 *p_value) {
	Basis *self = (Basis *)p_self;
	const Vector3 *value = (const Vector3 *)p_value;
	self->set_row(p_row, *value);
}

godot_bool GDAPI godot_basis_operator_equal(const godot_basis *p_self, const godot_basis *p_b) {
	const Basis *self = (const Basis *)p_self;
	const Basis *b = (const Basis *)p_b;
	return *self == *b;
}

godot_basis GDAPI godot_basis_operator_add(const godot_basis *p_self, const godot_basis *p_b) {
	godot_basis raw_dest;
	Basis *dest = (Basis *)&raw_dest;
	const Basis *self = (const Basis *)p_self;
	const Basis *b = (const Basis *)p_b;
	*dest = *self + *b;
	return raw_dest;
}

godot_basis GDAPI godot_basis_operator_substract(const godot_basis *p_self, const godot_basis *p_b) {
	godot_basis raw_dest;
	Basis *dest = (Basis *)&raw_dest;
	const Basis *self = (const Basis *)p_self;
	const Basis *b = (const Basis *)p_b;
	*dest = *self - *b;
	return raw_dest;
}

godot_basis GDAPI godot_basis_operator_multiply_vector(const godot_basis *p_self, const godot_basis *p_b) {
	godot_basis raw_dest;
	Basis *dest = (Basis *)&raw_dest;
	const Basis *self = (const Basis *)p_self;
	const Basis *b = (const Basis *)p_b;
	*dest = *self * *b;
	return raw_dest;
}

godot_basis GDAPI godot_basis_operator_multiply_scalar(const godot_basis *p_self, const godot_real p_b) {
	godot_basis raw_dest;
	Basis *dest = (Basis *)&raw_dest;
	const Basis *self = (const Basis *)p_self;
	*dest = *self * p_b;
	return raw_dest;
}

#ifdef __cplusplus
}
#endif
