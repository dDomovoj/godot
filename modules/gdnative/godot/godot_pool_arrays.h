/*************************************************************************/
/*  godot_pool_arrays.h                                                  */
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
#ifndef GODOT_POOL_ARRAYS_H
#define GODOT_POOL_ARRAYS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/////// PoolByteArray

#ifndef GODOT_CORE_API_GODOT_POOL_BYTE_ARRAY_TYPE_DEFINED
typedef struct godot_pool_byte_array {
	uint8_t _dont_touch_that[8];
} godot_pool_byte_array;
#endif

/////// PoolIntArray

#ifndef GODOT_CORE_API_GODOT_POOL_INT_ARRAY_TYPE_DEFINED
typedef struct godot_pool_int_array {
	uint8_t _dont_touch_that[8];
} godot_pool_int_array;
#endif

/////// PoolRealArray

#ifndef GODOT_CORE_API_GODOT_POOL_REAL_ARRAY_TYPE_DEFINED
typedef struct godot_pool_real_array {
	uint8_t _dont_touch_that[8];
} godot_pool_real_array;
#endif

/////// PoolStringArray

#ifndef GODOT_CORE_API_GODOT_POOL_STRING_ARRAY_TYPE_DEFINED
typedef struct godot_pool_string_array {
	uint8_t _dont_touch_that[8];
} godot_pool_string_array;
#endif

/////// PoolVector2Array

#ifndef GODOT_CORE_API_GODOT_POOL_VECTOR2_ARRAY_TYPE_DEFINED
typedef struct godot_pool_vector2_array {
	uint8_t _dont_touch_that[8];
} godot_pool_vector2_array;
#endif

/////// PoolVector3Array

#ifndef GODOT_CORE_API_GODOT_POOL_VECTOR3_ARRAY_TYPE_DEFINED
typedef struct godot_pool_vector3_array {
	uint8_t _dont_touch_that[8];
} godot_pool_vector3_array;
#endif

/////// PoolColorArray

#ifndef GODOT_CORE_API_GODOT_POOL_COLOR_ARRAY_TYPE_DEFINED
typedef struct godot_pool_color_array {
	uint8_t _dont_touch_that[8];
} godot_pool_color_array;
#endif

#include "godot_array.h"
#include "godot_color.h"
#include "godot_vector2.h"
#include "godot_vector3.h"

#include "../godot.h"

// byte

void GDAPI godot_pool_byte_array_new(godot_pool_byte_array *r_dest);
void GDAPI godot_pool_byte_array_new_copy(godot_pool_byte_array *r_dest, const godot_pool_byte_array *p_src);
void GDAPI godot_pool_byte_array_new_with_array(godot_pool_byte_array *r_dest, const godot_array *p_a);

void GDAPI godot_pool_byte_array_append(godot_pool_byte_array *p_self, const uint8_t p_data);

void GDAPI godot_pool_byte_array_append_array(godot_pool_byte_array *p_self, const godot_pool_byte_array *p_array);

godot_error GDAPI godot_pool_byte_array_insert(godot_pool_byte_array *p_self, const godot_int p_idx, const uint8_t p_data);

void GDAPI godot_pool_byte_array_invert(godot_pool_byte_array *p_self);

void GDAPI godot_pool_byte_array_push_back(godot_pool_byte_array *p_self, const uint8_t p_data);

void GDAPI godot_pool_byte_array_remove(godot_pool_byte_array *p_self, const godot_int p_idx);

void GDAPI godot_pool_byte_array_resize(godot_pool_byte_array *p_self, const godot_int p_size);

void GDAPI godot_pool_byte_array_set(godot_pool_byte_array *p_self, const godot_int p_idx, const uint8_t p_data);
uint8_t GDAPI godot_pool_byte_array_get(const godot_pool_byte_array *p_self, const godot_int p_idx);

godot_int GDAPI godot_pool_byte_array_size(const godot_pool_byte_array *p_self);

void GDAPI godot_pool_byte_array_destroy(godot_pool_byte_array *p_self);

// int

void GDAPI godot_pool_int_array_new(godot_pool_int_array *r_dest);
void GDAPI godot_pool_int_array_new_copy(godot_pool_int_array *r_dest, const godot_pool_int_array *p_src);
void GDAPI godot_pool_int_array_new_with_array(godot_pool_int_array *r_dest, const godot_array *p_a);

void GDAPI godot_pool_int_array_append(godot_pool_int_array *p_self, const godot_int p_data);

void GDAPI godot_pool_int_array_append_array(godot_pool_int_array *p_self, const godot_pool_int_array *p_array);

godot_error GDAPI godot_pool_int_array_insert(godot_pool_int_array *p_self, const godot_int p_idx, const godot_int p_data);

void GDAPI godot_pool_int_array_invert(godot_pool_int_array *p_self);

void GDAPI godot_pool_int_array_push_back(godot_pool_int_array *p_self, const godot_int p_data);

void GDAPI godot_pool_int_array_remove(godot_pool_int_array *p_self, const godot_int p_idx);

void GDAPI godot_pool_int_array_resize(godot_pool_int_array *p_self, const godot_int p_size);

void GDAPI godot_pool_int_array_set(godot_pool_int_array *p_self, const godot_int p_idx, const godot_int p_data);
godot_int GDAPI godot_pool_int_array_get(const godot_pool_int_array *p_self, const godot_int p_idx);

godot_int GDAPI godot_pool_int_array_size(const godot_pool_int_array *p_self);

void GDAPI godot_pool_int_array_destroy(godot_pool_int_array *p_self);

// real

void GDAPI godot_pool_real_array_new(godot_pool_real_array *r_dest);
void GDAPI godot_pool_real_array_new_copy(godot_pool_real_array *r_dest, const godot_pool_real_array *p_src);
void GDAPI godot_pool_real_array_new_with_array(godot_pool_real_array *r_dest, const godot_array *p_a);

void GDAPI godot_pool_real_array_append(godot_pool_real_array *p_self, const godot_real p_data);

void GDAPI godot_pool_real_array_append_array(godot_pool_real_array *p_self, const godot_pool_real_array *p_array);

godot_error GDAPI godot_pool_real_array_insert(godot_pool_real_array *p_self, const godot_int p_idx, const godot_real p_data);

void GDAPI godot_pool_real_array_invert(godot_pool_real_array *p_self);

void GDAPI godot_pool_real_array_push_back(godot_pool_real_array *p_self, const godot_real p_data);

void GDAPI godot_pool_real_array_remove(godot_pool_real_array *p_self, const godot_int p_idx);

void GDAPI godot_pool_real_array_resize(godot_pool_real_array *p_self, const godot_int p_size);

void GDAPI godot_pool_real_array_set(godot_pool_real_array *p_self, const godot_int p_idx, const godot_real p_data);
godot_real GDAPI godot_pool_real_array_get(const godot_pool_real_array *p_self, const godot_int p_idx);

godot_int GDAPI godot_pool_real_array_size(const godot_pool_real_array *p_self);

void GDAPI godot_pool_real_array_destroy(godot_pool_real_array *p_self);

// string

void GDAPI godot_pool_string_array_new(godot_pool_string_array *r_dest);
void GDAPI godot_pool_string_array_new_copy(godot_pool_string_array *r_dest, const godot_pool_string_array *p_src);
void GDAPI godot_pool_string_array_new_with_array(godot_pool_string_array *r_dest, const godot_array *p_a);

void GDAPI godot_pool_string_array_append(godot_pool_string_array *p_self, const godot_string *p_data);

void GDAPI godot_pool_string_array_append_array(godot_pool_string_array *p_self, const godot_pool_string_array *p_array);

godot_error GDAPI godot_pool_string_array_insert(godot_pool_string_array *p_self, const godot_int p_idx, const godot_string *p_data);

void GDAPI godot_pool_string_array_invert(godot_pool_string_array *p_self);

void GDAPI godot_pool_string_array_push_back(godot_pool_string_array *p_self, const godot_string *p_data);

void GDAPI godot_pool_string_array_remove(godot_pool_string_array *p_self, const godot_int p_idx);

void GDAPI godot_pool_string_array_resize(godot_pool_string_array *p_self, const godot_int p_size);

void GDAPI godot_pool_string_array_set(godot_pool_string_array *p_self, const godot_int p_idx, const godot_string *p_data);
godot_string GDAPI godot_pool_string_array_get(const godot_pool_string_array *p_self, const godot_int p_idx);

godot_int GDAPI godot_pool_string_array_size(const godot_pool_string_array *p_self);

void GDAPI godot_pool_string_array_destroy(godot_pool_string_array *p_self);

// vector2

void GDAPI godot_pool_vector2_array_new(godot_pool_vector2_array *r_dest);
void GDAPI godot_pool_vector2_array_new_copy(godot_pool_vector2_array *r_dest, const godot_pool_vector2_array *p_src);
void GDAPI godot_pool_vector2_array_new_with_array(godot_pool_vector2_array *r_dest, const godot_array *p_a);

void GDAPI godot_pool_vector2_array_append(godot_pool_vector2_array *p_self, const godot_vector2 *p_data);

void GDAPI godot_pool_vector2_array_append_array(godot_pool_vector2_array *p_self, const godot_pool_vector2_array *p_array);

godot_error GDAPI godot_pool_vector2_array_insert(godot_pool_vector2_array *p_self, const godot_int p_idx, const godot_vector2 *p_data);

void GDAPI godot_pool_vector2_array_invert(godot_pool_vector2_array *p_self);

void GDAPI godot_pool_vector2_array_push_back(godot_pool_vector2_array *p_self, const godot_vector2 *p_data);

void GDAPI godot_pool_vector2_array_remove(godot_pool_vector2_array *p_self, const godot_int p_idx);

void GDAPI godot_pool_vector2_array_resize(godot_pool_vector2_array *p_self, const godot_int p_size);

void GDAPI godot_pool_vector2_array_set(godot_pool_vector2_array *p_self, const godot_int p_idx, const godot_vector2 *p_data);
godot_vector2 GDAPI godot_pool_vector2_array_get(const godot_pool_vector2_array *p_self, const godot_int p_idx);

godot_int GDAPI godot_pool_vector2_array_size(const godot_pool_vector2_array *p_self);

void GDAPI godot_pool_vector2_array_destroy(godot_pool_vector2_array *p_self);

// vector3

void GDAPI godot_pool_vector3_array_new(godot_pool_vector3_array *r_dest);
void GDAPI godot_pool_vector3_array_new_copy(godot_pool_vector3_array *r_dest, const godot_pool_vector3_array *p_src);
void GDAPI godot_pool_vector3_array_new_with_array(godot_pool_vector3_array *r_dest, const godot_array *p_a);

void GDAPI godot_pool_vector3_array_append(godot_pool_vector3_array *p_self, const godot_vector3 *p_data);

void GDAPI godot_pool_vector3_array_append_array(godot_pool_vector3_array *p_self, const godot_pool_vector3_array *p_array);

godot_error GDAPI godot_pool_vector3_array_insert(godot_pool_vector3_array *p_self, const godot_int p_idx, const godot_vector3 *p_data);

void GDAPI godot_pool_vector3_array_invert(godot_pool_vector3_array *p_self);

void GDAPI godot_pool_vector3_array_push_back(godot_pool_vector3_array *p_self, const godot_vector3 *p_data);

void GDAPI godot_pool_vector3_array_remove(godot_pool_vector3_array *p_self, const godot_int p_idx);

void GDAPI godot_pool_vector3_array_resize(godot_pool_vector3_array *p_self, const godot_int p_size);

void GDAPI godot_pool_vector3_array_set(godot_pool_vector3_array *p_self, const godot_int p_idx, const godot_vector3 *p_data);
godot_vector3 GDAPI godot_pool_vector3_array_get(const godot_pool_vector3_array *p_self, const godot_int p_idx);

godot_int GDAPI godot_pool_vector3_array_size(const godot_pool_vector3_array *p_self);

void GDAPI godot_pool_vector3_array_destroy(godot_pool_vector3_array *p_self);

// color

void GDAPI godot_pool_color_array_new(godot_pool_color_array *r_dest);
void GDAPI godot_pool_color_array_new_copy(godot_pool_color_array *r_dest, const godot_pool_color_array *p_src);
void GDAPI godot_pool_color_array_new_with_array(godot_pool_color_array *r_dest, const godot_array *p_a);

void GDAPI godot_pool_color_array_append(godot_pool_color_array *p_self, const godot_color *p_data);

void GDAPI godot_pool_color_array_append_array(godot_pool_color_array *p_self, const godot_pool_color_array *p_array);

godot_error GDAPI godot_pool_color_array_insert(godot_pool_color_array *p_self, const godot_int p_idx, const godot_color *p_data);

void GDAPI godot_pool_color_array_invert(godot_pool_color_array *p_self);

void GDAPI godot_pool_color_array_push_back(godot_pool_color_array *p_self, const godot_color *p_data);

void GDAPI godot_pool_color_array_remove(godot_pool_color_array *p_self, const godot_int p_idx);

void GDAPI godot_pool_color_array_resize(godot_pool_color_array *p_self, const godot_int p_size);

void GDAPI godot_pool_color_array_set(godot_pool_color_array *p_self, const godot_int p_idx, const godot_color *p_data);
godot_color GDAPI godot_pool_color_array_get(const godot_pool_color_array *p_self, const godot_int p_idx);

godot_int GDAPI godot_pool_color_array_size(const godot_pool_color_array *p_self);

void GDAPI godot_pool_color_array_destroy(godot_pool_color_array *p_self);

#ifdef __cplusplus
}
#endif

#endif // GODOT_POOL_ARRAYS_H
