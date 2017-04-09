/*************************************************************************/
/*  godot_basis.h                                                        */
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
#ifndef GODOT_BASIS_H
#define GODOT_BASIS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#ifndef GODOT_CORE_API_GODOT_BASIS_TYPE_DEFINED
typedef struct godot_basis {
	uint8_t _dont_touch_that[36];
} godot_basis;
#endif

#include "../godot.h"

void GDAPI godot_basis_new(godot_basis *p_basis);
void GDAPI godot_basis_new_with_euler_quat(godot_basis *p_basis, const godot_quat *p_euler);
void GDAPI godot_basis_new_with_euler(godot_basis *p_basis, const godot_vector3 *p_euler);

godot_quat GDAPI godot_basis_as_quat(const godot_basis *p_basis);
godot_vector3 GDAPI godot_basis_get_euler(const godot_basis *p_basis);

/*
 * p_elements is a pointer to an array of 3 (!!) vector3
 */
void GDAPI godot_basis_get_elements(godot_basis *p_basis, godot_vector3 *p_elements);

#ifdef __cplusplus
}
#endif

#endif // GODOT_BASIS_H
