/*************************************************************************/
/*  circle_shape_2d.cpp                                                  */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                    http://www.godotengine.org                         */
/*************************************************************************/
/* Copyright (c) 2007-2015 Juan Linietsky, Ariel Manzur.                 */
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
#include "circle_shape_2d.h"

#include "servers/physics_2d_server.h"
#include "servers/visual_server.h"
void CircleShape2D::_update_shape() {

	Physics2DServer::get_singleton()->shape_set_data(get_rid(),radius);
	emit_changed();
}


void CircleShape2D::set_radius(real_t p_radius) {

	radius=p_radius;
	_update_shape();
}

real_t CircleShape2D::get_radius() const {

	return radius;
}


void CircleShape2D::_bind_methods() {

	ObjectTypeDB::bind_method(_MD("set_radius","radius"),&CircleShape2D::set_radius);
	ObjectTypeDB::bind_method(_MD("get_radius"),&CircleShape2D::get_radius);

	ADD_PROPERTY( PropertyInfo(Variant::REAL,"radius",PROPERTY_HINT_RANGE,"0.01,16384,0.5"),_SCS("set_radius"),_SCS("get_radius") );

}

Rect2 CircleShape2D::get_rect() const {
	Rect2 rect;
	rect.pos=-Point2(get_radius(),get_radius());
	rect.size=Point2(get_radius(),get_radius())*2.0;
	return rect;
}

void CircleShape2D::draw(const RID& p_to_rid,const Color& p_color) {

	Vector<Vector2> points;
	for(int i=0;i<24;i++) {

		points.push_back(Vector2(Math::cos(i*Math_PI*2/24.0),Math::sin(i*Math_PI*2/24.0))*get_radius());
	}

	Vector<Color> col;
	col.push_back(p_color);
	VisualServer::get_singleton()->canvas_item_add_polygon(p_to_rid,points,col);

}

CircleShape2D::CircleShape2D() : Shape2D( Physics2DServer::get_singleton()->shape_create(Physics2DServer::SHAPE_CIRCLE)) {

	radius=10;
	_update_shape();
}
