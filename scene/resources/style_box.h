/*************************************************************************/
/*  style_box.h                                                          */
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
#ifndef STYLE_BOX_H
#define STYLE_BOX_H

#include "resource.h"
#include "scene/resources/texture.h"
#include "servers/visual_server.h"
/**
	@author Juan Linietsky <reduzio@gmail.com>
*/
class StyleBox : public Resource {

	GDCLASS(StyleBox, Resource);
	RES_BASE_EXTENSION("stylebox");
	OBJ_SAVE_TYPE(StyleBox);
	float margin[4];

protected:
	virtual float get_style_margin(Margin p_margin) const = 0;
	static void _bind_methods();

public:
	virtual bool test_mask(const Point2 &p_point, const Rect2 &p_rect) const;

	void set_default_margin(Margin p_margin, float p_value);
	float get_default_margin(Margin p_margin) const;
	float get_margin(Margin p_margin) const;
	virtual Size2 get_center_size() const;

	virtual void draw(RID p_canvas_item, const Rect2 &p_rect) const = 0;

	Size2 get_minimum_size() const;
	Point2 get_offset() const;

	StyleBox();
};

class StyleBoxEmpty : public StyleBox {

	GDCLASS(StyleBoxEmpty, StyleBox);
	virtual float get_style_margin(Margin p_margin) const { return 0; }

public:
	virtual void draw(RID p_canvas_item, const Rect2 &p_rect) const {}
	StyleBoxEmpty() {}
};

class StyleBoxTexture : public StyleBox {

	GDCLASS(StyleBoxTexture, StyleBox);

	float expand_margin[4];
	float margin[4];
	Rect2 region_rect;
	Ref<Texture> texture;
	Ref<Texture> normal_map;
	bool draw_center;
	Color modulate;

protected:
	virtual float get_style_margin(Margin p_margin) const;
	static void _bind_methods();

public:
	void set_expand_margin_size(Margin p_expand_margin, float p_size);
	float get_expand_margin_size(Margin p_expand_margin) const;

	void set_margin_size(Margin p_margin, float p_size);
	float get_margin_size(Margin p_margin) const;

	void set_region_rect(const Rect2 &p_region_rect);
	Rect2 get_region_rect() const;

	void set_texture(RES p_texture);
	RES get_texture() const;

	void set_normal_map(RES p_normal_map);
	RES get_normal_map() const;

	void set_draw_center(bool p_draw);
	bool get_draw_center() const;
	virtual Size2 get_center_size() const;

	void set_modulate(const Color &p_modulate);
	Color get_modulate() const;

	virtual void draw(RID p_canvas_item, const Rect2 &p_rect) const;

	StyleBoxTexture();
	~StyleBoxTexture();
};

class StyleBoxFlat : public StyleBox {

	GDCLASS(StyleBoxFlat, StyleBox);

	Color bg_color;
	PoolVector<Color> border_color;

	int border_width[4];
	int expand_margin[4];
	int corner_radius[4];

	bool filled;
	bool blend_border;
	int corner_detail;

protected:
	virtual float get_style_margin(Margin p_margin) const;
	static void _bind_methods();

public:
	//Color
	void set_bg_color(const Color &p_color);
	Color get_bg_color() const;

	void set_light_color(const Color &p_color);
	Color get_light_color() const;

	void set_dark_color(const Color &p_color);
	Color get_dark_color() const;

	//Border Color
	void set_border_color_all(const Color &p_color);
	void set_border_color(Margin p_border, const Color &p_color);
	Color get_border_color(Margin p_border) const;

	//BORDER
	//width
	void set_border_width_all(int p_size);
	int get_border_width_min() const;

	void set_border_width(Margin p_margin, int p_size);
	int get_border_width(Margin p_margin) const;

	//blend
	void set_border_blend(bool p_blend);
	bool get_border_blend() const;

	//CORNER_RADIUS
	void set_corner_radius_all(int radius);
	void set_corner_radius_individual(const int radius_top_left, const int radius_top_right, const int radius_botton_right, const int radius_bottom_left);
	int get_corner_radius_min() const;

	//EXPANDS
	void set_expand_margin_size(Margin p_expand_margin, float p_size);
	float get_expand_margin_size(Margin p_expand_margin) const;

	//FILLED
	void set_filled(bool p_draw);
	bool is_filled() const;

	virtual Size2 get_center_size() const;

	virtual void draw(RID p_canvas_item, const Rect2 &p_rect) const;

	StyleBoxFlat();
	~StyleBoxFlat();
};

// just used to draw lines.
class StyleBoxLine : public StyleBox {

	GDCLASS(StyleBoxLine, StyleBox);
	Color color;
	int thickness;
	bool vertical;
	float grow;

protected:
	virtual float get_style_margin(Margin p_margin) const;
	static void _bind_methods();

public:
	void set_color(const Color &p_color);
	Color get_color() const;

	void set_thickness(int p_thickness);
	int get_thickness() const;

	void set_vertical(bool p_vertical);
	bool is_vertical() const;

	void set_grow(float p_grow);
	float get_grow() const;

	virtual Size2 get_center_size() const;

	virtual void draw(RID p_canvas_item, const Rect2 &p_rect) const;

	StyleBoxLine();
	~StyleBoxLine();
};

#endif
