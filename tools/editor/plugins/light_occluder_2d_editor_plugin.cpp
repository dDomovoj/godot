/*************************************************************************/
/*  light_occluder_2d_editor_plugin.cpp                                  */
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
#include "light_occluder_2d_editor_plugin.h"

#include "canvas_item_editor_plugin.h"
#include "os/file_access.h"
#include "tools/editor/editor_settings.h"

void LightOccluder2DEditor::_notification(int p_what) {

	switch(p_what) {

		case NOTIFICATION_READY: {

			button_create->set_icon( get_icon("Edit","EditorIcons"));
			button_edit->set_icon( get_icon("MovePoint","EditorIcons"));
			button_edit->set_pressed(true);
			get_tree()->connect("node_removed",this,"_node_removed");
			create_poly->connect("confirmed",this,"_create_poly");

		} break;
		case NOTIFICATION_FIXED_PROCESS: {


		} break;
	}

}
void LightOccluder2DEditor::_node_removed(Node *p_node) {

	if(p_node==node) {
		node=NULL;
		hide();
		canvas_item_editor->get_viewport_control()->update();
	}

}


void LightOccluder2DEditor::_menu_option(int p_option) {

	switch(p_option) {

		case MODE_CREATE: {

			mode=MODE_CREATE;
			button_create->set_pressed(true);
			button_edit->set_pressed(false);
		} break;
		case MODE_EDIT: {

			mode=MODE_EDIT;
			button_create->set_pressed(false);
			button_edit->set_pressed(true);
		} break;

	}
}

void LightOccluder2DEditor::_wip_close(bool p_closed) {

	undo_redo->create_action(TTR("Create Poly"));
	undo_redo->add_undo_method(node->get_occluder_polygon().ptr(),"set_polygon",node->get_occluder_polygon()->get_polygon());
	undo_redo->add_do_method(node->get_occluder_polygon().ptr(),"set_polygon",wip);
	undo_redo->add_undo_method(node->get_occluder_polygon().ptr(),"set_closed",node->get_occluder_polygon()->is_closed());
	undo_redo->add_do_method(node->get_occluder_polygon().ptr(),"set_closed",p_closed);

	undo_redo->add_do_method(canvas_item_editor->get_viewport_control(),"update");
	undo_redo->add_undo_method(canvas_item_editor->get_viewport_control(),"update");
	undo_redo->commit_action();
	wip.clear();
	wip_active=false;
	mode=MODE_EDIT;
	button_edit->set_pressed(true);
	button_create->set_pressed(false);
	edited_point=-1;
}

bool LightOccluder2DEditor::forward_gui_input(const InputEvent& p_event) {


	if (!node)
		return false;

	if (node->get_occluder_polygon().is_null()) {
		if (p_event.type==InputEvent::MOUSE_BUTTON && p_event.mouse_button.button_index==1 && p_event.mouse_button.pressed) {
			create_poly->set_text("No OccluderPolygon2D resource on this node.\nCreate and assign one?");
			create_poly->popup_centered_minsize();
		}
		return (p_event.type==InputEvent::MOUSE_BUTTON && p_event.mouse_button.button_index==1);
	}
	switch(p_event.type) {

		case InputEvent::MOUSE_BUTTON: {

			const InputEventMouseButton &mb=p_event.mouse_button;

			Matrix32 xform = canvas_item_editor->get_canvas_transform() * node->get_global_transform();


			Vector2 gpoint = Point2(mb.x,mb.y);
			Vector2 cpoint = canvas_item_editor->get_canvas_transform().affine_inverse().xform(gpoint);
			cpoint=canvas_item_editor->snap_point(cpoint);
			cpoint = node->get_global_transform().affine_inverse().xform(cpoint);

			Vector<Vector2> poly = Variant(node->get_occluder_polygon()->get_polygon());

			//first check if a point is to be added (segment split)
			real_t grab_treshold=EDITOR_DEF("editors/poly_editor/point_grab_radius",8);

			switch(mode) {


				case MODE_CREATE: {

					if (mb.button_index==BUTTON_LEFT && mb.pressed) {


						if (!wip_active) {

							wip.clear();
							wip.push_back( cpoint );
							wip_active=true;
							edited_point_pos=cpoint;
							canvas_item_editor->get_viewport_control()->update();
							edited_point=1;
							return true;
						} else {


							if (wip.size()>1 && xform.xform(wip[0]).distance_to(gpoint)<grab_treshold) {
								//wip closed
								_wip_close(true);

								return true;
							} else if (wip.size()>1 && xform.xform(wip[wip.size()-1]).distance_to(gpoint)<grab_treshold) {
									//wip closed
									_wip_close(false);
									return true;

							} else {

								wip.push_back( cpoint );
								edited_point=wip.size();
								canvas_item_editor->get_viewport_control()->update();
								return true;

								//add wip point
							}
						}
					} else if (mb.button_index==BUTTON_RIGHT && mb.pressed && wip_active) {
						_wip_close(true);
					}



				} break;

				case MODE_EDIT: {

					if (mb.button_index==BUTTON_LEFT) {
						if (mb.pressed) {

							if (mb.mod.control) {


								if (poly.size() < 3) {

									undo_redo->create_action(TTR("Edit Poly"));
									undo_redo->add_undo_method(node->get_occluder_polygon().ptr(),"set_polygon",poly);
									poly.push_back(cpoint);
									undo_redo->add_do_method(node->get_occluder_polygon().ptr(),"set_polygon",poly);
									undo_redo->add_do_method(canvas_item_editor->get_viewport_control(),"update");
									undo_redo->add_undo_method(canvas_item_editor->get_viewport_control(),"update");
									undo_redo->commit_action();
									return true;
								}

								//search edges
								int closest_idx=-1;
								Vector2 closest_pos;
								real_t closest_dist=1e10;
								for(int i=0;i<poly.size();i++) {

									Vector2 points[2] ={ xform.xform(poly[i]),
										xform.xform(poly[(i+1)%poly.size()]) };

									Vector2 cp = Geometry::get_closest_point_to_segment_2d(gpoint,points);
									if (cp.distance_squared_to(points[0])<CMP_EPSILON2 || cp.distance_squared_to(points[1])<CMP_EPSILON2)
										continue; //not valid to reuse point

									real_t d = cp.distance_to(gpoint);
									if (d<closest_dist && d<grab_treshold) {
										closest_dist=d;
										closest_pos=cp;
										closest_idx=i;
									}


								}

								if (closest_idx>=0) {

									pre_move_edit=poly;
									poly.insert(closest_idx+1,xform.affine_inverse().xform(closest_pos));
									edited_point=closest_idx+1;
									edited_point_pos=xform.affine_inverse().xform(closest_pos);
									node->get_occluder_polygon()->set_polygon(Variant(poly));
									canvas_item_editor->get_viewport_control()->update();
									return true;
								}
							} else {

								//look for points to move

								int closest_idx=-1;
								Vector2 closest_pos;
								real_t closest_dist=1e10;
								for(int i=0;i<poly.size();i++) {

									Vector2 cp =xform.xform(poly[i]);

									real_t d = cp.distance_to(gpoint);
									if (d<closest_dist && d<grab_treshold) {
										closest_dist=d;
										closest_pos=cp;
										closest_idx=i;
									}

								}

								if (closest_idx>=0) {

									pre_move_edit=poly;
									edited_point=closest_idx;
									edited_point_pos=xform.affine_inverse().xform(closest_pos);
									canvas_item_editor->get_viewport_control()->update();
									return true;
								}
							}
						} else {

							if (edited_point!=-1) {

								//apply

								ERR_FAIL_INDEX_V(edited_point,poly.size(),false);
								poly[edited_point]=edited_point_pos;
								undo_redo->create_action(TTR("Edit Poly"));
								undo_redo->add_do_method(node->get_occluder_polygon().ptr(),"set_polygon",poly);
								undo_redo->add_undo_method(node->get_occluder_polygon().ptr(),"set_polygon",pre_move_edit);
								undo_redo->add_do_method(canvas_item_editor->get_viewport_control(),"update");
								undo_redo->add_undo_method(canvas_item_editor->get_viewport_control(),"update");
								undo_redo->commit_action();

								edited_point=-1;
								return true;
							}
						}
					} if (mb.button_index==BUTTON_RIGHT && mb.pressed && edited_point==-1) {



						int closest_idx=-1;
						Vector2 closest_pos;
						real_t closest_dist=1e10;
						for(int i=0;i<poly.size();i++) {

							Vector2 cp =xform.xform(poly[i]);

							real_t d = cp.distance_to(gpoint);
							if (d<closest_dist && d<grab_treshold) {
								closest_dist=d;
								closest_pos=cp;
								closest_idx=i;
							}

						}

						if (closest_idx>=0) {


							undo_redo->create_action(TTR("Edit Poly (Remove Point)"));
							undo_redo->add_undo_method(node->get_occluder_polygon().ptr(),"set_polygon",poly);
							poly.remove(closest_idx);
							undo_redo->add_do_method(node->get_occluder_polygon().ptr(),"set_polygon",poly);
							undo_redo->add_do_method(canvas_item_editor->get_viewport_control(),"update");
							undo_redo->add_undo_method(canvas_item_editor->get_viewport_control(),"update");
							undo_redo->commit_action();
							return true;
						}

					}



				} break;
			}



		} break;
		case InputEvent::MOUSE_MOTION: {

			const InputEventMouseMotion &mm=p_event.mouse_motion;

			if (edited_point!=-1 && (wip_active || mm.button_mask&BUTTON_MASK_LEFT)) {

				Vector2 gpoint = Point2(mm.x,mm.y);
				Vector2 cpoint = canvas_item_editor->get_canvas_transform().affine_inverse().xform(gpoint);
				cpoint=canvas_item_editor->snap_point(cpoint);
				edited_point_pos = node->get_global_transform().affine_inverse().xform(cpoint);

				canvas_item_editor->get_viewport_control()->update();

			}

		} break;
	}

	return false;
}
void LightOccluder2DEditor::_canvas_draw() {

	if (!node || !node->get_occluder_polygon().is_valid())
		return;

	Control *vpc = canvas_item_editor->get_viewport_control();

	Vector<Vector2> poly;

	if (wip_active)
		poly=wip;
	else
		poly=Variant(node->get_occluder_polygon()->get_polygon());


	Matrix32 xform = canvas_item_editor->get_canvas_transform() * node->get_global_transform();
	Ref<Texture> handle= get_icon("EditorHandle","EditorIcons");

	for(int i=0;i<poly.size();i++) {


		Vector2 p,p2;
		p = i==edited_point ? edited_point_pos : poly[i];
		if ((wip_active && i==poly.size()-1) || (((i+1)%poly.size())==edited_point))
			p2=edited_point_pos;
		else
			p2 = poly[(i+1)%poly.size()];

		Vector2 point = xform.xform(p);
		Vector2 next_point = xform.xform(p2);

		Color col=Color(1,0.3,0.1,0.8);

		if (i==poly.size()-1 && (!node->get_occluder_polygon()->is_closed() || wip_active)) {

		} else {
			vpc->draw_line(point,next_point,col,2);
		}
		vpc->draw_texture(handle,point-handle->get_size()*0.5);
	}
}



void LightOccluder2DEditor::edit(Node *p_collision_polygon) {

	if (!canvas_item_editor) {
		canvas_item_editor=CanvasItemEditor::get_singleton();
	}

	if (p_collision_polygon) {

		node=p_collision_polygon->cast_to<LightOccluder2D>();
		if (!canvas_item_editor->get_viewport_control()->is_connected("draw",this,"_canvas_draw"))
			canvas_item_editor->get_viewport_control()->connect("draw",this,"_canvas_draw");
		wip.clear();
		wip_active=false;
		edited_point=-1;
		canvas_item_editor->get_viewport_control()->update();
	} else {
		node=NULL;

		if (canvas_item_editor->get_viewport_control()->is_connected("draw",this,"_canvas_draw"))
			canvas_item_editor->get_viewport_control()->disconnect("draw",this,"_canvas_draw");

	}

}

void LightOccluder2DEditor::_create_poly()  {

	if (!node)
		return;
	undo_redo->create_action(TTR("Create Occluder Polygon"));
	undo_redo->add_do_method(node,"set_occluder_polygon",Ref<OccluderPolygon2D>(memnew( OccluderPolygon2D)));
	undo_redo->add_undo_method(node,"set_occluder_polygon",Variant(REF()));
	undo_redo->commit_action();
}

void LightOccluder2DEditor::_bind_methods() {

	ClassDB::bind_method(_MD("_menu_option"),&LightOccluder2DEditor::_menu_option);
	ClassDB::bind_method(_MD("_canvas_draw"),&LightOccluder2DEditor::_canvas_draw);
	ClassDB::bind_method(_MD("_node_removed"),&LightOccluder2DEditor::_node_removed);
	ClassDB::bind_method(_MD("_create_poly"),&LightOccluder2DEditor::_create_poly);

}


LightOccluder2DEditor::LightOccluder2DEditor(EditorNode *p_editor) {

	node=NULL;
	canvas_item_editor=NULL;
	editor=p_editor;
	undo_redo = editor->get_undo_redo();

	add_child( memnew( VSeparator ));
	button_create = memnew( ToolButton );
	add_child(button_create);
	button_create->connect("pressed",this,"_menu_option",varray(MODE_CREATE));
	button_create->set_toggle_mode(true);
	button_create->set_tooltip(TTR("Create a new polygon from scratch."));

	button_edit = memnew( ToolButton );
	add_child(button_edit);
	button_edit->connect("pressed",this,"_menu_option",varray(MODE_EDIT));
	button_edit->set_toggle_mode(true);
	button_edit->set_tooltip(TTR("Edit existing polygon:")+"\n"+TTR("LMB: Move Point.")+"\n"+TTR("Ctrl+LMB: Split Segment.")+"\n"+TTR("RMB: Erase Point."));

	create_poly = memnew( ConfirmationDialog );
	add_child(create_poly);
	create_poly->get_ok()->set_text(TTR("Create"));


	//add_constant_override("separation",0);

#if 0
	options = memnew( MenuButton );
	add_child(options);
	options->set_area_as_parent_rect();
	options->set_text("Polygon");
	//options->get_popup()->add_item("Parse BBCode",PARSE_BBCODE);
	options->get_popup()->connect("item_pressed", this,"_menu_option");
#endif

	mode = MODE_EDIT;
	wip_active=false;

}


void LightOccluder2DEditorPlugin::edit(Object *p_object) {

	collision_polygon_editor->edit(p_object->cast_to<Node>());
}

bool LightOccluder2DEditorPlugin::handles(Object *p_object) const {

	return p_object->is_class("LightOccluder2D");
}

void LightOccluder2DEditorPlugin::make_visible(bool p_visible) {

	if (p_visible) {
		collision_polygon_editor->show();
	} else {

		collision_polygon_editor->hide();
		collision_polygon_editor->edit(NULL);
	}

}

LightOccluder2DEditorPlugin::LightOccluder2DEditorPlugin(EditorNode *p_node) {

	editor=p_node;
	collision_polygon_editor = memnew( LightOccluder2DEditor(p_node) );
	CanvasItemEditor::get_singleton()->add_control_to_menu_panel(collision_polygon_editor);

	collision_polygon_editor->hide();



}


LightOccluder2DEditorPlugin::~LightOccluder2DEditorPlugin()
{
}

