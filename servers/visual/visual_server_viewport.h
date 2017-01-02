#ifndef VISUALSERVERVIEWPORT_H
#define VISUALSERVERVIEWPORT_H

#include "servers/visual_server.h"
#include "rasterizer.h"
#include "self_list.h"

class VisualServerViewport  {
public:

	struct CanvasBase : public RID_Data {


	};



	struct Viewport : public RID_Data {

		RID self;
		RID parent;

		Size2i size;
		RID camera;
		RID scenario;

		VS::ViewportUpdateMode update_mode;
		RID render_target;
		RID render_target_texture;

		int viewport_to_screen;
		Rect2 viewport_to_screen_rect;

		bool hide_scenario;
		bool hide_canvas;
		bool disable_environment;
		bool disable_3d;

		RID shadow_atlas;
		int shadow_atlas_size;


		VS::ViewportClearMode clear_mode;

		bool rendered_in_prev_frame;

		struct CanvasKey {

			int layer;
			RID canvas;
			bool operator<(const CanvasKey& p_canvas) const { if (layer==p_canvas.layer) return canvas < p_canvas.canvas; return layer<p_canvas.layer; }
			CanvasKey() { layer=0; }
			CanvasKey(const RID& p_canvas, int p_layer) { canvas=p_canvas; layer=p_layer; }
		};

		struct CanvasData {

			CanvasBase *canvas;
			Matrix32 transform;
			int layer;
		};

		Matrix32 global_transform;

		Map<RID,CanvasData> canvas_map;

		Viewport() {
			update_mode=VS::VIEWPORT_UPDATE_WHEN_VISIBLE;
			clear_mode=VS::VIEWPORT_CLEAR_ALWAYS;
			rendered_in_prev_frame=false;
			disable_environment=false;
			viewport_to_screen=0;
			shadow_atlas_size=0;
			disable_3d=false;

		}
	};

	mutable RID_Owner<Viewport> viewport_owner;


	struct ViewportSort {
		_FORCE_INLINE_ bool operator()(const Viewport*p_left,const Viewport* p_right) const {

			bool left_to_screen = p_left->viewport_to_screen_rect.size!=Size2();
			bool right_to_screen = p_right->viewport_to_screen_rect.size!=Size2();

			if (left_to_screen==right_to_screen) {

				return p_left->parent==p_right->self;
			} else {
				return right_to_screen;
			}
		}
	};


	Vector<Viewport*> active_viewports;
private:
	Color clear_color;
	void _draw_viewport(Viewport *p_viewport);
public:


	RID viewport_create();

	void viewport_set_size(RID p_viewport,int p_width,int p_height);

	void viewport_attach_to_screen(RID p_viewport,const Rect2& p_rect=Rect2(),int p_screen=0);
	void viewport_detach(RID p_viewport);

	void viewport_set_active(RID p_viewport,bool p_active);
	void viewport_set_parent_viewport(RID p_viewport,RID p_parent_viewport);
	void viewport_set_update_mode(RID p_viewport,VS::ViewportUpdateMode p_mode);
	void viewport_set_vflip(RID p_viewport,bool p_enable);


	void viewport_set_clear_mode(RID p_viewport,VS::ViewportClearMode p_clear_mode);

	RID viewport_get_texture(RID p_viewport) const;

	void viewport_set_hide_scenario(RID p_viewport,bool p_hide);
	void viewport_set_hide_canvas(RID p_viewport,bool p_hide);
	void viewport_set_disable_environment(RID p_viewport,bool p_disable);
	void viewport_set_disable_3d(RID p_viewport,bool p_disable);

	void viewport_attach_camera(RID p_viewport,RID p_camera);
	void viewport_set_scenario(RID p_viewport,RID p_scenario);
	void viewport_attach_canvas(RID p_viewport,RID p_canvas);
	void viewport_remove_canvas(RID p_viewport,RID p_canvas);
	void viewport_set_canvas_transform(RID p_viewport,RID p_canvas,const Matrix32& p_offset);
	void viewport_set_transparent_background(RID p_viewport,bool p_enabled);

	void viewport_set_global_canvas_transform(RID p_viewport,const Matrix32& p_transform);
	void viewport_set_canvas_layer(RID p_viewport,RID p_canvas,int p_layer);

	void viewport_set_shadow_atlas_size(RID p_viewport,int p_size);
	void viewport_set_shadow_atlas_quadrant_subdivision(RID p_viewport,int p_quadrant,int p_subdiv);

	void viewport_set_msaa(RID p_viewport,VS::ViewportMSAA p_msaa);
	void viewport_set_hdr(RID p_viewport,bool p_enabled);

	void draw_viewports();

	bool free(RID p_rid);

	VisualServerViewport();
};

#endif // VISUALSERVERVIEWPORT_H
