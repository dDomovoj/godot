#ifndef CCBA3375_365D_4A90_9BFD_977B89BA3507
#define CCBA3375_365D_4A90_9BFD_977B89BA3507

#include "core/math/face3.h"
#include "core/math/triangle_mesh.h"
#include "core/resource.h"
#include "scene/resources/material.h"
#include "scene/resources/shape.h"
#include "servers/visual_server.h"
#include "scene/resources/mesh.h"

class VoxelMesh : public Mesh {

	GDCLASS(VoxelMesh, Mesh);
	RES_BASE_EXTENSION("mesh");

private:
	struct Surface {
		String name;
		AABB aabb;
		Ref<Material> material;
		bool is_2d;
	};
	Vector<Surface> surfaces;
	RID mesh;
	AABB aabb;
	BlendShapeMode blend_shape_mode;
	Vector<StringName> blend_shapes;
	AABB custom_aabb;

	void _recompute_aabb();

protected:
	virtual bool _is_generated() const { return false; }

	bool _set(const StringName &p_name, const Variant &p_value);
	bool _get(const StringName &p_name, Variant &r_ret) const;
	void _get_property_list(List<PropertyInfo> *p_list) const;

	static void _bind_methods();

public:
	void add_surface_from_arrays(PrimitiveType p_primitive, const Array &p_arrays, const Array &p_blend_shapes = Array(), uint32_t p_flags = ARRAY_COMPRESS_DEFAULT);
	void add_surface(uint32_t p_format, PrimitiveType p_primitive, const PoolVector<uint8_t> &p_array, int p_vertex_count, const PoolVector<uint8_t> &p_index_array, int p_index_count, const AABB &p_aabb, const Vector<PoolVector<uint8_t> > &p_blend_shapes = Vector<PoolVector<uint8_t> >(), const Vector<AABB> &p_bone_aabbs = Vector<AABB>());

	Array surface_get_arrays(int p_surface) const;
	Array surface_get_blend_shape_arrays(int p_surface) const;

	void add_blend_shape(const StringName &p_name);
	int get_blend_shape_count() const;
	StringName get_blend_shape_name(int p_index) const;
	void clear_blend_shapes();

	void set_blend_shape_mode(BlendShapeMode p_mode);
	BlendShapeMode get_blend_shape_mode() const;

	void surface_update_region(int p_surface, int p_offset, const PoolVector<uint8_t> &p_data);

	int get_surface_count() const;
	void surface_remove(int p_idx);

	void surface_set_custom_aabb(int p_idx, const AABB &p_aabb); //only recognized by driver

	int surface_get_array_len(int p_idx) const;
	int surface_get_array_index_len(int p_idx) const;
	uint32_t surface_get_format(int p_idx) const;
	PrimitiveType surface_get_primitive_type(int p_idx) const;
	bool surface_is_alpha_sorting_enabled(int p_idx) const;

	virtual void surface_set_material(int p_idx, const Ref<Material> &p_material);
	virtual Ref<Material> surface_get_material(int p_idx) const;

	int surface_find_by_name(const String &p_name) const;
	void surface_set_name(int p_idx, const String &p_name);
	String surface_get_name(int p_idx) const;

	void add_surface_from_mesh_data(const Geometry::MeshData &p_mesh_data);

	void set_custom_aabb(const AABB &p_custom);
	AABB get_custom_aabb() const;

	AABB get_aabb() const;
	virtual RID get_rid() const;

	void regen_normalmaps();

	Error lightmap_unwrap(const Transform &p_base_transform = Transform(), float p_texel_size = 0.05);

	virtual void reload_from_file();

	VoxelMesh();

	~VoxelMesh();
};
#endif /* CCBA3375_365D_4A90_9BFD_977B89BA3507 */
