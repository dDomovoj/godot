#ifndef VOXEL_MESH_H
#define VOXEL_MESH_H

// #include "core/math/face3.h"
// #include "core/math/triangle_mesh.h"
#include "core/resource.h"
#include "scene/resources/material.h"
// #include "scene/resources/shape.h"
#include "servers/visual_server.h"
// #include "scene/resources/mesh.h"

class VoxelMesh : public Resource {

	GDCLASS(VoxelMesh, Resource);
	// RES_BASE_EXTENSION("vxlmesh");

private:
	struct VoxelSurface {
		String name;
		AABB aabb;
		Ref<Material> material;
	};
	Vector<VoxelSurface> surfaces;
	RID mesh;
	AABB aabb;

	void _recompute_aabb();

protected:

	bool _get(const StringName &p_name, Variant &r_ret) const;
	bool _set(const StringName &p_name, const Variant &p_value);
	
	void _get_property_list(List<PropertyInfo> *p_list) const;

	static void _bind_methods();

public:

/// MARK: - Types
enum VoxelArrayType {
	VOXEL_ARRAY_VERTEX = VisualServer::VOXEL_ARRAY_VERTEX,
	VOXEL_ARRAY_NORMAL = VisualServer::VOXEL_ARRAY_NORMAL,
	VOXEL_ARRAY_TEX_UV = VisualServer::VOXEL_ARRAY_TEX_UV,
	VOXEL_ARRAY_INDEX = VisualServer::VOXEL_ARRAY_INDEX,
	VOXEL_ARRAY_MAX = VisualServer::VOXEL_ARRAY_MAX
};

// enum VoxelArrayFormat {
// 	VOXEL_ARRAY_FORMAT_VERTEX = 1 << VOXEL_ARRAY_VERTEX,
// 	VOXEL_ARRAY_FORMAT_NORMAL = 1 << VOXEL_ARRAY_NORMAL,
// 	VOXEL_ARRAY_FORMAT_TEX_UV = 1 << VOXEL_ARRAY_TEX_UV,
// 	VOXEL_ARRAY_FORMAT_INDEX = 1 << VOXEL_ARRAY_INDEX,
// };

enum VoxelPrimitiveType {
	VOXEL_PRIMITIVE_TRIANGLES = VisualServer::VOXEL_PRIMITIVE_TRIANGLES,
};

/// MARK: - Init

	VoxelMesh();
	~VoxelMesh();

/// MARK: - Overrides

	virtual RID get_rid() const;
	virtual void reload_from_file();

/// MARK: - Public

	AABB get_aabb() const;

	int get_surface_count() const;
	Array surface_get_arrays(int p_surface) const;
	void add_surface_from_arrays(const Array &p_arrays, VoxelPrimitiveType p_primitive, const int p_uv_size);
	void surface_remove(int p_idx);

	void surface_update_region(int p_surface, int p_offset, const PoolVector<uint8_t> &p_data);

	virtual Ref<Material> surface_get_material(int p_idx) const;
	virtual void surface_set_material(int p_idx, const Ref<Material> &p_material);

	String surface_get_name(int p_idx) const;
	void surface_set_name(int p_idx, const String &p_name);
	int surface_find_by_name(const String &p_name) const;

	// void add_surface(VoxelPrimitiveType p_primitive, const PoolVector<uint8_t> &p_array, int p_vertex_count, const PoolVector<uint8_t> &p_index_array, int p_index_count, const AABB &p_aabb);

	// void surface_set_custom_aabb(int p_idx, const AABB &p_aabb); //only recognized by driver

	// int surface_get_array_len(int p_idx) const;
	// int surface_get_array_index_len(int p_idx) const;
	// uint32_t surface_get_format(int p_idx) const;
	// VoxelPrimitiveType surface_get_primitive_type(int p_idx) const;
	// bool surface_is_alpha_sorting_enabled(int p_idx) const;

	// void add_surface_from_mesh_data(const Geometry::MeshData &p_mesh_data);

	// void set_custom_aabb(const AABB &p_custom);
	// AABB get_custom_aabb() const;

	// void regen_normalmaps();

	// Error lightmap_unwrap(const Transform &p_base_transform = Transform(), float p_texel_size = 0.05);

};

VARIANT_ENUM_CAST(VoxelMesh::VoxelArrayType);
// VARIANT_ENUM_CAST(VoxelMesh::VoxelArrayFormat);
VARIANT_ENUM_CAST(VoxelMesh::VoxelPrimitiveType);

#endif /* VOXEL_MESH_H */
