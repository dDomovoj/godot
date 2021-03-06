#ifndef VOXEL_MESH_H
#define VOXEL_MESH_H

#include "core/resource.h"
#include "scene/resources/material.h"
#include "servers/visual_server.h"

class VoxelMesh : public Resource {

	GDCLASS(VoxelMesh, Resource);

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

	float get_voxel_size() const;
	void set_voxel_size(const float p_size);

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

	// Error lightmap_unwrap(const Transform &p_base_transform = Transform(), float p_texel_size = 0.05);

};

VARIANT_ENUM_CAST(VoxelMesh::VoxelArrayType);
VARIANT_ENUM_CAST(VoxelMesh::VoxelPrimitiveType);

#endif /* VOXEL_MESH_H */
