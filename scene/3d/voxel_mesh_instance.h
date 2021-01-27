#ifndef VOXEL_MESH_INSTANCE_H
#define VOXEL_MESH_INSTANCE_H

#include "scene/3d/visual_instance.h"
#include "scene/resources/voxel_mesh.h"

class VoxelMeshInstance : public GeometryInstance {
	GDCLASS(VoxelMeshInstance, GeometryInstance);

protected:
	Ref<VoxelMesh> mesh;

	Vector<Ref<Material> > materials;

	void _mesh_changed();

protected:
	bool _set(const StringName &p_name, const Variant &p_value);
	bool _get(const StringName &p_name, Variant &r_ret) const;
	void _get_property_list(List<PropertyInfo> *p_list) const;

	static void _bind_methods();

public:
    /// MARK: - Init

    VoxelMeshInstance();
	~VoxelMeshInstance();

    /// MARK: - Overrides

    virtual PoolVector<Face3> get_faces(uint32_t p_usage_flags) const;
	virtual AABB get_aabb() const;

    /// MARK: - Public

	void set_mesh(const Ref<VoxelMesh> &p_mesh);
	Ref<VoxelMesh> get_mesh() const;

	int get_surface_material_count() const;
	void set_surface_material(int p_surface, const Ref<Material> &p_material);
	Ref<Material> get_surface_material(int p_surface) const;

};

#endif /* VOXEL_MESH_INSTANCE */
