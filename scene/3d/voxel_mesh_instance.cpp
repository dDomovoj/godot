#include "voxel_mesh_instance.h"

#include "core/core_string_names.h"
#include "physics_body.h"
#include "scene/resources/material.h"
#include "scene/scene_string_names.h"

/// MARK: - Init

VoxelMeshInstance::VoxelMeshInstance() { }

VoxelMeshInstance::~VoxelMeshInstance() { }

/// MARK: - Private

void VoxelMeshInstance::_mesh_changed() {
	materials.resize(mesh->get_surface_count());
}

/// MARK: - Protected

void VoxelMeshInstance::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_mesh", "mesh"), &VoxelMeshInstance::set_mesh);
	ClassDB::bind_method(D_METHOD("get_mesh"), &VoxelMeshInstance::get_mesh);

	ClassDB::bind_method(D_METHOD("get_surface_material_count"), &VoxelMeshInstance::get_surface_material_count);
	ClassDB::bind_method(D_METHOD("set_surface_material", "surface", "material"), &VoxelMeshInstance::set_surface_material);
	ClassDB::bind_method(D_METHOD("get_surface_material", "surface"), &VoxelMeshInstance::get_surface_material);

	ClassDB::bind_method(D_METHOD("_mesh_changed"), &VoxelMeshInstance::_mesh_changed);

	ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "mesh", PROPERTY_HINT_RESOURCE_TYPE, "VoxelMesh"), "set_mesh", "get_mesh");
}

bool VoxelMeshInstance::_set(const StringName &p_name, const Variant &p_value) {
	//this is not _too_ bad performance wise, really. it only arrives here if the property was not set anywhere else.
	//add to it that it's probably found on first call to _set anyway.

	if (!get_instance().is_valid())
		return false;

	if (p_name.operator String().begins_with("material/")) {
		int idx = p_name.operator String().get_slicec('/', 1).to_int();
		if (idx >= materials.size() || idx < 0)
			return false;

		set_surface_material(idx, p_value);
		return true;
	}

	return false;
}

bool VoxelMeshInstance::_get(const StringName &p_name, Variant &r_ret) const {
	if (!get_instance().is_valid())
		return false;

	if (p_name.operator String().begins_with("material/")) {
		int idx = p_name.operator String().get_slicec('/', 1).to_int();
		if (idx >= materials.size() || idx < 0)
			return false;
		r_ret = materials[idx];
		return true;
	}
	return false;
}

void VoxelMeshInstance::_get_property_list(List<PropertyInfo> *p_list) const {
	if (mesh.is_valid()) {
		for (int i = 0; i < mesh->get_surface_count(); i++) {
			p_list->push_back(PropertyInfo(Variant::OBJECT, "material/" + itos(i), PROPERTY_HINT_RESOURCE_TYPE, "ShaderMaterial,SpatialMaterial"));
		}
	}
}

/// MARK: - Overrides

AABB VoxelMeshInstance::get_aabb() const {
	if (!mesh.is_null())
		return mesh->get_aabb();

	return AABB();
}

PoolVector<Face3> VoxelMeshInstance::get_faces(uint32_t p_usage_flags) const {
    return PoolVector<Face3>();
}

/// MARK: - Public

void VoxelMeshInstance::set_mesh(const Ref<VoxelMesh> &p_mesh) {
	if (mesh == p_mesh)
		return;

	if (mesh.is_valid()) {
		mesh->disconnect(CoreStringNames::get_singleton()->changed, this, SceneStringNames::get_singleton()->_mesh_changed);
		materials.clear();
	}

	mesh = p_mesh;

	if (mesh.is_valid()) {
		mesh->connect(CoreStringNames::get_singleton()->changed, this, SceneStringNames::get_singleton()->_mesh_changed);
		materials.resize(mesh->get_surface_count());
		set_base(mesh->get_rid());
	} else {
		set_base(RID());
	}

	update_gizmo();
	_change_notify();
}

Ref<VoxelMesh> VoxelMeshInstance::get_mesh() const {
	return mesh;
}

int VoxelMeshInstance::get_surface_material_count() const {
	return materials.size();
}

void VoxelMeshInstance::set_surface_material(int p_surface, const Ref<Material> &p_material) {
	ERR_FAIL_INDEX(p_surface, materials.size());

	materials.write[p_surface] = p_material;

	if (materials[p_surface].is_valid())
		VS::get_singleton()->instance_set_surface_material(get_instance(), p_surface, materials[p_surface]->get_rid());
	else
		VS::get_singleton()->instance_set_surface_material(get_instance(), p_surface, RID());
}

Ref<Material> VoxelMeshInstance::get_surface_material(int p_surface) const {
	ERR_FAIL_INDEX_V(p_surface, materials.size(), Ref<Material>());
	return materials[p_surface];
}
