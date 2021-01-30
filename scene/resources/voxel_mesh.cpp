#include "voxel_mesh.h"

// #include "core/pair.h"

#include <stdlib.h>

/// MARK: - Init

VoxelMesh::VoxelMesh() {
	mesh = VisualServer::get_singleton()->voxel_mesh_create();
}

VoxelMesh::~VoxelMesh() {
	VisualServer::get_singleton()->free(mesh);
}

/// MARK: - Lifecycle

void VoxelMesh::_bind_methods() {
	// ClassDB::bind_method(D_METHOD("set_lightmap_size_hint", "size"), &Mesh::set_lightmap_size_hint);
	// ClassDB::bind_method(D_METHOD("get_lightmap_size_hint"), &Mesh::get_lightmap_size_hint);
	
	// ADD_PROPERTY(PropertyInfo(Variant::VECTOR2, "lightmap_size_hint"), "set_lightmap_size_hint", "get_lightmap_size_hint");
	ClassDB::bind_method(D_METHOD("get_aabb"), &VoxelMesh::get_aabb);

	ClassDB::bind_method(D_METHOD("add_surface_from_arrays", "arrays", "primitive", "uv_size"), &VoxelMesh::add_surface_from_arrays);
	ClassDB::bind_method(D_METHOD("surface_remove", "surf_idx"), &VoxelMesh::surface_remove);
	ClassDB::bind_method(D_METHOD("get_surface_count"), &VoxelMesh::get_surface_count);
	ClassDB::bind_method(D_METHOD("surface_get_arrays", "surf_idx"), &VoxelMesh::surface_get_arrays);

	ClassDB::bind_method(D_METHOD("surface_set_material", "surf_idx", "material"), &VoxelMesh::surface_set_material);
	ClassDB::bind_method(D_METHOD("surface_get_material", "surf_idx"), &VoxelMesh::surface_get_material);
	
	ClassDB::bind_method(D_METHOD("surface_update_region", "surf_idx", "offset", "data"), &VoxelMesh::surface_update_region);
	// ClassDB::bind_method(D_METHOD("surface_get_array_len", "surf_idx"), &VoxelMesh::surface_get_array_len);
	// ClassDB::bind_method(D_METHOD("surface_get_array_index_len", "surf_idx"), &VoxelMesh::surface_get_array_index_len);
	// ClassDB::bind_method(D_METHOD("surface_get_format", "surf_idx"), &VoxelMesh::surface_get_format);
	// ClassDB::bind_method(D_METHOD("surface_get_primitive_type", "surf_idx"), &VoxelMesh::surface_get_primitive_type);
	ClassDB::bind_method(D_METHOD("surface_find_by_name", "name"), &VoxelMesh::surface_find_by_name);
	ClassDB::bind_method(D_METHOD("surface_set_name", "surf_idx", "name"), &VoxelMesh::surface_set_name);
	ClassDB::bind_method(D_METHOD("surface_get_name", "surf_idx"), &VoxelMesh::surface_get_name);
	// ClassDB::bind_method(D_METHOD("regen_normalmaps"), &VoxelMesh::regen_normalmaps);
	// ClassDB::set_method_flags(get_class_static(), _scs_create("regen_normalmaps"), METHOD_FLAGS_DEFAULT | METHOD_FLAG_EDITOR);
	// ClassDB::bind_method(D_METHOD("lightmap_unwrap", "transform", "texel_size"), &VoxelMesh::lightmap_unwrap);
	// ClassDB::set_method_flags(get_class_static(), _scs_create("lightmap_unwrap"), METHOD_FLAGS_DEFAULT | METHOD_FLAG_EDITOR);

	// ClassDB::bind_method(D_METHOD("set_custom_aabb", "aabb"), &VoxelMesh::set_custom_aabb);
	// ClassDB::bind_method(D_METHOD("get_custom_aabb"), &VoxelMesh::get_custom_aabb);

	BIND_ENUM_CONSTANT(VOXEL_PRIMITIVE_TRIANGLES);

	BIND_ENUM_CONSTANT(VOXEL_ARRAY_VERTEX);
	BIND_ENUM_CONSTANT(VOXEL_ARRAY_NORMAL);
	BIND_ENUM_CONSTANT(VOXEL_ARRAY_TEX_UV);
	BIND_ENUM_CONSTANT(VOXEL_ARRAY_INDEX);
	BIND_ENUM_CONSTANT(VOXEL_ARRAY_MAX);

	// BIND_ENUM_CONSTANT(VOXEL_ARRAY_FORMAT_VERTEX);
	// BIND_ENUM_CONSTANT(VOXEL_ARRAY_FORMAT_NORMAL);
	// BIND_ENUM_CONSTANT(VOXEL_ARRAY_FORMAT_TEX_UV);
	// BIND_ENUM_CONSTANT(VOXEL_ARRAY_FORMAT_INDEX);
}

/// MARK: - Private

void VoxelMesh::_recompute_aabb() {
	aabb = AABB();
	for (int i = 0; i < surfaces.size(); i++) {
		if (i == 0)
			aabb = surfaces[i].aabb;
		else
			aabb.merge_with(surfaces[i].aabb);
	}
}

/// MARK: - Protected

bool VoxelMesh::_get(const StringName &p_name, Variant &r_ret) const {
	String sname = p_name;

	if (sname.begins_with("surface_")) {

		int sl = sname.find("/");
		if (sl == -1)
			return false;
		int idx = sname.substr(8, sl - 8).to_int() - 1;
		String what = sname.get_slicec('/', 1);
		if (what == "material")
			r_ret = surface_get_material(idx);
		else if (what == "name")
			r_ret = surface_get_name(idx);
		return true;
	} else if (!sname.begins_with("surfaces"))
		return false;

	int idx = sname.get_slicec('/', 1).to_int();
	ERR_FAIL_INDEX_V(idx, surfaces.size(), false);

	Dictionary d;

	d["uv_size"] = 1;
	d["array_data"] = VS::get_singleton()->voxel_mesh_surface_get_array(mesh, idx);
	d["vertex_count"] = VS::get_singleton()->voxel_mesh_surface_get_array_len(mesh, idx);
	d["array_index_data"] = VS::get_singleton()->voxel_mesh_surface_get_index_array(mesh, idx);
	d["index_count"] = VS::get_singleton()->voxel_mesh_surface_get_array_index_len(mesh, idx);
	d["primitive"] = VS::get_singleton()->voxel_mesh_surface_get_primitive_type(mesh, idx);
	d["aabb"] = VS::get_singleton()->voxel_mesh_surface_get_aabb(mesh, idx);

	Ref<Material> m = surface_get_material(idx);
	if (m.is_valid())
		d["material"] = m;
	String n = surface_get_name(idx);
	if (n != "")
		d["name"] = n;

	r_ret = d;

	return true; 
}

bool VoxelMesh::_set(const StringName &p_name, const Variant &p_value) {
	String sname = p_name;

	if (sname.begins_with("surface_")) {

		int sl = sname.find("/");
		if (sl == -1)
			return false;
		int idx = sname.substr(8, sl - 8).to_int() - 1;
		String what = sname.get_slicec('/', 1);
		if (what == "material")
			surface_set_material(idx, p_value);
		else if (what == "name")
			surface_set_name(idx, p_value);
		return true;
	}

	if (!sname.begins_with("surfaces"))
		return false;

	int idx = sname.get_slicec('/', 1).to_int();
	String what = sname.get_slicec('/', 2);

	if (idx == surfaces.size()) {

		//create
		Dictionary d = p_value;
		ERR_FAIL_COND_V(!d.has("primitive"), false);

		if (d.has("arrays")) {
			//old format
			ERR_FAIL_COND_V(!d.has("morph_arrays"), false);
			add_surface_from_arrays(d["arrays"], VoxelPrimitiveType(int(d["primitive"])), d["uv_size"]);

		} else if (d.has("array_data")) {

			PoolVector<uint8_t> array_data = d["array_data"];
			PoolVector<uint8_t> array_index_data;
			if (d.has("array_index_data"))
				array_index_data = d["array_index_data"];

			uint32_t primitive = d["primitive"];

			ERR_FAIL_COND_V(!d.has("vertex_count"), false);
			int vertex_count = d["vertex_count"];

			int index_count = 0;
			if (d.has("index_count"))
				index_count = d["index_count"];

			ERR_FAIL_COND_V(!d.has("aabb"), false);
			AABB aabb = d["aabb"];

			{
				// add_surface(VoxelPrimitiveType(primitive), array_data, vertex_count, array_index_data, index_count, aabb);

				VoxelSurface s;
				s.aabb = aabb;
				surfaces.push_back(s);
				_recompute_aabb();

				VisualServer::get_singleton()->voxel_mesh_add_surface(mesh, (VS::VoxelPrimitiveType)primitive, array_data, vertex_count, array_index_data, index_count, aabb);
			}
		} else {
			ERR_FAIL_V(false);
		}

		if (d.has("material")) {

			surface_set_material(idx, d["material"]);
		}
		if (d.has("name")) {
			surface_set_name(idx, d["name"]);
		}

		return true;
	}
	return false;
}

void VoxelMesh::_get_property_list(List<PropertyInfo> *p_list) const {
	for (int i = 0; i < surfaces.size(); i++) {
		p_list->push_back(PropertyInfo(Variant::DICTIONARY, "surfaces/" + itos(i), PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NOEDITOR | PROPERTY_USAGE_INTERNAL));
		p_list->push_back(PropertyInfo(Variant::STRING, "surface_" + itos(i + 1) + "/name", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR));
		p_list->push_back(PropertyInfo(Variant::OBJECT, "surface_" + itos(i + 1) + "/material", PROPERTY_HINT_RESOURCE_TYPE, "ShaderMaterial,SpatialMaterial", PROPERTY_USAGE_EDITOR));
	}
}

/// MARK: - Overrides

RID VoxelMesh::get_rid() const {
	return mesh;
}

void VoxelMesh::reload_from_file() {
	VisualServer::get_singleton()->voxel_mesh_clear(mesh);
	surfaces.clear();

	Resource::reload_from_file();

	_change_notify();
}

/// MARK: - Public

AABB VoxelMesh::get_aabb() const {
	return aabb;
}

int VoxelMesh::get_surface_count() const {
	return surfaces.size();
}

Array VoxelMesh::surface_get_arrays(int p_surface) const {
	ERR_FAIL_INDEX_V(p_surface, surfaces.size(), Array());
	return VisualServer::get_singleton()->voxel_mesh_surface_get_arrays(mesh, p_surface);
	return Array();
}

void VoxelMesh::add_surface_from_arrays(const Array &p_arrays, VoxelPrimitiveType p_primitive, const int p_uv_size) {
	ERR_FAIL_COND(p_arrays.size() != VOXEL_ARRAY_MAX);

	VoxelSurface s;
	VisualServer::get_singleton()->voxel_mesh_add_surface_from_arrays(mesh, (VisualServer::VoxelPrimitiveType)p_primitive, p_arrays, p_uv_size);

	{ // make aABB?

		Variant arr = p_arrays[VOXEL_ARRAY_VERTEX];
		PoolVector<Vector3> vertices = arr;
		int len = vertices.size();
		ERR_FAIL_COND(len == 0);
		PoolVector<Vector3>::Read r = vertices.read();
		const Vector3 *vtx = r.ptr();

		// check AABB
		AABB aabb;
		for (int i = 0; i < len; i++) {

			if (i == 0)
				aabb.position = vtx[i];
			else
				aabb.expand_to(vtx[i]);
		}

		s.aabb = aabb;
		surfaces.push_back(s);

		_recompute_aabb();
	}

	_change_notify();
	emit_changed();
}
void VoxelMesh::surface_remove(int p_idx) {
	ERR_FAIL_INDEX(p_idx, surfaces.size());
	VisualServer::get_singleton()->voxel_mesh_remove_surface(mesh, p_idx);
	surfaces.remove(p_idx);

	_recompute_aabb();
	_change_notify();
	emit_changed();
}

void VoxelMesh::surface_update_region(int p_surface, int p_offset, const PoolVector<uint8_t> &p_data) {
	ERR_FAIL_INDEX(p_surface, surfaces.size());
	VS::get_singleton()->voxel_mesh_surface_update_region(mesh, p_surface, p_offset, p_data);
	emit_changed();
}

/// MARK: Material

Ref<Material> VoxelMesh::surface_get_material(int p_idx) const {
	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), Ref<Material>());
	return surfaces[p_idx].material;
}

void VoxelMesh::surface_set_material(int p_idx, const Ref<Material> &p_material) {
	ERR_FAIL_INDEX(p_idx, surfaces.size());
	if (surfaces[p_idx].material == p_material)
		return;
	surfaces.write[p_idx].material = p_material;
	VisualServer::get_singleton()->voxel_mesh_surface_set_material(mesh, p_idx, p_material.is_null() ? RID() : p_material->get_rid());

	_change_notify("material");
	emit_changed();
}

/// MARK: Name

String VoxelMesh::surface_get_name(int p_idx) const {
	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), String());
	return surfaces[p_idx].name;
}

void VoxelMesh::surface_set_name(int p_idx, const String &p_name) {
	ERR_FAIL_INDEX(p_idx, surfaces.size());
	surfaces.write[p_idx].name = p_name;
	emit_changed();
}

int VoxelMesh::surface_find_by_name(const String &p_name) const {
	for (int i = 0; i < surfaces.size(); i++) {
		if (surfaces[i].name == p_name) {
			return i;
		}
	}
	return -1;
}

// void VoxelMesh::add_surface_from_mesh_data(const Geometry::MeshData &p_mesh_data) {
// 
// 	VisualServer::get_singleton()->mesh_add_surface_from_mesh_data(mesh, p_mesh_data);
// 	AABB aabb;
// 	for (int i = 0; i < p_mesh_data.vertices.size(); i++) {
// 
// 		if (i == 0)
// 			aabb.position = p_mesh_data.vertices[i];
// 		else
// 			aabb.expand_to(p_mesh_data.vertices[i]);
// 	}
// 
// 	Surface s;
// 	s.aabb = aabb;
// 	if (surfaces.size() == 0)
// 		aabb = s.aabb;
// 	else
// 		aabb.merge_with(s.aabb);
// 
// 	clear_cache();
// 
// 	surfaces.push_back(s);
// 	_change_notify();
// 
// 	emit_changed();
// }

// void VoxelMesh::regen_normalmaps() {
// 
// 	Vector<Ref<SurfaceTool> > surfs;
// 	for (int i = 0; i < get_surface_count(); i++) {
// 
// 		Ref<SurfaceTool> st = memnew(SurfaceTool);
// 		st->create_from(Ref<VoxelMesh>(this), i);
// 		surfs.push_back(st);
// 	}
// 
// 	while (get_surface_count()) {
// 		surface_remove(0);
// 	}
// 
// 	for (int i = 0; i < surfs.size(); i++) {
// 
// 		surfs.write[i]->generate_tangents();
// 		surfs.write[i]->commit(Ref<VoxelMesh>(this));
// 	}
// }

// //dirty hack
// bool (*voxel_mesh_lightmap_unwrap_callback)(float p_texel_size, const float *p_vertices, const float *p_normals, int p_vertex_count, const int *p_indices, const int *p_face_materials, int p_index_count, float **r_uv, int **r_vertex, int *r_vertex_count, int **r_index, int *r_index_count, int *r_size_hint_x, int *r_size_hint_y) = NULL;

// struct VoxelMeshLightmapSurface {

// 	Ref<Material> material;
// 	Vector<SurfaceTool::Vertex> vertices;
// 	Mesh::PrimitiveType primitive;
// 	uint32_t format;
// };

// Error VoxelMesh::lightmap_unwrap(const Transform &p_base_transform, float p_texel_size) {

// 	ERR_FAIL_COND_V(!voxel_mesh_lightmap_unwrap_callback, ERR_UNCONFIGURED);

// 	Vector<float> vertices;
// 	Vector<float> normals;
// 	Vector<int> indices;
// 	Vector<int> face_materials;
// 	Vector<float> uv;
// 	Vector<Pair<int, int> > uv_index;

// 	Vector<VoxelMeshLightmapSurface> surfaces;
// 	for (int i = 0; i < get_surface_count(); i++) {
// 		VoxelMeshLightmapSurface s;
// 		s.primitive = surface_get_primitive_type(i);

// 		ERR_FAIL_COND_V_MSG(s.primitive != Mesh::PRIMITIVE_TRIANGLES, ERR_UNAVAILABLE, "Only triangles are supported for lightmap unwrap.");
// 		s.format = surface_get_format(i);
// 		ERR_FAIL_COND_V_MSG(!(s.format & ARRAY_FORMAT_NORMAL), ERR_UNAVAILABLE, "Normals are required for lightmap unwrap.");

// 		Array arrays = surface_get_arrays(i);
// 		s.material = surface_get_material(i);
// 		s.vertices = SurfaceTool::create_vertex_array_from_triangle_arrays(arrays);

// 		PoolVector<Vector3> rvertices = arrays[Mesh::ARRAY_VERTEX];
// 		int vc = rvertices.size();
// 		PoolVector<Vector3>::Read r = rvertices.read();

// 		PoolVector<Vector3> rnormals = arrays[Mesh::ARRAY_NORMAL];
// 		PoolVector<Vector3>::Read rn = rnormals.read();

// 		int vertex_ofs = vertices.size() / 3;

// 		vertices.resize((vertex_ofs + vc) * 3);
// 		normals.resize((vertex_ofs + vc) * 3);
// 		uv_index.resize(vertex_ofs + vc);

// 		for (int j = 0; j < vc; j++) {

// 			Vector3 v = p_base_transform.xform(r[j]);
// 			Vector3 n = p_base_transform.basis.xform(rn[j]).normalized();

// 			vertices.write[(j + vertex_ofs) * 3 + 0] = v.x;
// 			vertices.write[(j + vertex_ofs) * 3 + 1] = v.y;
// 			vertices.write[(j + vertex_ofs) * 3 + 2] = v.z;
// 			normals.write[(j + vertex_ofs) * 3 + 0] = n.x;
// 			normals.write[(j + vertex_ofs) * 3 + 1] = n.y;
// 			normals.write[(j + vertex_ofs) * 3 + 2] = n.z;
// 			uv_index.write[j + vertex_ofs] = Pair<int, int>(i, j);
// 		}

// 		PoolVector<int> rindices = arrays[Mesh::ARRAY_INDEX];
// 		int ic = rindices.size();

// 		if (ic == 0) {

// 			for (int j = 0; j < vc / 3; j++) {
// 				if (Face3(r[j * 3 + 0], r[j * 3 + 1], r[j * 3 + 2]).is_degenerate())
// 					continue;

// 				indices.push_back(vertex_ofs + j * 3 + 0);
// 				indices.push_back(vertex_ofs + j * 3 + 1);
// 				indices.push_back(vertex_ofs + j * 3 + 2);
// 				face_materials.push_back(i);
// 			}

// 		} else {
// 			PoolVector<int>::Read ri = rindices.read();

// 			for (int j = 0; j < ic / 3; j++) {
// 				if (Face3(r[ri[j * 3 + 0]], r[ri[j * 3 + 1]], r[ri[j * 3 + 2]]).is_degenerate())
// 					continue;
// 				indices.push_back(vertex_ofs + ri[j * 3 + 0]);
// 				indices.push_back(vertex_ofs + ri[j * 3 + 1]);
// 				indices.push_back(vertex_ofs + ri[j * 3 + 2]);
// 				face_materials.push_back(i);
// 			}
// 		}

// 		surfaces.push_back(s);
// 	}

// 	//unwrap

// 	float *gen_uvs;
// 	int *gen_vertices;
// 	int *gen_indices;
// 	int gen_vertex_count;
// 	int gen_index_count;
// 	int size_x;
// 	int size_y;

// 	bool ok = voxel_mesh_lightmap_unwrap_callback(p_texel_size, vertices.ptr(), normals.ptr(), vertices.size() / 3, indices.ptr(), face_materials.ptr(), indices.size(), &gen_uvs, &gen_vertices, &gen_vertex_count, &gen_indices, &gen_index_count, &size_x, &size_y);

// 	if (!ok) {
// 		return ERR_CANT_CREATE;
// 	}

// 	//remove surfaces
// 	while (get_surface_count()) {
// 		surface_remove(0);
// 	}

// 	//create surfacetools for each surface..
// 	Vector<Ref<SurfaceTool> > surfaces_tools;

// 	for (int i = 0; i < surfaces.size(); i++) {
// 		Ref<SurfaceTool> st;
// 		st.instance();
// 		st->begin(Mesh::PRIMITIVE_TRIANGLES);
// 		st->set_material(surfaces[i].material);
// 		surfaces_tools.push_back(st); //stay there
// 	}

// 	print_verbose("Mesh: Gen indices: " + itos(gen_index_count));
// 	//go through all indices
// 	for (int i = 0; i < gen_index_count; i += 3) {

// 		ERR_FAIL_INDEX_V(gen_vertices[gen_indices[i + 0]], uv_index.size(), ERR_BUG);
// 		ERR_FAIL_INDEX_V(gen_vertices[gen_indices[i + 1]], uv_index.size(), ERR_BUG);
// 		ERR_FAIL_INDEX_V(gen_vertices[gen_indices[i + 2]], uv_index.size(), ERR_BUG);

// 		ERR_FAIL_COND_V(uv_index[gen_vertices[gen_indices[i + 0]]].first != uv_index[gen_vertices[gen_indices[i + 1]]].first || uv_index[gen_vertices[gen_indices[i + 0]]].first != uv_index[gen_vertices[gen_indices[i + 2]]].first, ERR_BUG);

// 		int surface = uv_index[gen_vertices[gen_indices[i + 0]]].first;

// 		for (int j = 0; j < 3; j++) {

// 			SurfaceTool::Vertex v = surfaces[surface].vertices[uv_index[gen_vertices[gen_indices[i + j]]].second];

// 			if (surfaces[surface].format & ARRAY_FORMAT_COLOR) {
// 				surfaces_tools.write[surface]->add_color(v.color);
// 			}
// 			if (surfaces[surface].format & ARRAY_FORMAT_TEX_UV) {
// 				surfaces_tools.write[surface]->add_uv(v.uv);
// 			}
// 			if (surfaces[surface].format & ARRAY_FORMAT_NORMAL) {
// 				surfaces_tools.write[surface]->add_normal(v.normal);
// 			}
// 			if (surfaces[surface].format & ARRAY_FORMAT_TANGENT) {
// 				Plane t;
// 				t.normal = v.tangent;
// 				t.d = v.binormal.dot(v.normal.cross(v.tangent)) < 0 ? -1 : 1;
// 				surfaces_tools.write[surface]->add_tangent(t);
// 			}
// 			if (surfaces[surface].format & ARRAY_FORMAT_BONES) {
// 				surfaces_tools.write[surface]->add_bones(v.bones);
// 			}
// 			if (surfaces[surface].format & ARRAY_FORMAT_WEIGHTS) {
// 				surfaces_tools.write[surface]->add_weights(v.weights);
// 			}

// 			Vector2 uv2(gen_uvs[gen_indices[i + j] * 2 + 0], gen_uvs[gen_indices[i + j] * 2 + 1]);
// 			surfaces_tools.write[surface]->add_uv2(uv2);

// 			surfaces_tools.write[surface]->add_vertex(v.vertex);
// 		}
// 	}

// 	//free stuff
// 	::free(gen_vertices);
// 	::free(gen_indices);
// 	::free(gen_uvs);

// 	//generate surfaces

// 	for (int i = 0; i < surfaces_tools.size(); i++) {
// 		surfaces_tools.write[i]->index();
// 		surfaces_tools.write[i]->commit(Ref<VoxelMesh>((VoxelMesh *)this), surfaces[i].format);
// 	}

// 	set_lightmap_size_hint(Size2(size_x, size_y));

// 	return OK;
// }

// void VoxelMesh::_bind_methods() {
// 	ClassDB::bind_method(D_METHOD("add_surface_from_arrays", "primitive", "arrays", "compress_flags"), &VoxelMesh::add_surface_from_arrays, DEFVAL(ARRAY_COMPRESS_DEFAULT));
// 	ClassDB::bind_method(D_METHOD("surface_remove", "surf_idx"), &VoxelMesh::surface_remove);
// 	ClassDB::bind_method(D_METHOD("surface_update_region", "surf_idx", "offset", "data"), &VoxelMesh::surface_update_region);
// 	ClassDB::bind_method(D_METHOD("surface_get_array_len", "surf_idx"), &VoxelMesh::surface_get_array_len);
// 	ClassDB::bind_method(D_METHOD("surface_get_array_index_len", "surf_idx"), &VoxelMesh::surface_get_array_index_len);
// 	ClassDB::bind_method(D_METHOD("surface_get_format", "surf_idx"), &VoxelMesh::surface_get_format);
// 	ClassDB::bind_method(D_METHOD("surface_get_primitive_type", "surf_idx"), &VoxelMesh::surface_get_primitive_type);
// 	ClassDB::bind_method(D_METHOD("surface_find_by_name", "name"), &VoxelMesh::surface_find_by_name);
// 	ClassDB::bind_method(D_METHOD("surface_set_name", "surf_idx", "name"), &VoxelMesh::surface_set_name);
// 	ClassDB::bind_method(D_METHOD("surface_get_name", "surf_idx"), &VoxelMesh::surface_get_name);
// 	ClassDB::bind_method(D_METHOD("regen_normalmaps"), &VoxelMesh::regen_normalmaps);
// 	ClassDB::set_method_flags(get_class_static(), _scs_create("regen_normalmaps"), METHOD_FLAGS_DEFAULT | METHOD_FLAG_EDITOR);
// 	ClassDB::bind_method(D_METHOD("lightmap_unwrap", "transform", "texel_size"), &VoxelMesh::lightmap_unwrap);
// 	ClassDB::set_method_flags(get_class_static(), _scs_create("lightmap_unwrap"), METHOD_FLAGS_DEFAULT | METHOD_FLAG_EDITOR);

// 	ClassDB::bind_method(D_METHOD("set_custom_aabb", "aabb"), &VoxelMesh::set_custom_aabb);
// 	ClassDB::bind_method(D_METHOD("get_custom_aabb"), &VoxelMesh::get_custom_aabb);

// 	ADD_PROPERTY(PropertyInfo(Variant::AABB, "custom_aabb", PROPERTY_HINT_NONE, ""), "set_custom_aabb", "get_custom_aabb");

// 	BIND_CONSTANT(NO_INDEX_ARRAY);
// 	BIND_CONSTANT(ARRAY_WEIGHTS_SIZE);

// 	BIND_ENUM_CONSTANT(ARRAY_VERTEX);
// 	BIND_ENUM_CONSTANT(ARRAY_NORMAL);
// 	BIND_ENUM_CONSTANT(ARRAY_TANGENT);
// 	BIND_ENUM_CONSTANT(ARRAY_COLOR);
// 	BIND_ENUM_CONSTANT(ARRAY_TEX_UV);
// 	BIND_ENUM_CONSTANT(ARRAY_TEX_UV2);
// 	BIND_ENUM_CONSTANT(ARRAY_BONES);
// 	BIND_ENUM_CONSTANT(ARRAY_WEIGHTS);
// 	BIND_ENUM_CONSTANT(ARRAY_INDEX);
// 	BIND_ENUM_CONSTANT(ARRAY_MAX);

// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_VERTEX);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_NORMAL);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_TANGENT);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_COLOR);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_TEX_UV);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_TEX_UV2);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_BONES);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_WEIGHTS);
// 	BIND_ENUM_CONSTANT(ARRAY_FORMAT_INDEX);
// }
