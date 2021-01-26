#include "voxel_mesh.h"

#include "core/pair.h"
#include "scene/resources/concave_polygon_shape.h"
#include "scene/resources/convex_polygon_shape.h"
#include "surface_tool.h"

#include <stdlib.h>

bool VoxelMesh::_set(const StringName &p_name, const Variant &p_value) {

	String sname = p_name;

	if (p_name == "blend_shape/names") {

		PoolVector<String> sk = p_value;
		int sz = sk.size();
		PoolVector<String>::Read r = sk.read();
		for (int i = 0; i < sz; i++)
			add_blend_shape(r[i]);
		return true;
	}

	if (p_name == "blend_shape/mode") {

		set_blend_shape_mode(BlendShapeMode(int(p_value)));
		return true;
	}

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
			add_surface_from_arrays(PrimitiveType(int(d["primitive"])), d["arrays"], d["morph_arrays"]);

		} else if (d.has("array_data")) {

			PoolVector<uint8_t> array_data = d["array_data"];
			PoolVector<uint8_t> array_index_data;
			if (d.has("array_index_data"))
				array_index_data = d["array_index_data"];

			ERR_FAIL_COND_V(!d.has("format"), false);
			uint32_t format = d["format"];

			uint32_t primitive = d["primitive"];

			ERR_FAIL_COND_V(!d.has("vertex_count"), false);
			int vertex_count = d["vertex_count"];

			int index_count = 0;
			if (d.has("index_count"))
				index_count = d["index_count"];

			Vector<PoolVector<uint8_t> > blend_shapes;

			if (d.has("blend_shape_data")) {
				Array blend_shape_data = d["blend_shape_data"];
				for (int i = 0; i < blend_shape_data.size(); i++) {
					PoolVector<uint8_t> shape = blend_shape_data[i];
					blend_shapes.push_back(shape);
				}
			}

			ERR_FAIL_COND_V(!d.has("aabb"), false);
			AABB aabb = d["aabb"];

			Vector<AABB> bone_aabb;
			if (d.has("skeleton_aabb")) {
				Array baabb = d["skeleton_aabb"];
				bone_aabb.resize(baabb.size());

				for (int i = 0; i < baabb.size(); i++) {
					bone_aabb.write[i] = baabb[i];
				}
			}

			add_surface(format, PrimitiveType(primitive), array_data, vertex_count, array_index_data, index_count, aabb, blend_shapes, bone_aabb);
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

bool VoxelMesh::_get(const StringName &p_name, Variant &r_ret) const {

	if (_is_generated())
		return false;

	String sname = p_name;

	if (p_name == "blend_shape/names") {

		PoolVector<String> sk;
		for (int i = 0; i < blend_shapes.size(); i++)
			sk.push_back(blend_shapes[i]);
		r_ret = sk;
		return true;
	} else if (p_name == "blend_shape/mode") {

		r_ret = get_blend_shape_mode();
		return true;
	} else if (sname.begins_with("surface_")) {

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

	d["array_data"] = VS::get_singleton()->mesh_surface_get_array(mesh, idx);
	d["vertex_count"] = VS::get_singleton()->mesh_surface_get_array_len(mesh, idx);
	d["array_index_data"] = VS::get_singleton()->mesh_surface_get_index_array(mesh, idx);
	d["index_count"] = VS::get_singleton()->mesh_surface_get_array_index_len(mesh, idx);
	d["primitive"] = VS::get_singleton()->mesh_surface_get_primitive_type(mesh, idx);
	d["format"] = VS::get_singleton()->mesh_surface_get_format(mesh, idx);
	d["aabb"] = VS::get_singleton()->mesh_surface_get_aabb(mesh, idx);

	Vector<AABB> skel_aabb = VS::get_singleton()->mesh_surface_get_skeleton_aabb(mesh, idx);
	Array arr;
	arr.resize(skel_aabb.size());
	for (int i = 0; i < skel_aabb.size(); i++) {
		arr[i] = skel_aabb[i];
	}
	d["skeleton_aabb"] = arr;

	Vector<PoolVector<uint8_t> > blend_shape_data = VS::get_singleton()->mesh_surface_get_blend_shapes(mesh, idx);

	Array md;
	for (int i = 0; i < blend_shape_data.size(); i++) {
		md.push_back(blend_shape_data[i]);
	}

	d["blend_shape_data"] = md;

	Ref<Material> m = surface_get_material(idx);
	if (m.is_valid())
		d["material"] = m;
	String n = surface_get_name(idx);
	if (n != "")
		d["name"] = n;

	r_ret = d;

	return true;
}

void VoxelMesh::_get_property_list(List<PropertyInfo> *p_list) const {

	if (_is_generated())
		return;

	if (blend_shapes.size()) {
		p_list->push_back(PropertyInfo(Variant::POOL_STRING_ARRAY, "blend_shape/names", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NOEDITOR | PROPERTY_USAGE_INTERNAL));
		p_list->push_back(PropertyInfo(Variant::INT, "blend_shape/mode", PROPERTY_HINT_ENUM, "Normalized,Relative"));
	}

	for (int i = 0; i < surfaces.size(); i++) {

		p_list->push_back(PropertyInfo(Variant::DICTIONARY, "surfaces/" + itos(i), PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NOEDITOR | PROPERTY_USAGE_INTERNAL));
		p_list->push_back(PropertyInfo(Variant::STRING, "surface_" + itos(i + 1) + "/name", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR));
		if (surfaces[i].is_2d) {
			p_list->push_back(PropertyInfo(Variant::OBJECT, "surface_" + itos(i + 1) + "/material", PROPERTY_HINT_RESOURCE_TYPE, "ShaderMaterial,CanvasItemMaterial", PROPERTY_USAGE_EDITOR));
		} else {
			p_list->push_back(PropertyInfo(Variant::OBJECT, "surface_" + itos(i + 1) + "/material", PROPERTY_HINT_RESOURCE_TYPE, "ShaderMaterial,SpatialMaterial", PROPERTY_USAGE_EDITOR));
		}
	}
}

void VoxelMesh::_recompute_aabb() {

	// regenerate AABB
	aabb = AABB();

	for (int i = 0; i < surfaces.size(); i++) {

		if (i == 0)
			aabb = surfaces[i].aabb;
		else
			aabb.merge_with(surfaces[i].aabb);
	}
}

void VoxelMesh::add_surface(uint32_t p_format, PrimitiveType p_primitive, const PoolVector<uint8_t> &p_array, int p_vertex_count, const PoolVector<uint8_t> &p_index_array, int p_index_count, const AABB &p_aabb, const Vector<PoolVector<uint8_t> > &p_blend_shapes, const Vector<AABB> &p_bone_aabbs) {

	Surface s;
	s.aabb = p_aabb;
	s.is_2d = p_format & ARRAY_FLAG_USE_2D_VERTICES;
	surfaces.push_back(s);
	_recompute_aabb();

	VisualServer::get_singleton()->mesh_add_surface(mesh, p_format, (VS::PrimitiveType)p_primitive, p_array, p_vertex_count, p_index_array, p_index_count, p_aabb, p_blend_shapes, p_bone_aabbs);
}

void VoxelMesh::add_surface_from_arrays(PrimitiveType p_primitive, const Array &p_arrays, const Array &p_blend_shapes, uint32_t p_flags) {

	ERR_FAIL_COND(p_arrays.size() != ARRAY_MAX);

	Surface s;

	VisualServer::get_singleton()->mesh_add_surface_from_arrays(mesh, (VisualServer::PrimitiveType)p_primitive, p_arrays, p_blend_shapes, p_flags);

	/* make aABB? */ {

		Variant arr = p_arrays[ARRAY_VERTEX];
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
		s.is_2d = arr.get_type() == Variant::POOL_VECTOR2_ARRAY;
		surfaces.push_back(s);

		_recompute_aabb();
	}

	clear_cache();
	_change_notify();
	emit_changed();
}

Array VoxelMesh::surface_get_arrays(int p_surface) const {

	ERR_FAIL_INDEX_V(p_surface, surfaces.size(), Array());
	return VisualServer::get_singleton()->mesh_surface_get_arrays(mesh, p_surface);
}
Array VoxelMesh::surface_get_blend_shape_arrays(int p_surface) const {

	ERR_FAIL_INDEX_V(p_surface, surfaces.size(), Array());
	return VisualServer::get_singleton()->mesh_surface_get_blend_shape_arrays(mesh, p_surface);
}

int VoxelMesh::get_surface_count() const {

	return surfaces.size();
}

void VoxelMesh::add_blend_shape(const StringName &p_name) {

	ERR_FAIL_COND_MSG(surfaces.size(), "Can't add a shape key count if surfaces are already created.");

	StringName name = p_name;

	if (blend_shapes.find(name) != -1) {

		int count = 2;
		do {

			name = String(p_name) + " " + itos(count);
			count++;
		} while (blend_shapes.find(name) != -1);
	}

	blend_shapes.push_back(name);
	VS::get_singleton()->mesh_set_blend_shape_count(mesh, blend_shapes.size());
}

int VoxelMesh::get_blend_shape_count() const {

	return blend_shapes.size();
}
StringName VoxelMesh::get_blend_shape_name(int p_index) const {
	ERR_FAIL_INDEX_V(p_index, blend_shapes.size(), StringName());
	return blend_shapes[p_index];
}
void VoxelMesh::clear_blend_shapes() {

	ERR_FAIL_COND_MSG(surfaces.size(), "Can't set shape key count if surfaces are already created.");

	blend_shapes.clear();
}

void VoxelMesh::set_blend_shape_mode(BlendShapeMode p_mode) {

	blend_shape_mode = p_mode;
	VS::get_singleton()->mesh_set_blend_shape_mode(mesh, (VS::BlendShapeMode)p_mode);
}

VoxelMesh::BlendShapeMode VoxelMesh::get_blend_shape_mode() const {

	return blend_shape_mode;
}

void VoxelMesh::surface_remove(int p_idx) {

	ERR_FAIL_INDEX(p_idx, surfaces.size());
	VisualServer::get_singleton()->mesh_remove_surface(mesh, p_idx);
	surfaces.remove(p_idx);

	clear_cache();
	_recompute_aabb();
	_change_notify();
	emit_changed();
}

int VoxelMesh::surface_get_array_len(int p_idx) const {

	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), -1);
	return VisualServer::get_singleton()->mesh_surface_get_array_len(mesh, p_idx);
}

int VoxelMesh::surface_get_array_index_len(int p_idx) const {

	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), -1);
	return VisualServer::get_singleton()->mesh_surface_get_array_index_len(mesh, p_idx);
}

uint32_t VoxelMesh::surface_get_format(int p_idx) const {

	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), 0);
	return VisualServer::get_singleton()->mesh_surface_get_format(mesh, p_idx);
}

VoxelMesh::PrimitiveType VoxelMesh::surface_get_primitive_type(int p_idx) const {

	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), PRIMITIVE_LINES);
	return (PrimitiveType)VisualServer::get_singleton()->mesh_surface_get_primitive_type(mesh, p_idx);
}

void VoxelMesh::surface_set_material(int p_idx, const Ref<Material> &p_material) {

	ERR_FAIL_INDEX(p_idx, surfaces.size());
	if (surfaces[p_idx].material == p_material)
		return;
	surfaces.write[p_idx].material = p_material;
	VisualServer::get_singleton()->mesh_surface_set_material(mesh, p_idx, p_material.is_null() ? RID() : p_material->get_rid());

	_change_notify("material");
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

void VoxelMesh::surface_set_name(int p_idx, const String &p_name) {

	ERR_FAIL_INDEX(p_idx, surfaces.size());

	surfaces.write[p_idx].name = p_name;
	emit_changed();
}

String VoxelMesh::surface_get_name(int p_idx) const {

	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), String());
	return surfaces[p_idx].name;
}

void VoxelMesh::surface_update_region(int p_surface, int p_offset, const PoolVector<uint8_t> &p_data) {

	ERR_FAIL_INDEX(p_surface, surfaces.size());
	VS::get_singleton()->mesh_surface_update_region(mesh, p_surface, p_offset, p_data);
	emit_changed();
}

void VoxelMesh::surface_set_custom_aabb(int p_idx, const AABB &p_aabb) {

	ERR_FAIL_INDEX(p_idx, surfaces.size());
	surfaces.write[p_idx].aabb = p_aabb;
	// set custom aabb too?
	emit_changed();
}

Ref<Material> VoxelMesh::surface_get_material(int p_idx) const {

	ERR_FAIL_INDEX_V(p_idx, surfaces.size(), Ref<Material>());
	return surfaces[p_idx].material;
}

void VoxelMesh::add_surface_from_mesh_data(const Geometry::MeshData &p_mesh_data) {

	VisualServer::get_singleton()->mesh_add_surface_from_mesh_data(mesh, p_mesh_data);
	AABB aabb;
	for (int i = 0; i < p_mesh_data.vertices.size(); i++) {

		if (i == 0)
			aabb.position = p_mesh_data.vertices[i];
		else
			aabb.expand_to(p_mesh_data.vertices[i]);
	}

	Surface s;
	s.aabb = aabb;
	if (surfaces.size() == 0)
		aabb = s.aabb;
	else
		aabb.merge_with(s.aabb);

	clear_cache();

	surfaces.push_back(s);
	_change_notify();

	emit_changed();
}

RID VoxelMesh::get_rid() const {

	return mesh;
}
AABB VoxelMesh::get_aabb() const {

	return aabb;
}

void VoxelMesh::set_custom_aabb(const AABB &p_custom) {

	custom_aabb = p_custom;
	VS::get_singleton()->mesh_set_custom_aabb(mesh, custom_aabb);
	emit_changed();
}

AABB VoxelMesh::get_custom_aabb() const {

	return custom_aabb;
}

void VoxelMesh::regen_normalmaps() {

	Vector<Ref<SurfaceTool> > surfs;
	for (int i = 0; i < get_surface_count(); i++) {

		Ref<SurfaceTool> st = memnew(SurfaceTool);
		st->create_from(Ref<VoxelMesh>(this), i);
		surfs.push_back(st);
	}

	while (get_surface_count()) {
		surface_remove(0);
	}

	for (int i = 0; i < surfs.size(); i++) {

		surfs.write[i]->generate_tangents();
		surfs.write[i]->commit(Ref<VoxelMesh>(this));
	}
}

//dirty hack
bool (*voxel_mesh_lightmap_unwrap_callback)(float p_texel_size, const float *p_vertices, const float *p_normals, int p_vertex_count, const int *p_indices, const int *p_face_materials, int p_index_count, float **r_uv, int **r_vertex, int *r_vertex_count, int **r_index, int *r_index_count, int *r_size_hint_x, int *r_size_hint_y) = NULL;

struct VoxelMeshLightmapSurface {

	Ref<Material> material;
	Vector<SurfaceTool::Vertex> vertices;
	Mesh::PrimitiveType primitive;
	uint32_t format;
};

Error VoxelMesh::lightmap_unwrap(const Transform &p_base_transform, float p_texel_size) {

	ERR_FAIL_COND_V(!voxel_mesh_lightmap_unwrap_callback, ERR_UNCONFIGURED);
	ERR_FAIL_COND_V_MSG(blend_shapes.size() != 0, ERR_UNAVAILABLE, "Can't unwrap mesh with blend shapes.");

	Vector<float> vertices;
	Vector<float> normals;
	Vector<int> indices;
	Vector<int> face_materials;
	Vector<float> uv;
	Vector<Pair<int, int> > uv_index;

	Vector<VoxelMeshLightmapSurface> surfaces;
	for (int i = 0; i < get_surface_count(); i++) {
		VoxelMeshLightmapSurface s;
		s.primitive = surface_get_primitive_type(i);

		ERR_FAIL_COND_V_MSG(s.primitive != Mesh::PRIMITIVE_TRIANGLES, ERR_UNAVAILABLE, "Only triangles are supported for lightmap unwrap.");
		s.format = surface_get_format(i);
		ERR_FAIL_COND_V_MSG(!(s.format & ARRAY_FORMAT_NORMAL), ERR_UNAVAILABLE, "Normals are required for lightmap unwrap.");

		Array arrays = surface_get_arrays(i);
		s.material = surface_get_material(i);
		s.vertices = SurfaceTool::create_vertex_array_from_triangle_arrays(arrays);

		PoolVector<Vector3> rvertices = arrays[Mesh::ARRAY_VERTEX];
		int vc = rvertices.size();
		PoolVector<Vector3>::Read r = rvertices.read();

		PoolVector<Vector3> rnormals = arrays[Mesh::ARRAY_NORMAL];
		PoolVector<Vector3>::Read rn = rnormals.read();

		int vertex_ofs = vertices.size() / 3;

		vertices.resize((vertex_ofs + vc) * 3);
		normals.resize((vertex_ofs + vc) * 3);
		uv_index.resize(vertex_ofs + vc);

		for (int j = 0; j < vc; j++) {

			Vector3 v = p_base_transform.xform(r[j]);
			Vector3 n = p_base_transform.basis.xform(rn[j]).normalized();

			vertices.write[(j + vertex_ofs) * 3 + 0] = v.x;
			vertices.write[(j + vertex_ofs) * 3 + 1] = v.y;
			vertices.write[(j + vertex_ofs) * 3 + 2] = v.z;
			normals.write[(j + vertex_ofs) * 3 + 0] = n.x;
			normals.write[(j + vertex_ofs) * 3 + 1] = n.y;
			normals.write[(j + vertex_ofs) * 3 + 2] = n.z;
			uv_index.write[j + vertex_ofs] = Pair<int, int>(i, j);
		}

		PoolVector<int> rindices = arrays[Mesh::ARRAY_INDEX];
		int ic = rindices.size();

		if (ic == 0) {

			for (int j = 0; j < vc / 3; j++) {
				if (Face3(r[j * 3 + 0], r[j * 3 + 1], r[j * 3 + 2]).is_degenerate())
					continue;

				indices.push_back(vertex_ofs + j * 3 + 0);
				indices.push_back(vertex_ofs + j * 3 + 1);
				indices.push_back(vertex_ofs + j * 3 + 2);
				face_materials.push_back(i);
			}

		} else {
			PoolVector<int>::Read ri = rindices.read();

			for (int j = 0; j < ic / 3; j++) {
				if (Face3(r[ri[j * 3 + 0]], r[ri[j * 3 + 1]], r[ri[j * 3 + 2]]).is_degenerate())
					continue;
				indices.push_back(vertex_ofs + ri[j * 3 + 0]);
				indices.push_back(vertex_ofs + ri[j * 3 + 1]);
				indices.push_back(vertex_ofs + ri[j * 3 + 2]);
				face_materials.push_back(i);
			}
		}

		surfaces.push_back(s);
	}

	//unwrap

	float *gen_uvs;
	int *gen_vertices;
	int *gen_indices;
	int gen_vertex_count;
	int gen_index_count;
	int size_x;
	int size_y;

	bool ok = voxel_mesh_lightmap_unwrap_callback(p_texel_size, vertices.ptr(), normals.ptr(), vertices.size() / 3, indices.ptr(), face_materials.ptr(), indices.size(), &gen_uvs, &gen_vertices, &gen_vertex_count, &gen_indices, &gen_index_count, &size_x, &size_y);

	if (!ok) {
		return ERR_CANT_CREATE;
	}

	//remove surfaces
	while (get_surface_count()) {
		surface_remove(0);
	}

	//create surfacetools for each surface..
	Vector<Ref<SurfaceTool> > surfaces_tools;

	for (int i = 0; i < surfaces.size(); i++) {
		Ref<SurfaceTool> st;
		st.instance();
		st->begin(Mesh::PRIMITIVE_TRIANGLES);
		st->set_material(surfaces[i].material);
		surfaces_tools.push_back(st); //stay there
	}

	print_verbose("Mesh: Gen indices: " + itos(gen_index_count));
	//go through all indices
	for (int i = 0; i < gen_index_count; i += 3) {

		ERR_FAIL_INDEX_V(gen_vertices[gen_indices[i + 0]], uv_index.size(), ERR_BUG);
		ERR_FAIL_INDEX_V(gen_vertices[gen_indices[i + 1]], uv_index.size(), ERR_BUG);
		ERR_FAIL_INDEX_V(gen_vertices[gen_indices[i + 2]], uv_index.size(), ERR_BUG);

		ERR_FAIL_COND_V(uv_index[gen_vertices[gen_indices[i + 0]]].first != uv_index[gen_vertices[gen_indices[i + 1]]].first || uv_index[gen_vertices[gen_indices[i + 0]]].first != uv_index[gen_vertices[gen_indices[i + 2]]].first, ERR_BUG);

		int surface = uv_index[gen_vertices[gen_indices[i + 0]]].first;

		for (int j = 0; j < 3; j++) {

			SurfaceTool::Vertex v = surfaces[surface].vertices[uv_index[gen_vertices[gen_indices[i + j]]].second];

			if (surfaces[surface].format & ARRAY_FORMAT_COLOR) {
				surfaces_tools.write[surface]->add_color(v.color);
			}
			if (surfaces[surface].format & ARRAY_FORMAT_TEX_UV) {
				surfaces_tools.write[surface]->add_uv(v.uv);
			}
			if (surfaces[surface].format & ARRAY_FORMAT_NORMAL) {
				surfaces_tools.write[surface]->add_normal(v.normal);
			}
			if (surfaces[surface].format & ARRAY_FORMAT_TANGENT) {
				Plane t;
				t.normal = v.tangent;
				t.d = v.binormal.dot(v.normal.cross(v.tangent)) < 0 ? -1 : 1;
				surfaces_tools.write[surface]->add_tangent(t);
			}
			if (surfaces[surface].format & ARRAY_FORMAT_BONES) {
				surfaces_tools.write[surface]->add_bones(v.bones);
			}
			if (surfaces[surface].format & ARRAY_FORMAT_WEIGHTS) {
				surfaces_tools.write[surface]->add_weights(v.weights);
			}

			Vector2 uv2(gen_uvs[gen_indices[i + j] * 2 + 0], gen_uvs[gen_indices[i + j] * 2 + 1]);
			surfaces_tools.write[surface]->add_uv2(uv2);

			surfaces_tools.write[surface]->add_vertex(v.vertex);
		}
	}

	//free stuff
	::free(gen_vertices);
	::free(gen_indices);
	::free(gen_uvs);

	//generate surfaces

	for (int i = 0; i < surfaces_tools.size(); i++) {
		surfaces_tools.write[i]->index();
		surfaces_tools.write[i]->commit(Ref<VoxelMesh>((VoxelMesh *)this), surfaces[i].format);
	}

	set_lightmap_size_hint(Size2(size_x, size_y));

	return OK;
}

void VoxelMesh::_bind_methods() {

	ClassDB::bind_method(D_METHOD("add_blend_shape", "name"), &VoxelMesh::add_blend_shape);
	ClassDB::bind_method(D_METHOD("get_blend_shape_count"), &VoxelMesh::get_blend_shape_count);
	ClassDB::bind_method(D_METHOD("get_blend_shape_name", "index"), &VoxelMesh::get_blend_shape_name);
	ClassDB::bind_method(D_METHOD("clear_blend_shapes"), &VoxelMesh::clear_blend_shapes);
	ClassDB::bind_method(D_METHOD("set_blend_shape_mode", "mode"), &VoxelMesh::set_blend_shape_mode);
	ClassDB::bind_method(D_METHOD("get_blend_shape_mode"), &VoxelMesh::get_blend_shape_mode);

	ClassDB::bind_method(D_METHOD("add_surface_from_arrays", "primitive", "arrays", "blend_shapes", "compress_flags"), &VoxelMesh::add_surface_from_arrays, DEFVAL(Array()), DEFVAL(ARRAY_COMPRESS_DEFAULT));
	ClassDB::bind_method(D_METHOD("surface_remove", "surf_idx"), &VoxelMesh::surface_remove);
	ClassDB::bind_method(D_METHOD("surface_update_region", "surf_idx", "offset", "data"), &VoxelMesh::surface_update_region);
	ClassDB::bind_method(D_METHOD("surface_get_array_len", "surf_idx"), &VoxelMesh::surface_get_array_len);
	ClassDB::bind_method(D_METHOD("surface_get_array_index_len", "surf_idx"), &VoxelMesh::surface_get_array_index_len);
	ClassDB::bind_method(D_METHOD("surface_get_format", "surf_idx"), &VoxelMesh::surface_get_format);
	ClassDB::bind_method(D_METHOD("surface_get_primitive_type", "surf_idx"), &VoxelMesh::surface_get_primitive_type);
	ClassDB::bind_method(D_METHOD("surface_find_by_name", "name"), &VoxelMesh::surface_find_by_name);
	ClassDB::bind_method(D_METHOD("surface_set_name", "surf_idx", "name"), &VoxelMesh::surface_set_name);
	ClassDB::bind_method(D_METHOD("surface_get_name", "surf_idx"), &VoxelMesh::surface_get_name);
	// ClassDB::bind_method(D_METHOD("create_trimesh_shape"), &VoxelMesh::create_trimesh_shape);
	// ClassDB::bind_method(D_METHOD("create_convex_shape"), &VoxelMesh::create_convex_shape);
	// ClassDB::bind_method(D_METHOD("create_outline", "margin"), &VoxelMesh::create_outline);
	ClassDB::bind_method(D_METHOD("regen_normalmaps"), &VoxelMesh::regen_normalmaps);
	ClassDB::set_method_flags(get_class_static(), _scs_create("regen_normalmaps"), METHOD_FLAGS_DEFAULT | METHOD_FLAG_EDITOR);
	ClassDB::bind_method(D_METHOD("lightmap_unwrap", "transform", "texel_size"), &VoxelMesh::lightmap_unwrap);
	ClassDB::set_method_flags(get_class_static(), _scs_create("lightmap_unwrap"), METHOD_FLAGS_DEFAULT | METHOD_FLAG_EDITOR);
	// ClassDB::bind_method(D_METHOD("get_faces"), &VoxelMesh::get_faces);
	// ClassDB::bind_method(D_METHOD("generate_triangle_mesh"), &VoxelMesh::generate_triangle_mesh);

	ClassDB::bind_method(D_METHOD("set_custom_aabb", "aabb"), &VoxelMesh::set_custom_aabb);
	ClassDB::bind_method(D_METHOD("get_custom_aabb"), &VoxelMesh::get_custom_aabb);

	ADD_PROPERTY(PropertyInfo(Variant::INT, "blend_shape_mode", PROPERTY_HINT_ENUM, "Normalized,Relative", PROPERTY_USAGE_NOEDITOR), "set_blend_shape_mode", "get_blend_shape_mode");
	ADD_PROPERTY(PropertyInfo(Variant::AABB, "custom_aabb", PROPERTY_HINT_NONE, ""), "set_custom_aabb", "get_custom_aabb");

	BIND_CONSTANT(NO_INDEX_ARRAY);
	BIND_CONSTANT(ARRAY_WEIGHTS_SIZE);

	BIND_ENUM_CONSTANT(ARRAY_VERTEX);
	BIND_ENUM_CONSTANT(ARRAY_NORMAL);
	BIND_ENUM_CONSTANT(ARRAY_TANGENT);
	BIND_ENUM_CONSTANT(ARRAY_COLOR);
	BIND_ENUM_CONSTANT(ARRAY_TEX_UV);
	BIND_ENUM_CONSTANT(ARRAY_TEX_UV2);
	BIND_ENUM_CONSTANT(ARRAY_BONES);
	BIND_ENUM_CONSTANT(ARRAY_WEIGHTS);
	BIND_ENUM_CONSTANT(ARRAY_INDEX);
	BIND_ENUM_CONSTANT(ARRAY_MAX);

	BIND_ENUM_CONSTANT(ARRAY_FORMAT_VERTEX);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_NORMAL);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_TANGENT);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_COLOR);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_TEX_UV);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_TEX_UV2);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_BONES);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_WEIGHTS);
	BIND_ENUM_CONSTANT(ARRAY_FORMAT_INDEX);
}

void VoxelMesh::reload_from_file() {
	VisualServer::get_singleton()->mesh_clear(mesh);
	surfaces.clear();
	clear_blend_shapes();
	clear_cache();

	Resource::reload_from_file();

	_change_notify();
}

VoxelMesh::VoxelMesh() {

	mesh = VisualServer::get_singleton()->mesh_create();
	blend_shape_mode = BLEND_SHAPE_MODE_RELATIVE;
}

VoxelMesh::~VoxelMesh() {

	VisualServer::get_singleton()->free(mesh);
}