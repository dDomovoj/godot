/*************************************************************************/
/*  material.cpp                                                         */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                    http://www.godotengine.org                         */
/*************************************************************************/
/* Copyright (c) 2007-2016 Juan Linietsky, Ariel Manzur.                 */
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
#include "material.h"
#include "scene/scene_string_names.h"

RID Material::get_rid() const {

	return material;
}

Material::Material() {

	material=VisualServer::get_singleton()->material_create();
}

Material::~Material() {

	VisualServer::get_singleton()->free(material);
}


/////////////////////////////////

Mutex *FixedSpatialMaterial::material_mutex=NULL;
SelfList<FixedSpatialMaterial>::List FixedSpatialMaterial::dirty_materials;
Map<FixedSpatialMaterial::MaterialKey,FixedSpatialMaterial::ShaderData> FixedSpatialMaterial::shader_map;
FixedSpatialMaterial::ShaderNames* FixedSpatialMaterial::shader_names=NULL;

void FixedSpatialMaterial::init_shaders() {

#ifndef NO_THREADS
	material_mutex = Mutex::create();
#endif

	shader_names = memnew( ShaderNames );

	shader_names->albedo="albedo";
	shader_names->specular="specular";
	shader_names->roughness="roughness";
	shader_names->emission="emission";
	shader_names->normal_scale="normal_scale";
	shader_names->sheen="sheen";
	shader_names->sheen_color="sheen_color";
	shader_names->clearcoat="clearcoat";
	shader_names->clearcoat_gloss="clearcoat_gloss";
	shader_names->anisotropy="anisotropy";
	shader_names->height_scale="height_scale";
	shader_names->subsurface_scattering="subsurface_scattering";
	shader_names->refraction="refraction";
	shader_names->refraction_roughness="refraction_roughness";
	shader_names->point_size="point_size";

	shader_names->texture_names[TEXTURE_ALBEDO]="texture_albedo";
	shader_names->texture_names[TEXTURE_SPECULAR]="texture_specular";
	shader_names->texture_names[TEXTURE_EMISSION]="texture_emission";
	shader_names->texture_names[TEXTURE_NORMAL]="texture_normal";
	shader_names->texture_names[TEXTURE_SHEEN]="texture_sheen";
	shader_names->texture_names[TEXTURE_CLEARCOAT]="texture_clearcoat";
	shader_names->texture_names[TEXTURE_ANISOTROPY]="texture_anisotropy";
	shader_names->texture_names[TEXTURE_AMBIENT_OCCLUSION]="texture_ambient_occlusion";
	shader_names->texture_names[TEXTURE_HEIGHT]="texture_height";
	shader_names->texture_names[TEXTURE_SUBSURFACE_SCATTERING]="texture_subsurface_scattering";
	shader_names->texture_names[TEXTURE_REFRACTION]="texture_refraction";
	shader_names->texture_names[TEXTURE_REFRACTION_ROUGHNESS]="texture_refraction_roughness";
	shader_names->texture_names[TEXTURE_DETAIL_MASK]="texture_detail_mask";
	shader_names->texture_names[TEXTURE_DETAIL_ALBEDO]="texture_detail_albedo";
	shader_names->texture_names[TEXTURE_DETAIL_NORMAL]="texture_detail_normal";

}

void FixedSpatialMaterial::finish_shaders(){

#ifndef NO_THREADS
	memdelete( material_mutex );
#endif

	memdelete( shader_names );

}



void FixedSpatialMaterial::_update_shader() {

	dirty_materials.remove( &element );

	MaterialKey mk = _compute_key();
	if (mk.key==current_key.key)
		return; //no update required in the end

	if (shader_map.has(current_key)) {
		shader_map[current_key].users--;
		if (shader_map[current_key].users==0) {
			//deallocate shader, as it's no longer in use
			VS::get_singleton()->free(shader_map[current_key].shader);
			shader_map.erase(current_key);
		}
	}

	current_key=mk;

	if (shader_map.has(mk)) {

		VS::get_singleton()->material_set_shader(_get_material(),shader_map[mk].shader);
		shader_map[mk].users++;
		return;
	}

	//must create a shader!

	String code="render_mode ";
	switch(blend_mode) {
		case BLEND_MODE_MIX: code+="blend_mix"; break;
		case BLEND_MODE_ADD: code+="blend_add"; break;
		case BLEND_MODE_SUB: code+="blend_sub"; break;
		case BLEND_MODE_MUL: code+="blend_mul"; break;
	}

	switch(depth_draw_mode) {
		case DEPTH_DRAW_OPAQUE_ONLY: code+=",depth_draw_opaque"; break;
		case DEPTH_DRAW_ALWAYS: code+=",depth_draw_always"; break;
		case DEPTH_DRAW_DISABLED: code+=",depth_draw_never"; break;
		case DEPTH_DRAW_ALPHA_OPAQUE_PREPASS: code+=",depth_draw_alpha_prepass"; break;
	}

	switch(cull_mode) {
		case CULL_BACK: code+=",cull_back"; break;
		case CULL_FRONT: code+=",cull_front"; break;
		case CULL_DISABLED: code+=",cull_disabled"; break;

	}

	if (flags[FLAG_UNSHADED]) {
		code+=",unshaded";
	}
	if (flags[FLAG_ONTOP]) {
		code+=",ontop";
	}

	code+=";\n";


	code+="uniform vec4 albedo : hint_color;\n";
	code+="uniform sampler2D texture_albedo : hint_albedo;\n";
	code+="uniform vec4 specular : hint_color;\n";
	code+="uniform float roughness : hint_range(0,1);\n";
	code+="uniform float point_size : hint_range(0,128);\n";
	code+="uniform sampler2D texture_specular : hint_white;\n";
	code+="\n\n";

	code+="void vertex() {\n";

	if (flags[FLAG_SRGB_VERTEX_COLOR]) {

		code+="\tCOLOR.rgb = mix( pow((COLOR.rgb + vec3(0.055)) * (1.0 / (1.0 + 0.055)), vec3(2.4)), COLOR.rgb* (1.0 / 12.92), lessThan(COLOR.rgb,vec3(0.04045)) );\n";
	}
	if (flags[FLAG_USE_POINT_SIZE]) {

		code+="\tPOINT_SIZE=point_size;\n";
	}

	code+="}\n";
	code+="\n\n";
	code+="void fragment() {\n";

	if (flags[FLAG_USE_POINT_SIZE]) {
		code+="\tvec4 albedo_tex = texture(texture_albedo,POINT_COORD);\n";
	} else {
		code+="\tvec4 albedo_tex = texture(texture_albedo,UV);\n";
	}

	if (flags[FLAG_ALBEDO_FROM_VERTEX_COLOR]) {
		code+="\talbedo_tex *= COLOR;\n";
	}

	code+="\tALBEDO = albedo.rgb * albedo_tex.rgb;\n";
	if (features[FEATURE_TRANSPARENT]) {
		code+="\tALPHA = albedo.a * albedo_tex.a;\n";
	}
	code+="\tvec4 specular_tex = texture(texture_specular,UV);\n";
	code+="\tSPECULAR = specular.rgb * specular_tex.rgb;\n";
	code+="\tROUGHNESS = specular.a * roughness;\n";
	code+="}\n";

	ShaderData shader_data;
	shader_data.shader = VS::get_singleton()->shader_create(VS::SHADER_SPATIAL);
	shader_data.users=1;

	VS::get_singleton()->shader_set_code( shader_data.shader, code );

	shader_map[mk]=shader_data;

	VS::get_singleton()->material_set_shader(_get_material(),shader_data.shader);

}

void FixedSpatialMaterial::flush_changes() {

	if (material_mutex)
		material_mutex->lock();

	while (dirty_materials.first()) {

		dirty_materials.first()->self()->_update_shader();
	}

	if (material_mutex)
		material_mutex->unlock();
}

void FixedSpatialMaterial::_queue_shader_change() {

	if (material_mutex)
		material_mutex->lock();

	if (!element.in_list()) {
		dirty_materials.add(&element);
	}

	if (material_mutex)
		material_mutex->unlock();


}

bool FixedSpatialMaterial::_is_shader_dirty() const {

	bool dirty=false;

	if (material_mutex)
		material_mutex->lock();

	dirty=element.in_list();

	if (material_mutex)
		material_mutex->unlock();

	return dirty;
}
void FixedSpatialMaterial::set_albedo(const Color& p_albedo) {

	albedo=p_albedo;

	VS::get_singleton()->material_set_param(_get_material(),shader_names->albedo,p_albedo);
}

Color FixedSpatialMaterial::get_albedo() const{

	return albedo;
}

void FixedSpatialMaterial::set_specular(const Color& p_specular){

	specular=p_specular;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->specular,p_specular);

}
Color FixedSpatialMaterial::get_specular() const{

	return specular;
}

void FixedSpatialMaterial::set_roughness(float p_roughness){

	roughness=p_roughness;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->roughness,p_roughness);

}
float FixedSpatialMaterial::get_roughness() const{

	return roughness;
}

void FixedSpatialMaterial::set_emission(const Color& p_emission){

	emission=p_emission;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->emission,p_emission);

}
Color FixedSpatialMaterial::get_emission() const{

	return emission;
}

void FixedSpatialMaterial::set_normal_scale(float p_normal_scale){

	normal_scale=p_normal_scale;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->normal_scale,p_normal_scale);

}
float FixedSpatialMaterial::get_normal_scale() const{

	return normal_scale;
}

void FixedSpatialMaterial::set_sheen(float p_sheen){

	sheen=p_sheen;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->sheen,p_sheen);


}
float FixedSpatialMaterial::get_sheen() const{

	return sheen;
}

void FixedSpatialMaterial::set_sheen_color(const Color& p_sheen_color){

	sheen_color=p_sheen_color;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->sheen_color,p_sheen_color);

}
Color FixedSpatialMaterial::get_sheen_color() const{

	return sheen_color;
}

void FixedSpatialMaterial::set_clearcoat(float p_clearcoat){

	clearcoat=p_clearcoat;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->clearcoat,p_clearcoat);

}

float FixedSpatialMaterial::get_clearcoat() const{

	return clearcoat;
}

void FixedSpatialMaterial::set_clearcoat_gloss(float p_clearcoat_gloss){

	clearcoat_gloss=p_clearcoat_gloss;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->clearcoat_gloss,p_clearcoat_gloss);


}

float FixedSpatialMaterial::get_clearcoat_gloss() const{

	return clearcoat_gloss;
}

void FixedSpatialMaterial::set_anisotropy(float p_anisotropy){

	anisotropy=p_anisotropy;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->anisotropy,p_anisotropy);

}
float FixedSpatialMaterial::get_anisotropy() const{

	return anisotropy;
}

void FixedSpatialMaterial::set_height_scale(float p_height_scale){

	height_scale=p_height_scale;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->height_scale,p_height_scale);


}

float FixedSpatialMaterial::get_height_scale() const{

	return height_scale;
}

void FixedSpatialMaterial::set_subsurface_scattering(float p_subsurface_scattering){

	subsurface_scattering=p_subsurface_scattering;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->subsurface_scattering,subsurface_scattering);


}

float FixedSpatialMaterial::get_subsurface_scattering() const{

	return subsurface_scattering;
}

void FixedSpatialMaterial::set_refraction(float p_refraction){

	refraction=p_refraction;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->refraction,refraction);

}

float FixedSpatialMaterial::get_refraction() const {

	return refraction;
}

void FixedSpatialMaterial::set_refraction_roughness(float p_refraction_roughness) {

	refraction_roughness=p_refraction_roughness;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->refraction_roughness,refraction_roughness);


}
float FixedSpatialMaterial::get_refraction_roughness() const {

	return refraction_roughness;
}

void FixedSpatialMaterial::set_detail_uv(DetailUV p_detail_uv) {

	if (detail_uv==p_detail_uv)
		return;

	detail_uv=p_detail_uv;
	_queue_shader_change();
}
FixedSpatialMaterial::DetailUV FixedSpatialMaterial::get_detail_uv() const {

	return detail_uv;
}

void FixedSpatialMaterial::set_blend_mode(BlendMode p_mode) {

	if (blend_mode==p_mode)
		return;

	blend_mode=p_mode;
	_queue_shader_change();
}
FixedSpatialMaterial::BlendMode FixedSpatialMaterial::get_blend_mode() const {

	return blend_mode;
}

void FixedSpatialMaterial::set_detail_blend_mode(BlendMode p_mode) {

	detail_blend_mode=p_mode;
	_queue_shader_change();
}
FixedSpatialMaterial::BlendMode FixedSpatialMaterial::get_detail_blend_mode() const {

	return detail_blend_mode;
}

void FixedSpatialMaterial::set_depth_draw_mode(DepthDrawMode p_mode) {

	if (depth_draw_mode==p_mode)
		return;

	depth_draw_mode=p_mode;
	_queue_shader_change();
}
FixedSpatialMaterial::DepthDrawMode FixedSpatialMaterial::get_depth_draw_mode() const {

	return depth_draw_mode;
}

void FixedSpatialMaterial::set_cull_mode(CullMode p_mode) {

	if (cull_mode==p_mode)
		return;

	cull_mode=p_mode;
	_queue_shader_change();
}
FixedSpatialMaterial::CullMode FixedSpatialMaterial::get_cull_mode() const {

	return cull_mode;
}

void FixedSpatialMaterial::set_diffuse_mode(DiffuseMode p_mode) {

	if (diffuse_mode==p_mode)
		return;

	diffuse_mode=p_mode;
	_queue_shader_change();
}
FixedSpatialMaterial::DiffuseMode FixedSpatialMaterial::get_diffuse_mode() const {

	return diffuse_mode;
}

void FixedSpatialMaterial::set_flag(Flags p_flag,bool p_enabled) {

	ERR_FAIL_INDEX(p_flag,FLAG_MAX);

	if (flags[p_flag]==p_enabled)
		return;

	flags[p_flag]=p_enabled;
	_queue_shader_change();
}

bool FixedSpatialMaterial::get_flag(Flags p_flag) const {

	ERR_FAIL_INDEX_V(p_flag,FLAG_MAX,false);
	return flags[p_flag];
}

void FixedSpatialMaterial::set_feature(Feature p_feature,bool p_enabled) {

	ERR_FAIL_INDEX(p_feature,FEATURE_MAX);
	if (features[p_feature]==p_enabled)
		return;

	features[p_feature]=p_enabled;
	_change_notify();
	_queue_shader_change();


}

bool FixedSpatialMaterial::get_feature(Feature p_feature) const {

	ERR_FAIL_INDEX_V(p_feature,FEATURE_MAX,false);
	return features[p_feature];
}



void FixedSpatialMaterial::set_texture(TextureParam p_param, const Ref<Texture> &p_texture) {

	ERR_FAIL_INDEX(p_param,TEXTURE_MAX);
	textures[p_param]=p_texture;
	RID rid = p_texture.is_valid() ? p_texture->get_rid() : RID();
	VS::get_singleton()->material_set_param(_get_material(),shader_names->texture_names[p_param],rid);
}

Ref<Texture> FixedSpatialMaterial::get_texture(TextureParam p_param) const {

	ERR_FAIL_INDEX_V(p_param,TEXTURE_MAX,Ref<Texture>());
	return textures[p_param];
}


void FixedSpatialMaterial::_validate_feature(const String& text, Feature feature,PropertyInfo& property) const {
	if (property.name.begins_with(text) && property.name!=text+"/enabled" && !features[feature]) {
		property.usage=0;
	}

}

void FixedSpatialMaterial::_validate_property(PropertyInfo& property) const {
	_validate_feature("normal",FEATURE_NORMAL_MAPPING,property);
	_validate_feature("sheen",FEATURE_SHEEN,property);
	_validate_feature("clearcoat",FEATURE_CLEARCOAT,property);
	_validate_feature("anisotropy",FEATURE_ANISOTROPY,property);
	_validate_feature("ao",FEATURE_AMBIENT_OCCLUSION,property);
	_validate_feature("height",FEATURE_HEIGHT_MAPPING,property);
	_validate_feature("subsurf_scatter",FEATURE_SUBSURACE_SCATTERING,property);
	_validate_feature("refraction",FEATURE_REFRACTION,property);
	_validate_feature("detail",FEATURE_DETAIL,property);

}

void FixedSpatialMaterial::set_line_width(float p_line_width) {

	line_width=p_line_width;
	VS::get_singleton()->material_set_line_width(_get_material(),line_width);
}

float FixedSpatialMaterial::get_line_width() const {

	return line_width;
}

void FixedSpatialMaterial::set_point_size(float p_point_size) {

	point_size=p_point_size;
	VS::get_singleton()->material_set_param(_get_material(),shader_names->point_size,p_point_size);
}

float FixedSpatialMaterial::get_point_size() const {

	return point_size;
}

void FixedSpatialMaterial::_bind_methods() {


	ObjectTypeDB::bind_method(_MD("set_albedo","albedo"),&FixedSpatialMaterial::set_albedo);
	ObjectTypeDB::bind_method(_MD("get_albedo"),&FixedSpatialMaterial::get_albedo);

	ObjectTypeDB::bind_method(_MD("set_specular","specular"),&FixedSpatialMaterial::set_specular);
	ObjectTypeDB::bind_method(_MD("get_specular"),&FixedSpatialMaterial::get_specular);

	ObjectTypeDB::bind_method(_MD("set_roughness","roughness"),&FixedSpatialMaterial::set_roughness);
	ObjectTypeDB::bind_method(_MD("get_roughness"),&FixedSpatialMaterial::get_roughness);

	ObjectTypeDB::bind_method(_MD("set_emission","emission"),&FixedSpatialMaterial::set_emission);
	ObjectTypeDB::bind_method(_MD("get_emission"),&FixedSpatialMaterial::get_emission);

	ObjectTypeDB::bind_method(_MD("set_normal_scale","normal_scale"),&FixedSpatialMaterial::set_normal_scale);
	ObjectTypeDB::bind_method(_MD("get_normal_scale"),&FixedSpatialMaterial::get_normal_scale);

	ObjectTypeDB::bind_method(_MD("set_sheen","sheen"),&FixedSpatialMaterial::set_sheen);
	ObjectTypeDB::bind_method(_MD("get_sheen"),&FixedSpatialMaterial::get_sheen);

	ObjectTypeDB::bind_method(_MD("set_sheen_color","sheen_color"),&FixedSpatialMaterial::set_sheen_color);
	ObjectTypeDB::bind_method(_MD("get_sheen_color"),&FixedSpatialMaterial::get_sheen_color);

	ObjectTypeDB::bind_method(_MD("set_clearcoat","clearcoat"),&FixedSpatialMaterial::set_clearcoat);
	ObjectTypeDB::bind_method(_MD("get_clearcoat"),&FixedSpatialMaterial::get_clearcoat);

	ObjectTypeDB::bind_method(_MD("set_clearcoat_gloss","clearcoat_gloss"),&FixedSpatialMaterial::set_clearcoat_gloss);
	ObjectTypeDB::bind_method(_MD("get_clearcoat_gloss"),&FixedSpatialMaterial::get_clearcoat_gloss);

	ObjectTypeDB::bind_method(_MD("set_anisotropy","anisotropy"),&FixedSpatialMaterial::set_anisotropy);
	ObjectTypeDB::bind_method(_MD("get_anisotropy"),&FixedSpatialMaterial::get_anisotropy);

	ObjectTypeDB::bind_method(_MD("set_height_scale","height_scale"),&FixedSpatialMaterial::set_height_scale);
	ObjectTypeDB::bind_method(_MD("get_height_scale"),&FixedSpatialMaterial::get_height_scale);

	ObjectTypeDB::bind_method(_MD("set_subsurface_scattering","subsurface_scattering"),&FixedSpatialMaterial::set_subsurface_scattering);
	ObjectTypeDB::bind_method(_MD("get_subsurface_scattering"),&FixedSpatialMaterial::get_subsurface_scattering);

	ObjectTypeDB::bind_method(_MD("set_refraction","refraction"),&FixedSpatialMaterial::set_refraction);
	ObjectTypeDB::bind_method(_MD("get_refraction"),&FixedSpatialMaterial::get_refraction);

	ObjectTypeDB::bind_method(_MD("set_refraction_roughness","refraction_roughness"),&FixedSpatialMaterial::set_refraction_roughness);
	ObjectTypeDB::bind_method(_MD("get_refraction_roughness"),&FixedSpatialMaterial::get_refraction_roughness);

	ObjectTypeDB::bind_method(_MD("set_line_width","line_width"),&FixedSpatialMaterial::set_line_width);
	ObjectTypeDB::bind_method(_MD("get_line_width"),&FixedSpatialMaterial::get_line_width);

	ObjectTypeDB::bind_method(_MD("set_point_size","point_size"),&FixedSpatialMaterial::set_point_size);
	ObjectTypeDB::bind_method(_MD("get_point_size"),&FixedSpatialMaterial::get_point_size);

	ObjectTypeDB::bind_method(_MD("set_detail_uv","detail_uv"),&FixedSpatialMaterial::set_detail_uv);
	ObjectTypeDB::bind_method(_MD("get_detail_uv"),&FixedSpatialMaterial::get_detail_uv);

	ObjectTypeDB::bind_method(_MD("set_blend_mode","blend_mode"),&FixedSpatialMaterial::set_blend_mode);
	ObjectTypeDB::bind_method(_MD("get_blend_mode"),&FixedSpatialMaterial::get_blend_mode);

	ObjectTypeDB::bind_method(_MD("set_depth_draw_mode","depth_draw_mode"),&FixedSpatialMaterial::set_depth_draw_mode);
	ObjectTypeDB::bind_method(_MD("get_depth_draw_mode"),&FixedSpatialMaterial::get_depth_draw_mode);

	ObjectTypeDB::bind_method(_MD("set_cull_mode","cull_mode"),&FixedSpatialMaterial::set_cull_mode);
	ObjectTypeDB::bind_method(_MD("get_cull_mode"),&FixedSpatialMaterial::get_cull_mode);

	ObjectTypeDB::bind_method(_MD("set_diffuse_mode","diffuse_mode"),&FixedSpatialMaterial::set_diffuse_mode);
	ObjectTypeDB::bind_method(_MD("get_diffuse_mode"),&FixedSpatialMaterial::get_diffuse_mode);

	ObjectTypeDB::bind_method(_MD("set_flag","flag","enable"),&FixedSpatialMaterial::set_flag);
	ObjectTypeDB::bind_method(_MD("get_flag"),&FixedSpatialMaterial::get_flag);

	ObjectTypeDB::bind_method(_MD("set_feature","feature","enable"),&FixedSpatialMaterial::set_feature);
	ObjectTypeDB::bind_method(_MD("get_feature","feature"),&FixedSpatialMaterial::get_feature);

	ObjectTypeDB::bind_method(_MD("set_texture","param:Texture","texture"),&FixedSpatialMaterial::set_texture);
	ObjectTypeDB::bind_method(_MD("get_texture:Texture","param:Texture"),&FixedSpatialMaterial::get_texture);

	ObjectTypeDB::bind_method(_MD("set_detail_blend_mode","detail_blend_mode"),&FixedSpatialMaterial::set_detail_blend_mode);
	ObjectTypeDB::bind_method(_MD("get_detail_blend_mode"),&FixedSpatialMaterial::get_detail_blend_mode);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"flags/transparent"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_TRANSPARENT);
	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"flags/unshaded"),_SCS("set_flag"),_SCS("get_flag"),FLAG_UNSHADED);
	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"flags/on_top"),_SCS("set_flag"),_SCS("get_flag"),FLAG_ONTOP);
	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"flags/use_point_size"),_SCS("set_flag"),_SCS("get_flag"),FLAG_USE_POINT_SIZE);
	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"vertex_color/use_as_albedo"),_SCS("set_flag"),_SCS("get_flag"),FLAG_ALBEDO_FROM_VERTEX_COLOR);
	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"vertex_color/is_srgb"),_SCS("set_flag"),_SCS("get_flag"),FLAG_SRGB_VERTEX_COLOR);

	ADD_PROPERTY(PropertyInfo(Variant::INT,"params/diffuse_mode",PROPERTY_HINT_ENUM,"Labert,Lambert Wrap,Oren Nayar,Burley"),_SCS("set_diffuse_mode"),_SCS("get_diffuse_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::INT,"params/blend_mode",PROPERTY_HINT_ENUM,"Mix,Add,Sub,Mul"),_SCS("set_blend_mode"),_SCS("get_blend_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::INT,"params/cull_mode",PROPERTY_HINT_ENUM,"Back,Front,Disabled"),_SCS("set_cull_mode"),_SCS("get_cull_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::INT,"params/depth_draw_mode",PROPERTY_HINT_ENUM,"Opaque Only,Always,Never,Opaque Pre-Pass"),_SCS("set_depth_draw_mode"),_SCS("get_depth_draw_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"params/line_width",PROPERTY_HINT_RANGE,"0.1,128,0.1"),_SCS("set_line_width"),_SCS("get_line_width"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"params/point_size",PROPERTY_HINT_RANGE,"0.1,128,0.1"),_SCS("set_point_size"),_SCS("get_point_size"));

	ADD_PROPERTY(PropertyInfo(Variant::COLOR,"albedo/color"),_SCS("set_albedo"),_SCS("get_albedo"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"albedo/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_ALBEDO);

	ADD_PROPERTY(PropertyInfo(Variant::COLOR,"specular/color"),_SCS("set_specular"),_SCS("get_specular"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"specular/roughness",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_roughness"),_SCS("get_roughness"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"specular/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_SPECULAR);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"normal/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_NORMAL_MAPPING);
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"normal/scale",PROPERTY_HINT_RANGE,"-16,16,0.01"),_SCS("set_normal_scale"),_SCS("get_normal_scale"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"normal/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_NORMAL);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"sheen/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_SHEEN);
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"sheen/amount",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_sheen"),_SCS("get_sheen"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"sheen/color"),_SCS("set_sheen_color"),_SCS("get_sheen_color"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"sheen/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_SHEEN);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"clearcoat/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_CLEARCOAT);
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"clearcoat/amount",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_clearcoat"),_SCS("get_clearcoat"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"clearcoat/gloss"),_SCS("set_clearcoat_gloss"),_SCS("get_clearcoat_gloss"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"clearcoat/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_CLEARCOAT);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"anisotropy/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_ANISOTROPY);
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"anisotropy/amount",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_anisotropy"),_SCS("get_anisotropy"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"anisotropy/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_ANISOTROPY);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"ao/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_AMBIENT_OCCLUSION);
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"ao/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_AMBIENT_OCCLUSION);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"height/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_HEIGHT_MAPPING);
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"height/scale",PROPERTY_HINT_RANGE,"-16,16,0.01"),_SCS("set_height_scale"),_SCS("get_height_scale"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"height/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_HEIGHT);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"subsurf_scatter/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_SUBSURACE_SCATTERING);
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"subsurf_scatter/amount",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_subsurface_scattering"),_SCS("get_subsurface_scattering"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"subsurf_scatter/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_SUBSURFACE_SCATTERING);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"refraction/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_REFRACTION);
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"refraction/displacement",PROPERTY_HINT_RANGE,"-1,1,0.01"),_SCS("set_refraction"),_SCS("get_refraction"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"refraction/roughness",PROPERTY_HINT_RANGE,"0,1,0.01"),_SCS("set_refraction_roughness"),_SCS("get_refraction_roughness"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"refraction/texture",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_REFRACTION);

	ADD_PROPERTYI(PropertyInfo(Variant::BOOL,"detail/enabled"),_SCS("set_feature"),_SCS("get_feature"),FEATURE_DETAIL);
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"detail/mask",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_DETAIL_MASK);
	ADD_PROPERTY(PropertyInfo(Variant::INT,"detail/blend_mode",PROPERTY_HINT_ENUM,"Mix,Add,Sub,Mul"),_SCS("set_detail_blend_mode"),_SCS("get_detail_blend_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::INT,"detail/uv_layer",PROPERTY_HINT_ENUM,"UV1,UV2"),_SCS("set_detail_uv"),_SCS("get_detail_uv"));
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"detail/albedo",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_DETAIL_ALBEDO);
	ADD_PROPERTYI(PropertyInfo(Variant::OBJECT,"detail/normal",PROPERTY_HINT_RESOURCE_TYPE,"Texture"),_SCS("set_texture"),_SCS("get_texture"),TEXTURE_DETAIL_NORMAL);

	BIND_CONSTANT( TEXTURE_ALBEDO );
	BIND_CONSTANT( TEXTURE_SPECULAR );
	BIND_CONSTANT( TEXTURE_EMISSION );
	BIND_CONSTANT( TEXTURE_NORMAL );
	BIND_CONSTANT( TEXTURE_SHEEN );
	BIND_CONSTANT( TEXTURE_CLEARCOAT );
	BIND_CONSTANT( TEXTURE_ANISOTROPY );
	BIND_CONSTANT( TEXTURE_AMBIENT_OCCLUSION );
	BIND_CONSTANT( TEXTURE_HEIGHT );
	BIND_CONSTANT( TEXTURE_SUBSURFACE_SCATTERING );
	BIND_CONSTANT( TEXTURE_REFRACTION );
	BIND_CONSTANT( TEXTURE_REFRACTION_ROUGHNESS );
	BIND_CONSTANT( TEXTURE_DETAIL_MASK );
	BIND_CONSTANT( TEXTURE_DETAIL_ALBEDO );
	BIND_CONSTANT( TEXTURE_DETAIL_NORMAL );
	BIND_CONSTANT( TEXTURE_MAX );


	BIND_CONSTANT( DETAIL_UV_1 );
	BIND_CONSTANT( DETAIL_UV_2 );

	BIND_CONSTANT( FEATURE_TRANSPARENT );
	BIND_CONSTANT( FEATURE_EMISSION );
	BIND_CONSTANT( FEATURE_NORMAL_MAPPING );
	BIND_CONSTANT( FEATURE_SHEEN );
	BIND_CONSTANT( FEATURE_CLEARCOAT );
	BIND_CONSTANT( FEATURE_ANISOTROPY );
	BIND_CONSTANT( FEATURE_AMBIENT_OCCLUSION );
	BIND_CONSTANT( FEATURE_HEIGHT_MAPPING );
	BIND_CONSTANT( FEATURE_SUBSURACE_SCATTERING );
	BIND_CONSTANT( FEATURE_REFRACTION );
	BIND_CONSTANT( FEATURE_DETAIL );
	BIND_CONSTANT( FEATURE_MAX );

	BIND_CONSTANT( BLEND_MODE_MIX );
	BIND_CONSTANT( BLEND_MODE_ADD );
	BIND_CONSTANT( BLEND_MODE_SUB );
	BIND_CONSTANT( BLEND_MODE_MUL );

	BIND_CONSTANT( DEPTH_DRAW_OPAQUE_ONLY );
	BIND_CONSTANT( DEPTH_DRAW_ALWAYS );
	BIND_CONSTANT( DEPTH_DRAW_DISABLED );
	BIND_CONSTANT( DEPTH_DRAW_ALPHA_OPAQUE_PREPASS );


	BIND_CONSTANT( CULL_BACK );
	BIND_CONSTANT( CULL_FRONT );
	BIND_CONSTANT( CULL_DISABLED );

	BIND_CONSTANT( FLAG_UNSHADED );
	BIND_CONSTANT( FLAG_ONTOP );
	BIND_CONSTANT( FLAG_ALBEDO_FROM_VERTEX_COLOR );
	BIND_CONSTANT( FLAG_SRGB_VERTEX_COLOR )
	BIND_CONSTANT( FLAG_USE_POINT_SIZE )
	BIND_CONSTANT( FLAG_MAX );

	BIND_CONSTANT( DIFFUSE_LAMBERT );
	BIND_CONSTANT( DIFFUSE_LAMBERT_WRAP );
	BIND_CONSTANT( DIFFUSE_OREN_NAYAR );
	BIND_CONSTANT( DIFFUSE_BURLEY );
}


FixedSpatialMaterial::FixedSpatialMaterial() : element(this) {

	//initialize to right values
	set_albedo(Color(0.7,0.7,0.7,1.0));
	set_specular(Color(0.1,0.1,0.1));
	set_roughness(0.0);
	set_emission(Color(0,0,0));
	set_normal_scale(1);
	set_sheen(0);
	set_sheen_color(Color(1,1,1,1));
	set_clearcoat(0);
	set_clearcoat_gloss(0.5);
	set_anisotropy(0);
	set_height_scale(1);
	set_subsurface_scattering(0);
	set_refraction(0);
	set_refraction_roughness(0);
	set_line_width(1);
	set_point_size(1);

	detail_uv=DETAIL_UV_1;
	blend_mode=BLEND_MODE_MIX;
	detail_blend_mode=BLEND_MODE_MIX;
	depth_draw_mode=DEPTH_DRAW_OPAQUE_ONLY;
	cull_mode=CULL_BACK;
	for(int i=0;i<FLAG_MAX;i++) {
		flags[i]=0;
	}
	diffuse_mode=DIFFUSE_LAMBERT;

	for(int i=0;i<FEATURE_MAX;i++) {
		features[i]=false;
	}

	current_key.key=0;
	current_key.invalid_key=1;
	_queue_shader_change();
}

FixedSpatialMaterial::~FixedSpatialMaterial() {

	if (material_mutex)
		material_mutex->lock();

	if (shader_map.has(current_key)) {
		shader_map[current_key].users--;
		if (shader_map[current_key].users==0) {
			//deallocate shader, as it's no longer in use
			VS::get_singleton()->free(shader_map[current_key].shader);
			shader_map.erase(current_key);
		}

		VS::get_singleton()->material_set_shader(_get_material(),RID());
	}


	if (material_mutex)
		material_mutex->unlock();

}
