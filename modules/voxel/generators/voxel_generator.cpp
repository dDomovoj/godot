#include "voxel_generator.h"
#include "../voxel_string_names.h"

VoxelGenerator::VoxelGenerator() {
}

void VoxelGenerator::generate_block(VoxelBlockRequest &input) {
	ERR_FAIL_COND(input.voxel_buffer.is_null());
}

//bool VoxelGenerator::is_thread_safe() const {
//	return false;
//}

//bool VoxelGenerator::is_cloneable() const {
//	return false;
//}

void VoxelGenerator::emerge_block(Ref<VoxelBuffer> out_buffer, Vector3i origin_in_voxels, int lod) {
	VoxelBlockRequest r = { out_buffer, Vector3i(origin_in_voxels), lod };
	generate_block(r);
}

void VoxelGenerator::_b_generate_block(Ref<VoxelBuffer> out_buffer, Vector3 origin_in_voxels, int lod) {
	ERR_FAIL_COND(lod < 0);
	VoxelBlockRequest r = { out_buffer, Vector3i(origin_in_voxels), lod };
	generate_block(r);
}

void VoxelGenerator::_bind_methods() {
	ClassDB::bind_method(D_METHOD("generate_block", "out_buffer", "origin_in_voxels", "lod"),
			&VoxelGenerator::_b_generate_block);
}
