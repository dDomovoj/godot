#include "editor_vcs_interface.h"

EditorVCSInterface *EditorVCSInterface::singleton = NULL;

void EditorVCSInterface::_bind_methods() {

	// Proxy end points that act as fallbacks to unavailability of a function in the VCS addon
	ClassDB::bind_method(D_METHOD("_initialize", "project_root_path"), &EditorVCSInterface::_initialize);
	ClassDB::bind_method(D_METHOD("_get_is_vcs_intialized"), &EditorVCSInterface::_get_is_vcs_intialized);
	ClassDB::bind_method(D_METHOD("_get_vcs_name"), &EditorVCSInterface::_get_vcs_name);
	ClassDB::bind_method(D_METHOD("_shut_down"), &EditorVCSInterface::_shut_down);
	ClassDB::bind_method(D_METHOD("_get_project_name"), &EditorVCSInterface::_get_project_name);
	ClassDB::bind_method(D_METHOD("_get_modified_files_data"), &EditorVCSInterface::_get_modified_files_data);
	ClassDB::bind_method(D_METHOD("_commit", "msg"), &EditorVCSInterface::_commit);
	ClassDB::bind_method(D_METHOD("_get_file_diff", "file_path"), &EditorVCSInterface::_get_file_diff);
	ClassDB::bind_method(D_METHOD("_stage_file", "file_path"), &EditorVCSInterface::_stage_file);
	ClassDB::bind_method(D_METHOD("_unstage_file", "file_path"), &EditorVCSInterface::_unstage_file);

	ClassDB::bind_method(D_METHOD("is_addon_ready"), &EditorVCSInterface::is_addon_ready);

	// API methods that redirect calls to the proxy end points
	ClassDB::bind_method(D_METHOD("initialize", "project_root_path"), &EditorVCSInterface::initialize);
	ClassDB::bind_method(D_METHOD("get_is_vcs_intialized"), &EditorVCSInterface::get_is_vcs_intialized);
	ClassDB::bind_method(D_METHOD("get_modified_files_data"), &EditorVCSInterface::get_modified_files_data);
	ClassDB::bind_method(D_METHOD("stage_file", "file_path"), &EditorVCSInterface::stage_file);
	ClassDB::bind_method(D_METHOD("unstage_file", "file_path"), &EditorVCSInterface::unstage_file);
	ClassDB::bind_method(D_METHOD("commit", "msg"), &EditorVCSInterface::commit);
	ClassDB::bind_method(D_METHOD("get_file_diff", "file_path"), &EditorVCSInterface::get_file_diff);
	ClassDB::bind_method(D_METHOD("shut_down"), &EditorVCSInterface::shut_down);
	ClassDB::bind_method(D_METHOD("get_project_name"), &EditorVCSInterface::get_project_name);
	ClassDB::bind_method(D_METHOD("get_vcs_name"), &EditorVCSInterface::get_vcs_name);
}

bool EditorVCSInterface::_initialize(String p_project_root_path) {

	WARN_PRINT("Selected VCS addon does not implement an initialization function. This warning will be suppressed.")
	return true;
}

bool EditorVCSInterface::_get_is_vcs_intialized() {

	return false;
}

Dictionary EditorVCSInterface::_get_modified_files_data() {

	return Dictionary();
}

void EditorVCSInterface::_stage_file(String p_file_path) {

	return;
}

void EditorVCSInterface::_unstage_file(String p_file_path) {

	return;
}

void EditorVCSInterface::_commit(String p_msg) {

	return;
}

Array EditorVCSInterface::_get_file_diff(String p_file_path) {

	return Array();
}

bool EditorVCSInterface::_shut_down() {

	return false;
}

String EditorVCSInterface::_get_project_name() {

	return String();
}

String EditorVCSInterface::_get_vcs_name() {

	return "";
}

bool EditorVCSInterface::initialize(String p_project_root_path) {

	is_initialized = call("_initialize", p_project_root_path);
	return is_initialized;
}

bool EditorVCSInterface::get_is_vcs_intialized() {

	return call("_get_is_vcs_intialized");
}

Dictionary EditorVCSInterface::get_modified_files_data() {

	return call("_get_modified_files_data");
}

void EditorVCSInterface::stage_file(String p_file_path) {

	if (is_addon_ready()) {

		call("_stage_file", p_file_path);
	}
	return;
}

void EditorVCSInterface::unstage_file(String p_file_path) {

	if (is_addon_ready()) {

		call("_unstage_file", p_file_path);
	}
	return;
}

bool EditorVCSInterface::is_addon_ready() {

	return is_initialized;
}

void EditorVCSInterface::commit(String p_msg) {

	if (is_addon_ready()) {

		call("_commit", p_msg);
	}
	return;
}

Array EditorVCSInterface::get_file_diff(String p_file_path) {

	if (is_addon_ready()) {

		return call("_get_file_diff", p_file_path);
	}
	return Array();
}

bool EditorVCSInterface::shut_down() {

	return call("_shut_down");
}

String EditorVCSInterface::get_project_name() {

	return call("_get_project_name");
}

String EditorVCSInterface::get_vcs_name() {

	return call("_get_vcs_name");
}

EditorVCSInterface::EditorVCSInterface() {

	is_initialized = false;
}

EditorVCSInterface::~EditorVCSInterface() {
}

EditorVCSInterface *EditorVCSInterface::get_singleton() {

	return singleton;
}

void EditorVCSInterface::set_singleton(EditorVCSInterface *p_singleton) {

	singleton = p_singleton;
}
