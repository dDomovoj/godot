/*************************************************************************/
/*  filesystem_dock.cpp                                                  */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2018 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2018 Godot Engine contributors (cf. AUTHORS.md)    */
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

#include "filesystem_dock.h"

#include "core/io/resource_loader.h"
#include "core/os/dir_access.h"
#include "core/os/file_access.h"
#include "core/os/keyboard.h"
#include "core/os/os.h"
#include "core/project_settings.h"
#include "editor_node.h"
#include "editor_settings.h"
#include "scene/main/viewport.h"

Ref<Texture> FileSystemDock::_get_tree_item_icon(EditorFileSystemDirectory *p_dir, int p_idx) {
	Ref<Texture> file_icon;
	if (!p_dir->get_file_import_is_valid(p_idx)) {
		file_icon = get_icon("ImportFail", "EditorIcons");
	} else {
		String file_type = p_dir->get_file_type(p_idx);
		file_icon = (has_icon(file_type, "EditorIcons")) ? get_icon(file_type, "EditorIcons") : get_icon("File", "EditorIcons");
	}
	return file_icon;
}

bool FileSystemDock::_create_tree(TreeItem *p_parent, EditorFileSystemDirectory *p_dir, Vector<String> &uncollapsed_paths) {

	bool parent_should_expand = false;

	// Create a tree item for the subdirectory
	TreeItem *subdirectory_item = tree->create_item(p_parent);
	String dname = p_dir->get_name();
	if (dname == "")
		dname = "res://";

	subdirectory_item->set_text(0, dname);
	subdirectory_item->set_icon(0, get_icon("Folder", "EditorIcons"));
	subdirectory_item->set_selectable(0, true);
	String lpath = p_dir->get_path();
	subdirectory_item->set_metadata(0, lpath);
	if (path == lpath || ((display_mode_setting == DISPLAY_MODE_SETTING_SPLIT) && path.get_base_dir() == lpath)) {
		subdirectory_item->select(0);
	}

	if ((path.begins_with(lpath) && path != lpath)) {
		subdirectory_item->set_collapsed(false);
	} else {
		bool is_collapsed = true;
		for (int i = 0; i < uncollapsed_paths.size(); i++) {
			if (lpath == uncollapsed_paths[i]) {
				is_collapsed = false;
				break;
			}
		}
		subdirectory_item->set_collapsed(is_collapsed);
	}
	if (searched_string.length() > 0 && dname.to_lower().find(searched_string) >= 0) {
		parent_should_expand = true;
	}

	// Create items for all subdirectories
	for (int i = 0; i < p_dir->get_subdir_count(); i++)
		parent_should_expand = (_create_tree(subdirectory_item, p_dir->get_subdir(i), uncollapsed_paths) || parent_should_expand);

	// Create all items for the files in the subdirectory
	if (display_mode_setting == DISPLAY_MODE_SETTING_TREE_ONLY) {
		for (int i = 0; i < p_dir->get_file_count(); i++) {
			String file_name = p_dir->get_file(i);

			if (searched_string.length() > 0) {
				if (file_name.to_lower().find(searched_string) < 0) {
					// The seached string is not in the file name, we skip it
					continue;
				} else {
					// We expand all parents
					parent_should_expand = true;
				}
			}

			TreeItem *file_item = tree->create_item(subdirectory_item);
			file_item->set_text(0, file_name);
			file_item->set_icon(0, _get_tree_item_icon(p_dir, i));
			String file_metadata = lpath.plus_file(file_name);
			file_item->set_metadata(0, file_metadata);
			if (path == file_metadata) {
				file_item->select(0);
				file_item->set_as_cursor(0);
			}
			Array udata;
			udata.push_back(tree_update_id);
			udata.push_back(file_item);
			EditorResourcePreview::get_singleton()->queue_resource_preview(file_metadata, this, "_tree_thumbnail_done", udata);
		}
	}

	if (searched_string.length() > 0) {
		if (parent_should_expand) {
			subdirectory_item->set_collapsed(false);
		} else if (dname != "res://") {
			subdirectory_item->get_parent()->remove_child(subdirectory_item);
		}
	}

	return parent_should_expand;
}

Vector<String> FileSystemDock::_compute_uncollapsed_paths() {
	// Register currently collapsed paths
	Vector<String> uncollapsed_paths;
	TreeItem *root = tree->get_root();
	if (root) {
		TreeItem *resTree = root->get_children()->get_next();
		if (resTree) {
			Vector<TreeItem *> needs_check;
			needs_check.push_back(resTree);

			while (needs_check.size()) {
				if (!needs_check[0]->is_collapsed()) {
					uncollapsed_paths.push_back(needs_check[0]->get_metadata(0));
					TreeItem *child = needs_check[0]->get_children();
					while (child) {
						needs_check.push_back(child);
						child = child->get_next();
					}
				}
				needs_check.remove(0);
			}
		}
	}
	return uncollapsed_paths;
}

void FileSystemDock::_update_tree(const Vector<String> p_uncollapsed_paths, bool p_uncollapse_root) {

	// Recreate the tree
	tree->clear();
	tree_update_id++;
	updating_tree = true;
	TreeItem *root = tree->create_item();

	// Handles the favorites
	TreeItem *favorites = tree->create_item(root);
	favorites->set_icon(0, get_icon("Favorites", "EditorIcons"));
	favorites->set_text(0, TTR("Favorites:"));
	favorites->set_selectable(0, false);

	Vector<String> favorite_paths = EditorSettings::get_singleton()->get_favorite_dirs();
	for (int i = 0; i < favorite_paths.size(); i++) {
		String fave = favorite_paths[i];
		if (!fave.begins_with("res://"))
			continue;
		if (display_mode_setting == DISPLAY_MODE_SETTING_SPLIT && !fave.ends_with("/"))
			continue;

		Ref<Texture> folder_icon = get_icon("Folder", "EditorIcons");

		String text;
		Ref<Texture> icon;
		if (fave == "res://") {
			text = "/";
			icon = folder_icon;
		} else if (fave.ends_with("/")) {
			text = fave.substr(0, fave.length() - 1).get_file();
			icon = folder_icon;
		} else {
			text = fave.get_file();
			int index;
			EditorFileSystemDirectory *dir = EditorFileSystem::get_singleton()->find_file(fave, &index);
			if (dir) {
				icon = _get_tree_item_icon(dir, index);
			} else {
				icon = get_icon("File", "EditorIcons");
			}
		}

		if (searched_string.length() == 0 || text.to_lower().find(searched_string) >= 0) {
			TreeItem *ti = tree->create_item(favorites);
			ti->set_text(0, text);
			ti->set_icon(0, icon);
			ti->set_tooltip(0, fave);
			ti->set_selectable(0, true);
			ti->set_metadata(0, fave);
		}
	}

	Vector<String> uncollapsed_paths = p_uncollapsed_paths;
	if (p_uncollapse_root) {
		uncollapsed_paths.push_back("res://");
	}

	// Create the remaining of the tree
	_create_tree(root, EditorFileSystem::get_singleton()->get_filesystem(), uncollapsed_paths);
	tree->ensure_cursor_is_visible();
	updating_tree = false;
}

void FileSystemDock::_update_display_mode() {
	bool compact_mode = get_size().height < int(EditorSettings::get_singleton()->get("docks/filesystem/split_mode_minimum_height"));

	// Compute the new display mode
	DisplayMode new_display_mode;
	if ((display_mode_setting == DISPLAY_MODE_SETTING_TREE_ONLY) || compact_mode) {
		new_display_mode = file_list_view ? DISPLAY_MODE_FILE_LIST_ONLY : DISPLAY_MODE_TREE_ONLY;
	} else {
		new_display_mode = DISPLAY_MODE_SPLIT;
	}

	if (new_display_mode != display_mode || old_display_mode_setting != display_mode_setting) {
		display_mode = new_display_mode;
		old_display_mode_setting = display_mode_setting;
		button_toggle_display_mode->set_pressed(display_mode_setting == DISPLAY_MODE_SETTING_SPLIT ? true : false);
		switch (display_mode) {
			case DISPLAY_MODE_TREE_ONLY:
				tree->show();
				tree->set_v_size_flags(SIZE_EXPAND_FILL);
				if (display_mode_setting == DISPLAY_MODE_SETTING_TREE_ONLY) {
					tree_search_box->show();
				} else {
					tree_search_box->hide();
				}

				_update_tree(_compute_uncollapsed_paths());
				file_list_vb->hide();
				break;

			case DISPLAY_MODE_FILE_LIST_ONLY:
				tree->hide();
				tree_search_box->hide();
				button_tree->show();

				file_list_vb->show();
				_update_files(true);
				break;

			case DISPLAY_MODE_SPLIT:
				tree->show();
				tree->set_v_size_flags(SIZE_EXPAND_FILL);
				button_tree->hide();
				tree->ensure_cursor_is_visible();
				tree_search_box->hide();
				_update_tree(_compute_uncollapsed_paths());

				file_list_vb->show();
				_update_files(true);
				break;
		}
	}
}

void FileSystemDock::_notification(int p_what) {

	switch (p_what) {

		case NOTIFICATION_RESIZED: {
			_update_display_mode();
		} break;
		case NOTIFICATION_ENTER_TREE: {

			if (initialized)
				return;
			initialized = true;

			EditorFileSystem::get_singleton()->connect("filesystem_changed", this, "_fs_changed");
			EditorResourcePreview::get_singleton()->connect("preview_invalidated", this, "_preview_invalidated");

			String ei = "EditorIcons";
			button_reload->set_icon(get_icon("Reload", ei));
			button_toggle_display_mode->set_icon(get_icon("Panels2", ei));
			button_tree->set_icon(get_icon("Filesystem", ei));
			_update_file_list_display_mode_button();
			button_file_list_display_mode->connect("pressed", this, "_change_file_display");

			files->connect("item_activated", this, "_file_list_activate_file");
			button_hist_next->connect("pressed", this, "_fw_history");
			button_hist_prev->connect("pressed", this, "_bw_history");
			tree_search_box->set_right_icon(get_icon("Search", ei));
			tree_search_box->set_clear_button_enabled(true);
			file_list_search_box->set_right_icon(get_icon("Search", ei));
			file_list_search_box->set_clear_button_enabled(true);

			button_hist_next->set_icon(get_icon("Forward", ei));
			button_hist_prev->set_icon(get_icon("Back", ei));
			file_list_popup->connect("id_pressed", this, "_file_list_rmb_option");
			tree_popup->connect("id_pressed", this, "_tree_rmb_option");

			button_tree->connect("pressed", this, "_go_to_tree", varray(), CONNECT_DEFERRED);
			current_path->connect("text_entered", this, "navigate_to_path");

			_update_display_mode();

			always_show_folders = bool(EditorSettings::get_singleton()->get("docks/filesystem/always_show_folders"));

			if (EditorFileSystem::get_singleton()->is_scanning()) {
				_set_scanning_mode();
			} else {
				_update_tree(Vector<String>(), true);
			}

		} break;
		case NOTIFICATION_PROCESS: {
			if (EditorFileSystem::get_singleton()->is_scanning()) {
				scanning_progress->set_value(EditorFileSystem::get_singleton()->get_scanning_progress() * 100);
			}
		} break;
		case NOTIFICATION_EXIT_TREE: {

		} break;
		case NOTIFICATION_DRAG_BEGIN: {

			Dictionary dd = get_viewport()->gui_get_drag_data();
			if (tree->is_visible_in_tree() && dd.has("type")) {
				if ((String(dd["type"]) == "files") || (String(dd["type"]) == "files_and_dirs") || (String(dd["type"]) == "resource")) {
					tree->set_drop_mode_flags(Tree::DROP_MODE_ON_ITEM | Tree::DROP_MODE_INBETWEEN);
				} else if ((String(dd["type"]) == "favorite")) {
					tree->set_drop_mode_flags(Tree::DROP_MODE_INBETWEEN);
				}
			}

		} break;
		case NOTIFICATION_DRAG_END: {

			tree->set_drop_mode_flags(0);

		} break;
		case EditorSettings::NOTIFICATION_EDITOR_SETTINGS_CHANGED: {
			// Update icons
			String ei = "EditorIcons";
			button_reload->set_icon(get_icon("Reload", ei));
			button_toggle_display_mode->set_icon(get_icon("Panels2", ei));
			button_tree->set_icon(get_icon("Filesystem", ei));
			button_hist_next->set_icon(get_icon("Forward", ei));
			button_hist_prev->set_icon(get_icon("Back", ei));

			tree_search_box->set_right_icon(get_icon("Search", ei));
			tree_search_box->set_clear_button_enabled(true);
			file_list_search_box->set_right_icon(get_icon("Search", ei));
			file_list_search_box->set_clear_button_enabled(true);

			bool should_update_files = false;

			// Update file list display mode
			int new_file_list_mode = int(EditorSettings::get_singleton()->get("docks/filesystem/files_display_mode"));
			if (new_file_list_mode != file_list_display_mode) {
				set_file_list_display_mode(new_file_list_mode);
				_update_file_list_display_mode_button();
				should_update_files = true;
			}

			// Update display of files in tree
			display_mode_setting = DisplayModeSetting(int(EditorSettings::get_singleton()->get("docks/filesystem/display_mode")));

			// Update allways showfolders
			bool new_always_show_folders = bool(EditorSettings::get_singleton()->get("docks/filesystem/always_show_folders"));
			if (new_always_show_folders != always_show_folders) {
				always_show_folders = new_always_show_folders;
				should_update_files = true;
			}

			if (should_update_files) {
				_update_files(true);
			}

			// Change full tree mode
			_update_display_mode();

		} break;
	}
}

void FileSystemDock::_tree_multi_selected(Object *p_item, int p_column, bool p_selected) {

	// Return if we don't select something new
	if (!p_selected)
		return;

	// Tree item selected
	TreeItem *sel = tree->get_selected();
	if (!sel)
		return;
	path = sel->get_metadata(0);

	// Set the current path
	current_path->set_text(path);
	_push_to_history();

	// Update the file list
	if (!updating_tree && display_mode == DISPLAY_MODE_SPLIT) {
		_update_files(false);
	}
}

String FileSystemDock::get_selected_path() const {
	if (path.ends_with("/"))
		return path;
	else
		return path.get_base_dir();
}

String FileSystemDock::get_current_path() const {

	return path;
}

void FileSystemDock::navigate_to_path(const String &p_path) {

	String target_path = p_path;
	// If the path is a file, do not only go to the directory in the tree, also select the file in the file list.
	if (target_path.ends_with("/")) {
		target_path = target_path.substr(0, target_path.length() - 1);
	}
	DirAccess *dirAccess = DirAccess::open("res://");
	if (dirAccess->file_exists(p_path)) {
		path = target_path;
	} else if (dirAccess->dir_exists(p_path)) {
		path = target_path + "/";
	} else {
		ERR_EXPLAIN(vformat(TTR("Cannot navigate to '%s' as it has not been found in the file system!"), p_path));
		ERR_FAIL();
	}

	current_path->set_text(path);
	_push_to_history();

	if (display_mode == DISPLAY_MODE_SPLIT) {
		_update_tree(_compute_uncollapsed_paths());
		_update_files(false);
	} else if (display_mode == DISPLAY_MODE_TREE_ONLY) {
		if (path.ends_with("/")) {
			_go_to_file_list();
		} else {
			_update_tree(_compute_uncollapsed_paths());
		}
	} else { // DISPLAY_MODE_FILE_LIST_ONLY
		_update_files(true);
	}

	String file_name = p_path.get_file();
	if (!file_name.empty()) {
		for (int i = 0; i < files->get_item_count(); i++) {
			if (files->get_item_text(i) == file_name) {
				files->select(i, true);
				files->ensure_current_is_visible();
				break;
			}
		}
	}
}

void FileSystemDock::_file_list_thumbnail_done(const String &p_path, const Ref<Texture> &p_preview, const Ref<Texture> &p_small_preview, const Variant &p_udata) {

	if ((file_list_vb->is_visible_in_tree() || path == p_path.get_base_dir()) && p_preview.is_valid()) {
		Array uarr = p_udata;
		int idx = uarr[0];
		String file = uarr[1];
		if (idx < files->get_item_count() && files->get_item_text(idx) == file && files->get_item_metadata(idx) == p_path)
			files->set_item_icon(idx, p_preview);
	}
}

void FileSystemDock::_tree_thumbnail_done(const String &p_path, const Ref<Texture> &p_preview, const Ref<Texture> &p_small_preview, const Variant &p_udata) {
	if (p_small_preview.is_valid()) {
		Array uarr = p_udata;
		if (tree_update_id == (int)uarr[0]) {
			TreeItem *file_item = Object::cast_to<TreeItem>(uarr[1]);
			if (file_item) {
				file_item->set_icon(0, p_small_preview);

				// Update the favorite icon if needed
				TreeItem *favorite = tree->get_root()->get_children()->get_children();
				while (favorite) {
					if (favorite->get_metadata(0) == file_item->get_metadata(0)) {
						favorite->set_icon(0, p_small_preview);
					}
					favorite = favorite->get_next();
				}
			}
		}
	}
}

void FileSystemDock::_update_file_list_display_mode_button() {

	if (button_file_list_display_mode->is_pressed()) {
		file_list_display_mode = FILE_LIST_DISPLAY_LIST;
		button_file_list_display_mode->set_icon(get_icon("FileThumbnail", "EditorIcons"));
		button_file_list_display_mode->set_tooltip(TTR("View items as a grid of thumbnails."));
	} else {
		file_list_display_mode = FILE_LIST_DISPLAY_THUMBNAILS;
		button_file_list_display_mode->set_icon(get_icon("FileList", "EditorIcons"));
		button_file_list_display_mode->set_tooltip(TTR("View items as a list."));
	}
}

void FileSystemDock::_change_file_display() {

	_update_file_list_display_mode_button();

	EditorSettings::get_singleton()->set("docks/filesystem/files_display_mode", file_list_display_mode);

	_update_files(true);
}

void FileSystemDock::_search(EditorFileSystemDirectory *p_path, List<FileInfo> *matches, int p_max_items) {

	if (matches->size() > p_max_items)
		return;

	for (int i = 0; i < p_path->get_subdir_count(); i++) {
		_search(p_path->get_subdir(i), matches, p_max_items);
	}

	for (int i = 0; i < p_path->get_file_count(); i++) {
		String file = p_path->get_file(i);

		if (file.to_lower().find(searched_string) != -1) {

			FileInfo fi;
			fi.name = file;
			fi.type = p_path->get_file_type(i);
			fi.path = p_path->get_file_path(i);
			fi.import_broken = !p_path->get_file_import_is_valid(i);
			fi.import_status = 0;

			matches->push_back(fi);
			if (matches->size() > p_max_items)
				return;
		}
	}
}

void FileSystemDock::_update_files(bool p_keep_selection) {

	// Register the previously selected items
	Set<String> cselection;
	if (p_keep_selection) {
		for (int i = 0; i < files->get_item_count(); i++) {
			if (files->is_selected(i))
				cselection.insert(files->get_item_text(i));
		}
	}

	files->clear();

	current_path->set_text(path);

	String directory = path;
	if (directory.ends_with("/") && directory != "res://") {
		directory = directory.substr(0, directory.length() - 1);
	}
	String file = "";
	EditorFileSystemDirectory *efd = EditorFileSystem::get_singleton()->get_filesystem_path(directory);
	if (!efd) {
		directory = path.get_base_dir();
		file = path.get_file();
		efd = EditorFileSystem::get_singleton()->get_filesystem_path(directory);
	}
	if (!efd)
		return;

	String ei = "EditorIcons";
	int thumbnail_size = EditorSettings::get_singleton()->get("docks/filesystem/thumbnail_size");
	thumbnail_size *= EDSCALE;
	Ref<Texture> folder_thumbnail;
	Ref<Texture> file_thumbnail;
	Ref<Texture> file_thumbnail_broken;

	bool use_thumbnails = (file_list_display_mode == FILE_LIST_DISPLAY_THUMBNAILS);
	bool use_folders = searched_string.length() == 0 && ((display_mode == DISPLAY_MODE_FILE_LIST_ONLY || display_mode == DISPLAY_MODE_TREE_ONLY) || always_show_folders);

	if (use_thumbnails) {

		files->set_max_columns(0);
		files->set_icon_mode(ItemList::ICON_MODE_TOP);
		files->set_fixed_column_width(thumbnail_size * 3 / 2);
		files->set_max_text_lines(2);
		files->set_fixed_icon_size(Size2(thumbnail_size, thumbnail_size));

		if (thumbnail_size < 64) {
			folder_thumbnail = get_icon("FolderMediumThumb", ei);
			file_thumbnail = get_icon("FileMediumThumb", ei);
			file_thumbnail_broken = get_icon("FileDeadMediumThumb", ei);
		} else {
			folder_thumbnail = get_icon("FolderBigThumb", ei);
			file_thumbnail = get_icon("FileBigThumb", ei);
			file_thumbnail_broken = get_icon("FileDeadBigThumb", ei);
		}
	} else {

		files->set_icon_mode(ItemList::ICON_MODE_LEFT);
		files->set_max_columns(1);
		files->set_max_text_lines(1);
		files->set_fixed_column_width(0);
		files->set_fixed_icon_size(Size2());
	}

	if (use_folders) {
		Ref<Texture> folderIcon = (use_thumbnails) ? folder_thumbnail : get_icon("folder", "FileDialog");

		if (directory != "res://") {
			files->add_item("..", folderIcon, true);

			String bd = directory.get_base_dir();
			if (bd != "res://" && !bd.ends_with("/"))
				bd += "/";

			files->set_item_metadata(files->get_item_count() - 1, bd);
			files->set_item_selectable(files->get_item_count() - 1, false);
		}

		for (int i = 0; i < efd->get_subdir_count(); i++) {

			String dname = efd->get_subdir(i)->get_name();

			files->add_item(dname, folderIcon, true);
			files->set_item_metadata(files->get_item_count() - 1, directory.plus_file(dname) + "/");

			if (cselection.has(dname)) {
				files->select(files->get_item_count() - 1, false);
			}
		}
	}

	List<FileInfo> filelist;

	if (searched_string.length() > 0) {

		_search(EditorFileSystem::get_singleton()->get_filesystem(), &filelist, 128);
		filelist.sort();
	} else {

		for (int i = 0; i < efd->get_file_count(); i++) {

			FileInfo fi;
			fi.name = efd->get_file(i);
			fi.path = directory.plus_file(fi.name);
			fi.type = efd->get_file_type(i);
			fi.import_broken = !efd->get_file_import_is_valid(i);
			fi.import_status = 0;

			filelist.push_back(fi);
		}
		filelist.sort();
	}

	String oi = "Object";

	for (List<FileInfo>::Element *E = filelist.front(); E; E = E->next()) {
		FileInfo *finfo = &(E->get());
		String fname = finfo->name;
		String fpath = finfo->path;
		String ftype = finfo->type;

		Ref<Texture> type_icon;
		Ref<Texture> big_icon;

		String tooltip = fname;

		if (!finfo->import_broken) {
			type_icon = (has_icon(ftype, ei)) ? get_icon(ftype, ei) : get_icon(oi, ei);
			big_icon = file_thumbnail;
		} else {
			type_icon = get_icon("ImportFail", ei);
			big_icon = file_thumbnail_broken;
			tooltip += "\n" + TTR("Status: Import of file failed. Please fix file and reimport manually.");
		}

		int item_index;
		if (use_thumbnails) {
			files->add_item(fname, big_icon, true);
			item_index = files->get_item_count() - 1;
			files->set_item_metadata(item_index, fpath);
			files->set_item_tag_icon(item_index, type_icon);
			if (!finfo->import_broken) {
				Array udata;
				udata.resize(2);
				udata[0] = item_index;
				udata[1] = fname;
				EditorResourcePreview::get_singleton()->queue_resource_preview(fpath, this, "_file_list_thumbnail_done", udata);
			}
		} else {
			files->add_item(fname, type_icon, true);
			item_index = files->get_item_count() - 1;
			files->set_item_metadata(item_index, fpath);
		}

		if (cselection.has(fname))
			files->select(item_index, false);

		if (!p_keep_selection && file != "" && fname == file) {
			files->select(item_index, true);
			files->ensure_current_is_visible();
		}

		if (finfo->sources.size()) {
			for (int j = 0; j < finfo->sources.size(); j++) {
				tooltip += "\nSource: " + finfo->sources[j];
			}
		}
		files->set_item_tooltip(item_index, tooltip);
	}
}

void FileSystemDock::_select_file(const String p_path) {
	String fpath = p_path;
	if (fpath.ends_with("/")) {
		if (fpath != "res://") {
			fpath = fpath.substr(0, fpath.length() - 1);
		}
		navigate_to_path(fpath);
	} else {
		if (ResourceLoader::get_resource_type(fpath) == "PackedScene") {
			editor->open_request(fpath);
		} else {
			editor->load_resource(fpath);
		}
	}
}

void FileSystemDock::_tree_activate_file() {
	TreeItem *selected = tree->get_selected();
	if (selected) {
		call_deferred("_select_file", selected->get_metadata(0));
	}
}

void FileSystemDock::_file_list_activate_file(int p_idx) {
	_select_file(files->get_item_metadata(p_idx));
}

void FileSystemDock::_go_to_file_list() {

	if (display_mode == DISPLAY_MODE_TREE_ONLY) {

		file_list_view = true;
		_update_display_mode();
	} else {
		bool collapsed = tree->get_selected()->is_collapsed();
		tree->get_selected()->set_collapsed(!collapsed);
		_update_files(false);
	}
}

void FileSystemDock::_go_to_tree() {

	file_list_view = false;
	tree->grab_focus();
	_update_display_mode();
	tree->ensure_cursor_is_visible();
}

void FileSystemDock::_preview_invalidated(const String &p_path) {

	if (file_list_display_mode == FILE_LIST_DISPLAY_THUMBNAILS && p_path.get_base_dir() == path && searched_string.length() == 0 && file_list_vb->is_visible_in_tree()) {

		for (int i = 0; i < files->get_item_count(); i++) {

			if (files->get_item_metadata(i) == p_path) {
				//re-request preview
				Array udata;
				udata.resize(2);
				udata[0] = i;
				udata[1] = files->get_item_text(i);
				EditorResourcePreview::get_singleton()->queue_resource_preview(p_path, this, "_file_list_thumbnail_done", udata);
				break;
			}
		}
	}
}

void FileSystemDock::_fs_changed() {

	button_hist_prev->set_disabled(history_pos == 0);
	button_hist_next->set_disabled(history_pos == history.size() - 1);
	scanning_vb->hide();
	split_box->show();

	if (tree->is_visible()) {
		_update_tree(_compute_uncollapsed_paths());
	}

	if (file_list_vb->is_visible()) {
		_update_files(true);
	}

	set_process(false);
}

void FileSystemDock::_set_scanning_mode() {

	button_hist_prev->set_disabled(true);
	button_hist_next->set_disabled(true);
	split_box->hide();
	scanning_vb->show();
	set_process(true);
	if (EditorFileSystem::get_singleton()->is_scanning()) {
		scanning_progress->set_value(EditorFileSystem::get_singleton()->get_scanning_progress() * 100);
	} else {
		scanning_progress->set_value(0);
	}
}

void FileSystemDock::_fw_history() {

	if (history_pos < history.size() - 1)
		history_pos++;

	_update_history();
}

void FileSystemDock::_bw_history() {
	if (history_pos > 0)
		history_pos--;

	_update_history();
}

void FileSystemDock::_update_history() {
	path = history[history_pos];
	current_path->set_text(path);

	if (tree->is_visible()) {
		_update_tree(_compute_uncollapsed_paths());
		tree->grab_focus();
		tree->ensure_cursor_is_visible();
	}

	if (file_list_vb->is_visible()) {
		_update_files(false);
	}

	button_hist_prev->set_disabled(history_pos == 0);
	button_hist_next->set_disabled(history_pos == history.size() - 1);
}

void FileSystemDock::_push_to_history() {
	if (history[history_pos] != path) {
		history.resize(history_pos + 1);
		history.push_back(path);
		history_pos++;

		if (history.size() > history_max_size) {
			history.remove(0);
			history_pos = history_max_size - 1;
		}
	}

	button_hist_prev->set_disabled(history_pos == 0);
	button_hist_next->set_disabled(history_pos == history.size() - 1);
}

void FileSystemDock::_get_all_items_in_dir(EditorFileSystemDirectory *efsd, Vector<String> &files, Vector<String> &folders) const {
	if (efsd == NULL)
		return;

	for (int i = 0; i < efsd->get_subdir_count(); i++) {
		folders.push_back(efsd->get_subdir(i)->get_path());
		_get_all_items_in_dir(efsd->get_subdir(i), files, folders);
	}
	for (int i = 0; i < efsd->get_file_count(); i++) {
		files.push_back(efsd->get_file_path(i));
	}
}

void FileSystemDock::_find_remaps(EditorFileSystemDirectory *efsd, const Map<String, String> &renames, Vector<String> &to_remaps) const {
	for (int i = 0; i < efsd->get_subdir_count(); i++) {
		_find_remaps(efsd->get_subdir(i), renames, to_remaps);
	}
	for (int i = 0; i < efsd->get_file_count(); i++) {
		Vector<String> deps = efsd->get_file_deps(i);
		for (int j = 0; j < deps.size(); j++) {
			if (renames.has(deps[j])) {
				to_remaps.push_back(efsd->get_file_path(i));
				break;
			}
		}
	}
}

void FileSystemDock::_try_move_item(const FileOrFolder &p_item, const String &p_new_path,
		Map<String, String> &p_file_renames, Map<String, String> &p_folder_renames) const {
	//Ensure folder paths end with "/"
	String old_path = (p_item.is_file || p_item.path.ends_with("/")) ? p_item.path : (p_item.path + "/");
	String new_path = (p_item.is_file || p_new_path.ends_with("/")) ? p_new_path : (p_new_path + "/");

	if (new_path == old_path) {
		return;
	} else if (old_path == "res://") {
		EditorNode::get_singleton()->add_io_error(TTR("Cannot move/rename resources root."));
		return;
	} else if (!p_item.is_file && new_path.begins_with(old_path)) {
		//This check doesn't erroneously catch renaming to a longer name as folder paths always end with "/"
		EditorNode::get_singleton()->add_io_error(TTR("Cannot move a folder into itself.") + "\n" + old_path + "\n");
		return;
	}

	//Build a list of files which will have new paths as a result of this operation
	Vector<String> file_changed_paths;
	Vector<String> folder_changed_paths;
	if (p_item.is_file) {
		file_changed_paths.push_back(old_path);
	} else {
		folder_changed_paths.push_back(old_path);
		_get_all_items_in_dir(EditorFileSystem::get_singleton()->get_filesystem_path(old_path), file_changed_paths, folder_changed_paths);
	}

	DirAccess *da = DirAccess::create(DirAccess::ACCESS_RESOURCES);
	print_verbose("Moving " + old_path + " -> " + new_path);
	Error err = da->rename(old_path, new_path);
	if (err == OK) {
		//Move/Rename any corresponding import settings too
		if (p_item.is_file && FileAccess::exists(old_path + ".import")) {
			err = da->rename(old_path + ".import", new_path + ".import");
			if (err != OK) {
				EditorNode::get_singleton()->add_io_error(TTR("Error moving:") + "\n" + old_path + ".import\n");
			}
		}

		// update scene if it is open
		for (int i = 0; i < file_changed_paths.size(); ++i) {
			String new_item_path = p_item.is_file ? new_path : file_changed_paths[i].replace_first(old_path, new_path);
			if (ResourceLoader::get_resource_type(new_item_path) == "PackedScene" && editor->is_scene_open(file_changed_paths[i])) {
				EditorData *ed = &editor->get_editor_data();
				for (int j = 0; j < ed->get_edited_scene_count(); j++) {
					if (ed->get_scene_path(j) == file_changed_paths[i]) {
						ed->get_edited_scene_root(j)->set_filename(new_item_path);
						break;
					}
				}
			}
		}

		//Only treat as a changed dependency if it was successfully moved
		for (int i = 0; i < file_changed_paths.size(); ++i) {
			p_file_renames[file_changed_paths[i]] = file_changed_paths[i].replace_first(old_path, new_path);
			print_verbose("  Remap: " + file_changed_paths[i] + " -> " + p_file_renames[file_changed_paths[i]]);
		}
		for (int i = 0; i < folder_changed_paths.size(); ++i) {
			p_folder_renames[folder_changed_paths[i]] = folder_changed_paths[i].replace_first(old_path, new_path);
		}
	} else {
		EditorNode::get_singleton()->add_io_error(TTR("Error moving:") + "\n" + old_path + "\n");
	}
	memdelete(da);
}

void FileSystemDock::_try_duplicate_item(const FileOrFolder &p_item, const String &p_new_path) const {
	//Ensure folder paths end with "/"
	String old_path = (p_item.is_file || p_item.path.ends_with("/")) ? p_item.path : (p_item.path + "/");
	String new_path = (p_item.is_file || p_new_path.ends_with("/")) ? p_new_path : (p_new_path + "/");

	if (new_path == old_path) {
		return;
	} else if (old_path == "res://") {
		EditorNode::get_singleton()->add_io_error(TTR("Cannot move/rename resources root."));
		return;
	} else if (!p_item.is_file && new_path.begins_with(old_path)) {
		//This check doesn't erroneously catch renaming to a longer name as folder paths always end with "/"
		EditorNode::get_singleton()->add_io_error(TTR("Cannot move a folder into itself.") + "\n" + old_path + "\n");
		return;
	}

	DirAccess *da = DirAccess::create(DirAccess::ACCESS_RESOURCES);
	print_verbose("Duplicating " + old_path + " -> " + new_path);
	Error err = p_item.is_file ? da->copy(old_path, new_path) : da->copy_dir(old_path, new_path);
	if (err == OK) {
		//Move/Rename any corresponding import settings too
		if (p_item.is_file && FileAccess::exists(old_path + ".import")) {
			err = da->copy(old_path + ".import", new_path + ".import");
			if (err != OK) {
				EditorNode::get_singleton()->add_io_error(TTR("Error duplicating:") + "\n" + old_path + ".import\n");
			}
		}
	} else {
		EditorNode::get_singleton()->add_io_error(TTR("Error duplicating:") + "\n" + old_path + "\n");
	}
	memdelete(da);
}

void FileSystemDock::_update_resource_paths_after_move(const Map<String, String> &p_renames) const {

	//Rename all resources loaded, be it subresources or actual resources
	List<Ref<Resource> > cached;
	ResourceCache::get_cached_resources(&cached);

	for (List<Ref<Resource> >::Element *E = cached.front(); E; E = E->next()) {

		Ref<Resource> r = E->get();

		String base_path = r->get_path();
		String extra_path;
		int sep_pos = r->get_path().find("::");
		if (sep_pos >= 0) {
			extra_path = base_path.substr(sep_pos, base_path.length());
			base_path = base_path.substr(0, sep_pos);
		}

		if (p_renames.has(base_path)) {
			base_path = p_renames[base_path];
		}

		r->set_path(base_path + extra_path);
	}

	for (int i = 0; i < EditorNode::get_editor_data().get_edited_scene_count(); i++) {

		String path;
		if (i == EditorNode::get_editor_data().get_edited_scene()) {
			if (!get_tree()->get_edited_scene_root())
				continue;

			path = get_tree()->get_edited_scene_root()->get_filename();
		} else {

			path = EditorNode::get_editor_data().get_scene_path(i);
		}

		if (p_renames.has(path)) {
			path = p_renames[path];
		}

		if (i == EditorNode::get_editor_data().get_edited_scene()) {

			get_tree()->get_edited_scene_root()->set_filename(path);
		} else {

			EditorNode::get_editor_data().set_scene_path(i, path);
		}
	}
}

void FileSystemDock::_update_dependencies_after_move(const Map<String, String> &p_renames) const {
	//The following code assumes that the following holds:
	// 1) EditorFileSystem contains the old paths/folder structure from before the rename/move.
	// 2) ResourceLoader can use the new paths without needing to call rescan.
	Vector<String> remaps;
	_find_remaps(EditorFileSystem::get_singleton()->get_filesystem(), p_renames, remaps);
	for (int i = 0; i < remaps.size(); ++i) {
		//Because we haven't called a rescan yet the found remap might still be an old path itself.
		String file = p_renames.has(remaps[i]) ? p_renames[remaps[i]] : remaps[i];
		print_verbose("Remapping dependencies for: " + file);
		Error err = ResourceLoader::rename_dependencies(file, p_renames);
		if (err == OK) {
			if (ResourceLoader::get_resource_type(file) == "PackedScene")
				editor->reload_scene(file);
		} else {
			EditorNode::get_singleton()->add_io_error(TTR("Unable to update dependencies:") + "\n" + remaps[i] + "\n");
		}
	}
}

void FileSystemDock::_update_project_settings_after_move(const Map<String, String> &p_renames) const {

	// Find all project settings of type FILE and replace them if needed
	const Map<StringName, PropertyInfo> prop_info = ProjectSettings::get_singleton()->get_custom_property_info();
	for (const Map<StringName, PropertyInfo>::Element *E = prop_info.front(); E; E = E->next()) {
		if (E->get().hint == PROPERTY_HINT_FILE) {
			String old_path = GLOBAL_GET(E->key());
			if (p_renames.has(old_path)) {
				ProjectSettings::get_singleton()->set_setting(E->key(), p_renames[old_path]);
			}
		};
	}
	ProjectSettings::get_singleton()->save();
}

void FileSystemDock::_update_favorite_dirs_list_after_move(const Map<String, String> &p_renames) const {

	Vector<String> favorite_dirs = EditorSettings::get_singleton()->get_favorite_dirs();
	Vector<String> new_favorite_dirs;

	for (int i = 0; i < favorite_dirs.size(); i++) {
		String old_path = favorite_dirs[i] + "/";

		if (p_renames.has(old_path)) {
			String new_path = p_renames[old_path];
			new_favorite_dirs.push_back(new_path.substr(0, new_path.length() - 1));
		} else {
			new_favorite_dirs.push_back(favorite_dirs[i]);
		}
	}

	EditorSettings::get_singleton()->set_favorite_dirs(new_favorite_dirs);
}

void FileSystemDock::_make_dir_confirm() {
	String dir_name = make_dir_dialog_text->get_text().strip_edges();

	if (dir_name.length() == 0) {
		EditorNode::get_singleton()->show_warning(TTR("No name provided"));
		return;
	} else if (dir_name.find("/") != -1 || dir_name.find("\\") != -1 || dir_name.find(":") != -1 || dir_name.ends_with(".") || dir_name.ends_with(" ")) {
		EditorNode::get_singleton()->show_warning(TTR("Provided name contains invalid characters"));
		return;
	}

	String directory = path;
	if (!directory.ends_with("/")) {
		directory = directory.get_base_dir();
	}
	print_verbose("Making folder " + dir_name + " in " + directory);
	DirAccess *da = DirAccess::create(DirAccess::ACCESS_RESOURCES);
	Error err = da->change_dir(directory);
	if (err == OK) {
		err = da->make_dir(dir_name);
	}
	memdelete(da);

	if (err == OK) {
		print_verbose("FileSystem: calling rescan.");
		_rescan();
	} else {
		EditorNode::get_singleton()->show_warning(TTR("Could not create folder."));
	}
}

void FileSystemDock::_rename_operation_confirm() {

	String new_name = rename_dialog_text->get_text().strip_edges();
	if (new_name.length() == 0) {
		EditorNode::get_singleton()->show_warning(TTR("No name provided."));
		return;
	} else if (new_name.find("/") != -1 || new_name.find("\\") != -1 || new_name.find(":") != -1) {
		EditorNode::get_singleton()->show_warning(TTR("Name contains invalid characters."));
		return;
	}

	String old_path = to_rename.path.ends_with("/") ? to_rename.path.substr(0, to_rename.path.length() - 1) : to_rename.path;
	String new_path = old_path.get_base_dir().plus_file(new_name);
	if (old_path == new_path) {
		return;
	}

	//Present a more user friendly warning for name conflict
	DirAccess *da = DirAccess::create(DirAccess::ACCESS_RESOURCES);
#if defined(WINDOWS_ENABLED) || defined(UWP_ENABLED)
	// Workaround case insensitivity on Windows
	if ((da->file_exists(new_path) || da->dir_exists(new_path)) && new_path.to_lower() != old_path.to_lower()) {
#else
	if (da->file_exists(new_path) || da->dir_exists(new_path)) {
#endif
		EditorNode::get_singleton()->show_warning(TTR("A file or folder with this name already exists."));
		memdelete(da);
		return;
	}
	memdelete(da);

	Map<String, String> file_renames;
	Map<String, String> folder_renames;
	_try_move_item(to_rename, new_path, file_renames, folder_renames);
	_update_dependencies_after_move(file_renames);
	_update_resource_paths_after_move(file_renames);
	_update_project_settings_after_move(file_renames);
	_update_favorite_dirs_list_after_move(folder_renames);

	//Rescan everything
	print_verbose("FileSystem: calling rescan.");
	_rescan();
}

void FileSystemDock::_duplicate_operation_confirm() {

	String new_name = duplicate_dialog_text->get_text().strip_edges();
	if (new_name.length() == 0) {
		EditorNode::get_singleton()->show_warning(TTR("No name provided."));
		return;
	} else if (new_name.find("/") != -1 || new_name.find("\\") != -1 || new_name.find(":") != -1) {
		EditorNode::get_singleton()->show_warning(TTR("Name contains invalid characters."));
		return;
	}

	String new_path;
	String base_dir = to_duplicate.path.get_base_dir();
	if (to_duplicate.is_file) {
		new_path = base_dir.plus_file(new_name);
	} else {
		new_path = base_dir.substr(0, base_dir.find_last("/")) + "/" + new_name;
	}

	//Present a more user friendly warning for name conflict
	DirAccess *da = DirAccess::create(DirAccess::ACCESS_RESOURCES);
	if (da->file_exists(new_path) || da->dir_exists(new_path)) {
		EditorNode::get_singleton()->show_warning(TTR("A file or folder with this name already exists."));
		memdelete(da);
		return;
	}
	memdelete(da);

	_try_duplicate_item(to_duplicate, new_path);

	//Rescan everything
	print_verbose("FileSystem: calling rescan.");
	_rescan();
}

void FileSystemDock::_move_with_overwrite() {
	_move_operation_confirm(to_move_path, true);
}

bool FileSystemDock::_check_existing() {
	String &p_to_path = to_move_path;
	for (int i = 0; i < to_move.size(); i++) {
		String ol_pth = to_move[i].path.ends_with("/") ? to_move[i].path.substr(0, to_move[i].path.length() - 1) : to_move[i].path;
		String p_new_path = p_to_path.plus_file(ol_pth.get_file());
		FileOrFolder p_item = to_move[i];

		String old_path = (p_item.is_file || p_item.path.ends_with("/")) ? p_item.path : (p_item.path + "/");
		String new_path = (p_item.is_file || p_new_path.ends_with("/")) ? p_new_path : (p_new_path + "/");

		if (p_item.is_file && FileAccess::exists(new_path)) {
			return false;
		} else if (!p_item.is_file && DirAccess::exists(new_path)) {
			return false;
		}
	}
	return true;
}

void FileSystemDock::_move_operation_confirm(const String &p_to_path, bool overwrite) {
	if (!overwrite) {
		to_move_path = p_to_path;
		bool can_move = _check_existing();
		if (!can_move) {
			//ask to do something
			overwrite_dialog->popup_centered_minsize();
			overwrite_dialog->grab_focus();
			return;
		}
	}

	Map<String, String> file_renames;
	Map<String, String> folder_renames;
	bool is_moved = false;
	for (int i = 0; i < to_move.size(); i++) {
		String old_path = to_move[i].path.ends_with("/") ? to_move[i].path.substr(0, to_move[i].path.length() - 1) : to_move[i].path;
		String new_path = p_to_path.plus_file(old_path.get_file());
		if (old_path != new_path) {
			_try_move_item(to_move[i], new_path, file_renames, folder_renames);
			is_moved = true;
		}
	}

	if (is_moved) {
		_update_dependencies_after_move(file_renames);
		_update_resource_paths_after_move(file_renames);
		_update_project_settings_after_move(file_renames);
		_update_favorite_dirs_list_after_move(folder_renames);

		print_verbose("FileSystem: calling rescan.");
		_rescan();
	}
}

Vector<String> FileSystemDock::_tree_get_selected(bool remove_self_inclusion) {
	// Build a list of selected items with the active one at the first position
	Vector<String> selected_strings;

	TreeItem *active_selected = tree->get_selected();
	if (active_selected) {
		selected_strings.push_back(active_selected->get_metadata(0));
	}

	TreeItem *selected = tree->get_root();
	selected = tree->get_next_selected(selected);
	while (selected) {
		if (selected != active_selected) {
			selected_strings.push_back(selected->get_metadata(0));
		}
		selected = tree->get_next_selected(selected);
	}

	// Remove paths or files that are included into another
	if (remove_self_inclusion && selected_strings.size() > 1) {
		selected_strings.sort_custom<NaturalNoCaseComparator>();
		String last_path = "";
		for (int i = 0; i < selected_strings.size(); i++) {
			if (last_path != "" && selected_strings[i].begins_with(last_path)) {
				selected_strings.remove(i);
				i--;
			}
			if (selected_strings[i].ends_with("/")) {
				last_path = selected_strings[i];
			}
		}
	}
	return selected_strings;
}

void FileSystemDock::_tree_rmb_option(int p_option) {

	Vector<String> selected_strings = _tree_get_selected();

	// Execute the current option
	switch (p_option) {
		case FOLDER_EXPAND_ALL:
		case FOLDER_COLLAPSE_ALL: {
			// Expand or collapse the folder
			if (selected_strings.size() == 1) {
				bool is_collapsed = (p_option == FOLDER_COLLAPSE_ALL);

				Vector<TreeItem *> needs_check;
				needs_check.push_back(tree->get_selected());

				while (needs_check.size()) {
					needs_check[0]->set_collapsed(is_collapsed);

					TreeItem *child = needs_check[0]->get_children();
					while (child) {
						needs_check.push_back(child);
						child = child->get_next();
					}

					needs_check.remove(0);
				}
			}
		} break;
		default: {
			_file_option(p_option, selected_strings);
		} break;
	}
}

void FileSystemDock::_file_list_rmb_option(int p_option) {
	Vector<int> selected_id = files->get_selected_items();
	Vector<String> selected;
	for (int i = 0; i < selected_id.size(); i++) {
		selected.push_back(files->get_item_metadata(selected_id[i]));
	}
	_file_option(p_option, selected);
}

void FileSystemDock::_file_option(int p_option, const Vector<String> p_selected) {
	// The first one should be the active item

	switch (p_option) {
		case FILE_SHOW_IN_EXPLORER: {
			// Show the file / folder in the OS explorer
			String fpath = path;
			if (!fpath.ends_with("/")) {
				fpath = fpath.get_base_dir();
			}
			String dir = ProjectSettings::get_singleton()->globalize_path(fpath);
			OS::get_singleton()->shell_open(String("file://") + dir);
		} break;

		case FILE_OPEN: {
			// Open the file
			for (int i = 0; i < p_selected.size(); i++) {
				_select_file(p_selected[i]);
			}
		} break;

		case FILE_INSTANCE: {
			// Instance all selected scenes
			Vector<String> paths;
			for (int i = 0; i < p_selected.size(); i++) {
				String fpath = p_selected[i];
				if (EditorFileSystem::get_singleton()->get_file_type(fpath) == "PackedScene") {
					paths.push_back(fpath);
				}
			}
			if (!paths.empty()) {
				emit_signal("instance", paths);
			}
		} break;

		case FILE_ADD_FAVORITE: {
			// Add the files from favorites
			Vector<String> favorites = EditorSettings::get_singleton()->get_favorite_dirs();
			for (int i = 0; i < p_selected.size(); i++) {
				if (favorites.find(p_selected[i]) == -1) {
					favorites.push_back(p_selected[i]);
				}
			}
			EditorSettings::get_singleton()->set_favorite_dirs(favorites);
			_update_tree(_compute_uncollapsed_paths());
		} break;

		case FILE_REMOVE_FAVORITE: {
			// Remove the files from favorites
			Vector<String> favorites = EditorSettings::get_singleton()->get_favorite_dirs();
			for (int i = 0; i < p_selected.size(); i++) {
				favorites.erase(p_selected[i]);
			}
			EditorSettings::get_singleton()->set_favorite_dirs(favorites);
			_update_tree(_compute_uncollapsed_paths());
		} break;

		case FILE_DEPENDENCIES: {
			// Checkout the file dependencies
			if (!p_selected.empty()) {
				String fpath = p_selected[0];
				deps_editor->edit(fpath);
			}
		} break;

		case FILE_OWNERS: {
			// Checkout the file owners
			if (!p_selected.empty()) {
				String fpath = p_selected[0];
				owners_editor->show(fpath);
			}
		} break;

		case FILE_MOVE: {
			// Move the files to a given location
			to_move.clear();
			for (int i = 0; i < p_selected.size(); i++) {
				String fpath = p_selected[i];
				if (fpath != "res://") {
					to_move.push_back(FileOrFolder(fpath, !fpath.ends_with("/")));
				}
			}
			if (to_move.size() > 0) {
				move_dialog->popup_centered_ratio();
			}
		} break;

		case FILE_RENAME: {
			// Rename the active file
			if (!p_selected.empty()) {
				to_rename.path = p_selected[0];
				if (to_rename.path != "res://") {
					to_rename.is_file = !to_rename.path.ends_with("/");
					if (to_rename.is_file) {
						String name = to_rename.path.get_file();
						rename_dialog->set_title(TTR("Renaming file:") + " " + name);
						rename_dialog_text->set_text(name);
						rename_dialog_text->select(0, name.find_last("."));
					} else {
						String name = to_rename.path.substr(0, to_rename.path.length() - 1).get_file();
						rename_dialog->set_title(TTR("Renaming folder:") + " " + name);
						rename_dialog_text->set_text(name);
						rename_dialog_text->select(0, name.length());
					}
					rename_dialog->popup_centered_minsize(Size2(250, 80) * EDSCALE);
					rename_dialog_text->grab_focus();
				}
			}
		} break;

		case FILE_REMOVE: {
			// Remove the selected files
			Vector<String> remove_files;
			Vector<String> remove_folders;

			for (int i = 0; i < p_selected.size(); i++) {
				String fpath = p_selected[i];
				if (fpath != "res://") {
					if (fpath.ends_with("/")) {
						remove_folders.push_back(fpath);
					} else {
						remove_files.push_back(fpath);
					}
				}
			}

			if (remove_files.size() + remove_folders.size() > 0) {
				remove_dialog->show(remove_folders, remove_files);
			}
		} break;

		case FILE_DUPLICATE: {
			// Duplicate the selected files
			for (int i = 0; i < p_selected.size(); i++) {
				to_duplicate.path = p_selected[i];
				to_duplicate.is_file = !to_duplicate.path.ends_with("/");
				if (to_duplicate.is_file) {
					String name = to_duplicate.path.get_file();
					duplicate_dialog->set_title(TTR("Duplicating file:") + " " + name);
					duplicate_dialog_text->set_text(name);
					duplicate_dialog_text->select(0, name.find_last("."));
				} else {
					String name = to_duplicate.path.substr(0, to_duplicate.path.length() - 1).get_file();
					duplicate_dialog->set_title(TTR("Duplicating folder:") + " " + name);
					duplicate_dialog_text->set_text(name);
					duplicate_dialog_text->select(0, name.length());
				}
				duplicate_dialog->popup_centered_minsize(Size2(250, 80) * EDSCALE);
				duplicate_dialog_text->grab_focus();
			}
		} break;

		case FILE_INFO: {

		} break;

		case FILE_REIMPORT: {
			// Reimport all selected files
			Vector<String> reimport;
			for (int i = 0; i < p_selected.size(); i++) {
				reimport.push_back(p_selected[i]);
			}

			ERR_FAIL_COND(reimport.size() == 0);
		} break;

		case FILE_NEW_FOLDER: {
			// Create a new folder
			make_dir_dialog_text->set_text("new folder");
			make_dir_dialog_text->select_all();
			make_dir_dialog->popup_centered_minsize(Size2(250, 80) * EDSCALE);
			make_dir_dialog_text->grab_focus();
		} break;

		case FILE_NEW_SCRIPT: {
			// Create a new script
			String fpath = path;
			if (!fpath.ends_with("/")) {
				fpath = fpath.get_base_dir();
			}
			make_script_dialog_text->config("Node", fpath + "new_script.gd");
			make_script_dialog_text->popup_centered(Size2(300, 300) * EDSCALE);
		} break;

		case FILE_COPY_PATH: {
			// Copy the file path
			if (!p_selected.empty()) {
				String fpath = p_selected[0];
				OS::get_singleton()->set_clipboard(fpath);
			}
		} break;

		case FILE_NEW_RESOURCE: {
			// Create a new resource
			new_resource_dialog->popup_create(true);
		} break;
	}
}

void FileSystemDock::_resource_created() const {
	Object *c = new_resource_dialog->instance_selected();

	ERR_FAIL_COND(!c);
	Resource *r = Object::cast_to<Resource>(c);
	ERR_FAIL_COND(!r);

	REF res(r);
	editor->push_item(c);

	RES current_res = RES(r);

	String fpath = path;
	if (!fpath.ends_with("/")) {
		fpath = fpath.get_base_dir();
	}

	editor->save_resource_as(current_res, fpath);
}

void FileSystemDock::_search_changed(const String &p_text, const Control *p_from) {
	if (searched_string.length() == 0 && p_text.length() > 0) {
		// Register the uncollapsed paths before they change
		uncollapsed_paths_before_search = _compute_uncollapsed_paths();
	}

	searched_string = p_text.to_lower();

	if (p_from == tree_search_box)
		file_list_search_box->set_text(searched_string);
	else // file_list_search_box
		tree_search_box->set_text(searched_string);

	switch (display_mode) {
		case DISPLAY_MODE_FILE_LIST_ONLY: {
			_update_files(false);
		} break;
		case DISPLAY_MODE_TREE_ONLY: {
			_update_tree(searched_string.length() == 0 ? uncollapsed_paths_before_search : Vector<String>());
		} break;
		case DISPLAY_MODE_SPLIT: {
			_update_files(false);
			_update_tree(searched_string.length() == 0 ? uncollapsed_paths_before_search : Vector<String>());
		} break;
	}
}

void FileSystemDock::_rescan() {

	_set_scanning_mode();
	EditorFileSystem::get_singleton()->scan();
}

void FileSystemDock::_toggle_split_mode(bool p_active) {
	display_mode_setting = p_active ? DISPLAY_MODE_SETTING_SPLIT : DISPLAY_MODE_SETTING_TREE_ONLY;
	EditorSettings::get_singleton()->set("docks/filesystem/display_mode", int(display_mode_setting));
	_update_display_mode();
}

void FileSystemDock::fix_dependencies(const String &p_for_file) {
	deps_editor->edit(p_for_file);
}

void FileSystemDock::focus_on_filter() {

	if (display_mode == DISPLAY_MODE_FILE_LIST_ONLY && tree->is_visible()) {
		// Tree mode, switch to files list with search box
		tree->hide();
		file_list_vb->show();
	}

	file_list_search_box->grab_focus();
}

void FileSystemDock::set_file_list_display_mode(int p_mode) {

	if (p_mode == file_list_display_mode)
		return;

	button_file_list_display_mode->set_pressed(p_mode == FILE_LIST_DISPLAY_LIST);
	_change_file_display();
}

Variant FileSystemDock::get_drag_data_fw(const Point2 &p_point, Control *p_from) {
	bool all_favorites = true;
	bool all_not_favorites = true;

	Vector<String> paths;

	if (p_from == tree) {
		// Check if the first selected is in favorite
		TreeItem *selected = tree->get_next_selected(tree->get_root());
		while (selected) {
			bool is_favorite = selected->get_parent() != NULL && tree->get_root()->get_children() == selected->get_parent();
			all_favorites &= is_favorite;
			all_not_favorites &= !is_favorite;
			selected = tree->get_next_selected(selected);
		}
		if (all_favorites) {
			paths = _tree_get_selected(false);
		} else {
			paths = _tree_get_selected();
		}
	} else if (p_from == files) {
		for (int i = 0; i < files->get_item_count(); i++) {
			if (files->is_selected(i)) {
				paths.push_back(files->get_item_metadata(i));
			}
		}
		all_favorites = false;
		all_not_favorites = true;
	}

	if (paths.empty())
		return Variant();

	if (!all_favorites && !all_not_favorites)
		return Variant();

	Dictionary drag_data = EditorNode::get_singleton()->drag_files_and_dirs(paths, p_from);
	if (all_favorites) {
		drag_data["type"] = "favorite";
	}
	return drag_data;
}

bool FileSystemDock::can_drop_data_fw(const Point2 &p_point, const Variant &p_data, Control *p_from) const {

	Dictionary drag_data = p_data;

	if (drag_data.has("type") && String(drag_data["type"]) == "favorite") {

		// Moving favorite around
		TreeItem *ti = tree->get_item_at_position(p_point);
		if (!ti)
			return false;

		int drop_section = tree->get_drop_section_at_position(p_point);
		TreeItem *favorites_item = tree->get_root()->get_children();

		TreeItem *resources_item = favorites_item->get_next();

		if (ti == favorites_item) {
			return (drop_section == 1); // The parent, first fav
		}
		if (ti->get_parent() && favorites_item == ti->get_parent()) {
			return true; // A favorite
		}
		if (ti == resources_item) {
			return (drop_section == -1); // The tree, last fav
		}

		return false;
	}

	if (drag_data.has("type") && String(drag_data["type"]) == "resource") {
		// Move resources
		String to_dir;
		bool favorite;
		_get_drag_target_folder(to_dir, favorite, p_point, p_from);
		return !to_dir.empty();
	}

	if (drag_data.has("type") && (String(drag_data["type"]) == "files" || String(drag_data["type"]) == "files_and_dirs")) {
		// Move files or dir
		String to_dir;
		bool favorite;
		_get_drag_target_folder(to_dir, favorite, p_point, p_from);

		if (favorite)
			return true;

		if (to_dir.empty())
			return false;

		//Attempting to move a folder into itself will fail later
		//Rather than bring up a message don't try to do it in the first place
		to_dir = to_dir.ends_with("/") ? to_dir : (to_dir + "/");
		Vector<String> fnames = drag_data["files"];
		for (int i = 0; i < fnames.size(); ++i) {
			if (fnames[i].ends_with("/") && to_dir.begins_with(fnames[i]))
				return false;
		}

		return true;
	}

	return false;
}

void FileSystemDock::drop_data_fw(const Point2 &p_point, const Variant &p_data, Control *p_from) {

	if (!can_drop_data_fw(p_point, p_data, p_from))
		return;
	Dictionary drag_data = p_data;

	Vector<String> dirs = EditorSettings::get_singleton()->get_favorite_dirs();

	if (drag_data.has("type") && String(drag_data["type"]) == "favorite") {
		// Moving favorite around
		TreeItem *ti = tree->get_item_at_position(p_point);
		if (!ti)
			return;
		int drop_section = tree->get_drop_section_at_position(p_point);

		int drop_position;
		Vector<String> files = drag_data["files"];
		TreeItem *favorites_item = tree->get_root()->get_children();
		TreeItem *resources_item = favorites_item->get_next();

		if (ti == favorites_item) {
			// Drop on the favorite folder
			drop_position = 0;
		} else if (ti == resources_item) {
			// Drop on the resouce item
			drop_position = dirs.size();
		} else {
			// Drop in the list
			drop_position = dirs.find(ti->get_metadata(0));
			if (drop_section == 1) {
				drop_position++;
			}
		}

		// Remove dragged favorites
		Vector<int> to_remove;
		int offset = 0;
		for (int i = 0; i < files.size(); i++) {
			int to_remove_pos = dirs.find(files[i]);
			to_remove.push_back(to_remove_pos);
			if (to_remove_pos < drop_position) {
				offset++;
			}
		}
		drop_position -= offset;
		to_remove.sort();
		for (int i = 0; i < to_remove.size(); i++) {
			dirs.remove(to_remove[i] - i);
		}

		// Re-add them at the right position
		for (int i = 0; i < files.size(); i++) {
			dirs.insert(drop_position, files[i]);
			drop_position++;
		}

		EditorSettings::get_singleton()->set_favorite_dirs(dirs);
		_update_tree(_compute_uncollapsed_paths());
		return;
	}

	if (drag_data.has("type") && String(drag_data["type"]) == "resource") {
		// Moving resource
		Ref<Resource> res = drag_data["resource"];
		String to_dir;
		bool favorite;
		_get_drag_target_folder(to_dir, favorite, p_point, p_from);
		if (res.is_valid() && !to_dir.empty()) {
			EditorNode::get_singleton()->push_item(res.ptr());
			EditorNode::get_singleton()->save_resource_as(res, to_dir);
		}
	}

	if (drag_data.has("type") && (String(drag_data["type"]) == "files" || String(drag_data["type"]) == "files_and_dirs")) {
		// Move files or add to favorites
		String to_dir;
		bool favorite;
		_get_drag_target_folder(to_dir, favorite, p_point, p_from);
		if (!to_dir.empty()) {
			Vector<String> fnames = drag_data["files"];
			to_move.clear();
			for (int i = 0; i < fnames.size(); i++) {
				to_move.push_back(FileOrFolder(fnames[i], !fnames[i].ends_with("/")));
			}
			_move_operation_confirm(to_dir);
		} else if (favorite) {
			// Add the files from favorites
			Vector<String> fnames = drag_data["files"];
			Vector<String> favorites = EditorSettings::get_singleton()->get_favorite_dirs();
			for (int i = 0; i < fnames.size(); i++) {
				if (favorites.find(fnames[i]) == -1) {
					favorites.push_back(fnames[i]);
				}
			}
			EditorSettings::get_singleton()->set_favorite_dirs(favorites);
			_update_tree(_compute_uncollapsed_paths());
		}
	}
}

void FileSystemDock::_get_drag_target_folder(String &target, bool &target_favorites, const Point2 &p_point, Control *p_from) const {
	target = String();
	target_favorites = false;

	// In the file list
	if (p_from == files) {
		int pos = files->get_item_at_position(p_point, true);
		if (pos == -1) {
			return;
		}

		String ltarget = files->get_item_metadata(pos);
		target = ltarget.ends_with("/") ? target : path;
		return;
	}

	// In the tree
	if (p_from == tree) {
		TreeItem *ti = tree->get_item_at_position(p_point);
		int section = tree->get_drop_section_at_position(p_point);
		if (ti) {
			// Check the favorites first
			if (ti == tree->get_root()->get_children() && section >= 0) {
				target_favorites = true;
				return;
			} else if (ti->get_parent() == tree->get_root()->get_children()) {
				target_favorites = true;
				return;
			} else {
				String fpath = ti->get_metadata(0);
				if (section == 0) {
					if (fpath.ends_with("/")) {
						// We drop on a folder
						target = fpath;
						return;
					}
				} else {
					if (ti->get_parent() != tree->get_root()->get_children()) {
						// Not in the favorite section
						if (fpath != "res://") {
							// We drop between two files
							if (fpath.ends_with("/")) {
								fpath = fpath.substr(0, fpath.length() - 1);
							}
							target = fpath.get_base_dir();
							return;
						}
					}
				}
			}
		}
	}

	return;
}

void FileSystemDock::_file_and_folders_fill_popup(PopupMenu *p_popup, Vector<String> p_paths) {
	// Add options for files and folders
	ERR_FAIL_COND(p_paths.empty())

	Vector<String> filenames;
	Vector<String> foldernames;

	Vector<String> favorites = EditorSettings::get_singleton()->get_favorite_dirs();

	bool all_files = true;
	bool all_files_scenes = true;
	bool all_folders = true;
	bool all_favorites = true;
	bool all_not_favorites = true;
	for (int i = 0; i < p_paths.size(); i++) {
		String fpath = p_paths[i];
		if (fpath.ends_with("/")) {
			foldernames.push_back(fpath);
			all_files = false;
		} else {
			filenames.push_back(fpath);
			all_folders = false;
			all_files_scenes &= (EditorFileSystem::get_singleton()->get_file_type(fpath) == "PackedScene");
		}

		// Check if in favorites
		bool found = false;
		for (int j = 0; j < favorites.size(); j++) {
			if (favorites[j] == fpath) {
				found = true;
				break;
			}
		}
		if (found) {
			all_not_favorites = false;
		} else {
			all_favorites = false;
		}
	}

	if (all_files) {

		if (all_files_scenes && filenames.size() >= 1) {
			p_popup->add_item(TTR("Open Scene(s)"), FILE_OPEN);
			p_popup->add_item(TTR("Instance"), FILE_INSTANCE);
			p_popup->add_separator();
		}

		if (!all_files_scenes && filenames.size() == 1) {
			p_popup->add_item(TTR("Open"), FILE_OPEN);
			p_popup->add_separator();
		}
	}

	if (p_paths.size() >= 1) {
		if (!all_favorites) {
			p_popup->add_item(TTR("Add to favorites"), FILE_ADD_FAVORITE);
		}
		if (!all_not_favorites) {
			p_popup->add_item(TTR("Remove from favorites"), FILE_REMOVE_FAVORITE);
		}
		p_popup->add_separator();
	}

	if (all_files) {
		if (filenames.size() == 1) {
			p_popup->add_item(TTR("Edit Dependencies..."), FILE_DEPENDENCIES);
			p_popup->add_item(TTR("View Owners..."), FILE_OWNERS);
			p_popup->add_separator();
		}

	} else if (all_folders && foldernames.size() > 0) {
		p_popup->add_item(TTR("Open"), FILE_OPEN);
		p_popup->add_separator();
	}

	if (p_paths.size() == 1) {
		p_popup->add_item(TTR("Copy Path"), FILE_COPY_PATH);
		p_popup->add_item(TTR("Rename..."), FILE_RENAME);
		p_popup->add_item(TTR("Duplicate..."), FILE_DUPLICATE);
	}

	p_popup->add_item(TTR("Move To..."), FILE_MOVE);
	p_popup->add_item(TTR("Delete"), FILE_REMOVE);

	if (p_paths.size() == 1) {
		p_popup->add_separator();
		p_popup->add_item(TTR("New Folder..."), FILE_NEW_FOLDER);
		p_popup->add_item(TTR("New Script..."), FILE_NEW_SCRIPT);
		p_popup->add_item(TTR("New Resource..."), FILE_NEW_RESOURCE);

		String fpath = p_paths[0];
		String item_text = fpath.ends_with("/") ? TTR("Open In File Manager") : TTR("Show In File Manager");
		p_popup->add_item(item_text, FILE_SHOW_IN_EXPLORER);
	}
}

void FileSystemDock::_tree_rmb_select(const Vector2 &p_pos) {
	// Right click is pressed in the tree
	Vector<String> paths = _tree_get_selected();

	if (paths.size() == 1) {
		if (paths[0].ends_with("/")) {
			tree_popup->add_item(TTR("Expand all"), FOLDER_EXPAND_ALL);
			tree_popup->add_item(TTR("Collapse all"), FOLDER_COLLAPSE_ALL);
			tree_popup->add_separator();
		}
	}

	// Popup
	if (!paths.empty()) {
		tree_popup->clear();
		tree_popup->set_size(Size2(1, 1));
		_file_and_folders_fill_popup(tree_popup, paths);
		tree_popup->set_position(tree->get_global_position() + p_pos);
		tree_popup->popup();
	}
}

void FileSystemDock::_file_list_rmb_select(int p_item, const Vector2 &p_pos) {
	// Right click is pressed in the file list
	Vector<String> paths;
	for (int i = 0; i < files->get_item_count(); i++) {
		if (!files->is_selected(i))
			continue;
		if (files->get_item_text(p_item) == "..") {
			files->unselect(i);
			continue;
		}
		paths.push_back(files->get_item_metadata(i));
	}

	// Popup
	if (!paths.empty()) {
		file_list_popup->clear();
		file_list_popup->set_size(Size2(1, 1));
		_file_and_folders_fill_popup(file_list_popup, paths);
		file_list_popup->set_position(files->get_global_position() + p_pos);
		file_list_popup->popup();
	}
}

void FileSystemDock::_file_list_rmb_pressed(const Vector2 &p_pos) {
	// Right click on empty space for file list
	file_list_popup->clear();
	file_list_popup->set_size(Size2(1, 1));

	file_list_popup->add_item(TTR("New Folder..."), FILE_NEW_FOLDER);
	file_list_popup->add_item(TTR("New Script..."), FILE_NEW_SCRIPT);
	file_list_popup->add_item(TTR("New Resource..."), FILE_NEW_RESOURCE);
	file_list_popup->add_item(TTR("Show In File Manager"), FILE_SHOW_IN_EXPLORER);
	file_list_popup->set_position(files->get_global_position() + p_pos);
	file_list_popup->popup();
}

void FileSystemDock::select_file(const String &p_file) {

	navigate_to_path(p_file);
}

void FileSystemDock::_file_multi_selected(int p_index, bool p_selected) {

	// Set the path to the current focussed item
	int current = files->get_current();
	if (current == p_index) {
		String fpath = files->get_item_metadata(current);
		if (!fpath.ends_with("/")) {
			path = fpath;
			if (display_mode == DISPLAY_MODE_SPLIT) {
				_update_tree(_compute_uncollapsed_paths());
			}
		}
	}

	// Update the import dock
	import_dock_needs_update = true;
	call_deferred("_update_import_dock");
}

void FileSystemDock::_tree_gui_input(Ref<InputEvent> p_event) {

	if (get_viewport()->get_modal_stack_top())
		return; //ignore because of modal window

	Ref<InputEventKey> key = p_event;
	if (key.is_valid() && key->is_pressed() && !key->is_echo()) {
		if (ED_IS_SHORTCUT("filesystem_dock/duplicate", p_event)) {
			_tree_rmb_option(FILE_DUPLICATE);
		} else if (ED_IS_SHORTCUT("filesystem_dock/copy_path", p_event)) {
			_tree_rmb_option(FILE_COPY_PATH);
		} else if (ED_IS_SHORTCUT("filesystem_dock/delete", p_event)) {
			_tree_rmb_option(FILE_REMOVE);
		} else if (ED_IS_SHORTCUT("filesystem_dock/rename", p_event)) {
			_tree_rmb_option(FILE_RENAME);
		}
	}
}

void FileSystemDock::_file_list_gui_input(Ref<InputEvent> p_event) {

	if (get_viewport()->get_modal_stack_top())
		return; //ignore because of modal window

	Ref<InputEventKey> key = p_event;
	if (key.is_valid() && key->is_pressed() && !key->is_echo()) {
		if (ED_IS_SHORTCUT("filesystem_dock/duplicate", p_event)) {
			_file_list_rmb_option(FILE_DUPLICATE);
		} else if (ED_IS_SHORTCUT("filesystem_dock/copy_path", p_event)) {
			_file_list_rmb_option(FILE_COPY_PATH);
		} else if (ED_IS_SHORTCUT("filesystem_dock/delete", p_event)) {
			_file_list_rmb_option(FILE_REMOVE);
		} else if (ED_IS_SHORTCUT("filesystem_dock/rename", p_event)) {
			_file_list_rmb_option(FILE_RENAME);
		}
	}
}

void FileSystemDock::_update_import_dock() {

	if (!import_dock_needs_update)
		return;

	//check import
	Vector<String> imports;
	String import_type;

	for (int i = 0; i < files->get_item_count(); i++) {
		if (!files->is_selected(i))
			continue;

		String fpath = files->get_item_metadata(i);
		if (!FileAccess::exists(fpath + ".import")) {
			imports.clear();
			break;
		}
		Ref<ConfigFile> cf;
		cf.instance();
		Error err = cf->load(fpath + ".import");
		if (err != OK) {
			imports.clear();
			break;
		}

		String type = cf->get_value("remap", "type");
		if (import_type == "") {
			import_type = type;
		} else if (import_type != type) {
			//all should be the same type
			imports.clear();
			break;
		}
		imports.push_back(fpath);
	}

	if (imports.size() == 0) {
		EditorNode::get_singleton()->get_import_dock()->clear();
	} else if (imports.size() == 1) {
		EditorNode::get_singleton()->get_import_dock()->set_edit_path(imports[0]);
	} else {
		EditorNode::get_singleton()->get_import_dock()->set_edit_multiple_paths(imports);
	}

	import_dock_needs_update = false;
}

void FileSystemDock::_bind_methods() {

	ClassDB::bind_method(D_METHOD("_file_list_gui_input"), &FileSystemDock::_file_list_gui_input);
	ClassDB::bind_method(D_METHOD("_tree_gui_input"), &FileSystemDock::_tree_gui_input);

	ClassDB::bind_method(D_METHOD("_update_tree"), &FileSystemDock::_update_tree);
	ClassDB::bind_method(D_METHOD("_rescan"), &FileSystemDock::_rescan);

	ClassDB::bind_method(D_METHOD("_toggle_split_mode"), &FileSystemDock::_toggle_split_mode);

	ClassDB::bind_method(D_METHOD("_tree_rmb_option", "option"), &FileSystemDock::_tree_rmb_option);
	ClassDB::bind_method(D_METHOD("_file_list_rmb_option", "option"), &FileSystemDock::_file_list_rmb_option);

	ClassDB::bind_method(D_METHOD("_tree_rmb_select"), &FileSystemDock::_tree_rmb_select);
	ClassDB::bind_method(D_METHOD("_file_list_rmb_select"), &FileSystemDock::_file_list_rmb_select);
	ClassDB::bind_method(D_METHOD("_file_list_rmb_pressed"), &FileSystemDock::_file_list_rmb_pressed);

	ClassDB::bind_method(D_METHOD("_file_list_thumbnail_done"), &FileSystemDock::_file_list_thumbnail_done);
	ClassDB::bind_method(D_METHOD("_tree_thumbnail_done"), &FileSystemDock::_tree_thumbnail_done);
	ClassDB::bind_method(D_METHOD("_file_list_activate_file"), &FileSystemDock::_file_list_activate_file);
	ClassDB::bind_method(D_METHOD("_tree_activate_file"), &FileSystemDock::_tree_activate_file);
	ClassDB::bind_method(D_METHOD("_select_file"), &FileSystemDock::_select_file);
	ClassDB::bind_method(D_METHOD("_go_to_tree"), &FileSystemDock::_go_to_tree);
	ClassDB::bind_method(D_METHOD("navigate_to_path"), &FileSystemDock::navigate_to_path);
	ClassDB::bind_method(D_METHOD("_change_file_display"), &FileSystemDock::_change_file_display);
	ClassDB::bind_method(D_METHOD("_fw_history"), &FileSystemDock::_fw_history);
	ClassDB::bind_method(D_METHOD("_bw_history"), &FileSystemDock::_bw_history);
	ClassDB::bind_method(D_METHOD("_fs_changed"), &FileSystemDock::_fs_changed);
	ClassDB::bind_method(D_METHOD("_tree_multi_selected"), &FileSystemDock::_tree_multi_selected);
	ClassDB::bind_method(D_METHOD("_make_dir_confirm"), &FileSystemDock::_make_dir_confirm);
	ClassDB::bind_method(D_METHOD("_resource_created"), &FileSystemDock::_resource_created);
	ClassDB::bind_method(D_METHOD("_move_operation_confirm", "to_path", "overwrite"), &FileSystemDock::_move_operation_confirm, DEFVAL(false));
	ClassDB::bind_method(D_METHOD("_move_with_overwrite"), &FileSystemDock::_move_with_overwrite);
	ClassDB::bind_method(D_METHOD("_rename_operation_confirm"), &FileSystemDock::_rename_operation_confirm);
	ClassDB::bind_method(D_METHOD("_duplicate_operation_confirm"), &FileSystemDock::_duplicate_operation_confirm);

	ClassDB::bind_method(D_METHOD("_search_changed"), &FileSystemDock::_search_changed);

	ClassDB::bind_method(D_METHOD("get_drag_data_fw"), &FileSystemDock::get_drag_data_fw);
	ClassDB::bind_method(D_METHOD("can_drop_data_fw"), &FileSystemDock::can_drop_data_fw);
	ClassDB::bind_method(D_METHOD("drop_data_fw"), &FileSystemDock::drop_data_fw);

	ClassDB::bind_method(D_METHOD("_preview_invalidated"), &FileSystemDock::_preview_invalidated);
	ClassDB::bind_method(D_METHOD("_file_multi_selected"), &FileSystemDock::_file_multi_selected);
	ClassDB::bind_method(D_METHOD("_update_import_dock"), &FileSystemDock::_update_import_dock);

	ADD_SIGNAL(MethodInfo("instance", PropertyInfo(Variant::POOL_STRING_ARRAY, "files")));
	ADD_SIGNAL(MethodInfo("open"));
}

FileSystemDock::FileSystemDock(EditorNode *p_editor) {

	set_name("FileSystem");
	editor = p_editor;
	path = "res://";

	ED_SHORTCUT("filesystem_dock/copy_path", TTR("Copy Path"), KEY_MASK_CMD | KEY_C);
	ED_SHORTCUT("filesystem_dock/duplicate", TTR("Duplicate..."), KEY_MASK_CMD | KEY_D);
	ED_SHORTCUT("filesystem_dock/delete", TTR("Delete"), KEY_DELETE);
	ED_SHORTCUT("filesystem_dock/rename", TTR("Rename"));

	VBoxContainer *top_vbc = memnew(VBoxContainer);
	add_child(top_vbc);

	HBoxContainer *toolbar_hbc = memnew(HBoxContainer);
	toolbar_hbc->add_constant_override("separation", 0);
	top_vbc->add_child(toolbar_hbc);

	button_hist_prev = memnew(ToolButton);
	button_hist_prev->set_disabled(true);
	button_hist_prev->set_focus_mode(FOCUS_NONE);
	button_hist_prev->set_tooltip(TTR("Previous Directory"));
	toolbar_hbc->add_child(button_hist_prev);

	button_hist_next = memnew(ToolButton);
	button_hist_next->set_disabled(true);
	button_hist_next->set_focus_mode(FOCUS_NONE);
	button_hist_next->set_tooltip(TTR("Next Directory"));
	toolbar_hbc->add_child(button_hist_next);

	current_path = memnew(LineEdit);
	current_path->set_h_size_flags(SIZE_EXPAND_FILL);
	toolbar_hbc->add_child(current_path);

	button_reload = memnew(Button);
	button_reload->set_flat(true);
	button_reload->connect("pressed", this, "_rescan");
	button_reload->set_focus_mode(FOCUS_NONE);
	button_reload->set_tooltip(TTR("Re-Scan Filesystem"));
	button_reload->hide();
	toolbar_hbc->add_child(button_reload);

	button_toggle_display_mode = memnew(Button);
	button_toggle_display_mode->set_flat(true);
	button_toggle_display_mode->set_toggle_mode(true);
	button_toggle_display_mode->connect("toggled", this, "_toggle_split_mode");
	button_toggle_display_mode->set_focus_mode(FOCUS_NONE);
	button_toggle_display_mode->set_tooltip(TTR("Toggle split mode"));
	toolbar_hbc->add_child(button_toggle_display_mode);

	HBoxContainer *toolbar2_hbc = memnew(HBoxContainer);
	toolbar2_hbc->add_constant_override("separation", 0);
	top_vbc->add_child(toolbar2_hbc);

	tree_search_box = memnew(LineEdit);
	tree_search_box->set_h_size_flags(SIZE_EXPAND_FILL);
	tree_search_box->set_placeholder(TTR("Search files"));
	tree_search_box->connect("text_changed", this, "_search_changed", varray(tree_search_box));
	toolbar2_hbc->add_child(tree_search_box);

	//toolbar_hbc->add_spacer();

	//Control *spacer = memnew( Control);

	/*
	button_open = memnew( Button );
	button_open->set_flat(true);
	button_open->connect("pressed",this,"_go_to_file_list");
	toolbar_hbc->add_child(button_open);
	button_open->hide();
	button_open->set_focus_mode(FOCUS_NONE);
	button_open->set_tooltip("Open the selected file.\nOpen as scene if a scene, or as resource otherwise.");


	button_instance = memnew( Button );
	button_instance->set_flat(true);
	button_instance->connect("pressed",this,"_instance_pressed");
	toolbar_hbc->add_child(button_instance);
	button_instance->hide();
	button_instance->set_focus_mode(FOCUS_NONE);
	button_instance->set_tooltip(TTR("Instance the selected scene(s) as child of the selected node."));

*/
	file_list_popup = memnew(PopupMenu);
	file_list_popup->set_hide_on_window_lose_focus(true);
	add_child(file_list_popup);

	tree_popup = memnew(PopupMenu);
	tree_popup->set_hide_on_window_lose_focus(true);
	add_child(tree_popup);

	split_box = memnew(VSplitContainer);
	split_box->set_v_size_flags(SIZE_EXPAND_FILL);
	add_child(split_box);

	tree = memnew(Tree);

	tree->set_hide_root(true);
	tree->set_drag_forwarding(this);
	tree->set_allow_rmb_select(true);
	tree->set_select_mode(Tree::SELECT_MULTI);
	tree->set_custom_minimum_size(Size2(0, 200 * EDSCALE));
	split_box->add_child(tree);

	tree->connect("item_edited", this, "_favorite_toggled");
	tree->connect("item_activated", this, "_tree_activate_file");
	tree->connect("multi_selected", this, "_tree_multi_selected");
	tree->connect("item_rmb_selected", this, "_tree_rmb_select");
	tree->connect("gui_input", this, "_tree_gui_input");

	file_list_vb = memnew(VBoxContainer);
	file_list_vb->set_v_size_flags(SIZE_EXPAND_FILL);
	split_box->add_child(file_list_vb);

	path_hb = memnew(HBoxContainer);
	file_list_vb->add_child(path_hb);

	button_tree = memnew(ToolButton);
	button_tree->set_tooltip(TTR("Enter tree-view."));
	button_tree->hide();
	path_hb->add_child(button_tree);

	file_list_search_box = memnew(LineEdit);
	file_list_search_box->set_h_size_flags(SIZE_EXPAND_FILL);
	file_list_search_box->set_placeholder(TTR("Search files"));
	file_list_search_box->connect("text_changed", this, "_search_changed", varray(file_list_search_box));
	path_hb->add_child(file_list_search_box);

	button_file_list_display_mode = memnew(ToolButton);
	button_file_list_display_mode->set_toggle_mode(true);
	path_hb->add_child(button_file_list_display_mode);

	files = memnew(ItemList);
	files->set_v_size_flags(SIZE_EXPAND_FILL);
	files->set_select_mode(ItemList::SELECT_MULTI);
	files->set_drag_forwarding(this);
	files->connect("item_rmb_selected", this, "_file_list_rmb_select");
	files->connect("gui_input", this, "_file_list_gui_input");
	files->connect("multi_selected", this, "_file_multi_selected");
	files->connect("rmb_clicked", this, "_file_list_rmb_pressed");
	files->set_allow_rmb_select(true);
	file_list_vb->add_child(files);

	scanning_vb = memnew(VBoxContainer);
	scanning_vb->hide();
	add_child(scanning_vb);

	Label *slabel = memnew(Label);
	slabel->set_text(TTR("Scanning Files,\nPlease Wait..."));
	slabel->set_align(Label::ALIGN_CENTER);
	scanning_vb->add_child(slabel);

	scanning_progress = memnew(ProgressBar);
	scanning_vb->add_child(scanning_progress);

	deps_editor = memnew(DependencyEditor);
	add_child(deps_editor);

	owners_editor = memnew(DependencyEditorOwners(editor));
	add_child(owners_editor);

	remove_dialog = memnew(DependencyRemoveDialog);
	add_child(remove_dialog);

	move_dialog = memnew(EditorDirDialog);
	move_dialog->get_ok()->set_text(TTR("Move"));
	add_child(move_dialog);
	move_dialog->connect("dir_selected", this, "_move_operation_confirm");

	rename_dialog = memnew(ConfirmationDialog);
	VBoxContainer *rename_dialog_vb = memnew(VBoxContainer);
	rename_dialog->add_child(rename_dialog_vb);

	rename_dialog_text = memnew(LineEdit);
	rename_dialog_vb->add_margin_child(TTR("Name:"), rename_dialog_text);
	rename_dialog->get_ok()->set_text(TTR("Rename"));
	add_child(rename_dialog);
	rename_dialog->register_text_enter(rename_dialog_text);
	rename_dialog->connect("confirmed", this, "_rename_operation_confirm");

	overwrite_dialog = memnew(ConfirmationDialog);
	overwrite_dialog->set_text(TTR("There is already file or folder with the same name in this location."));
	overwrite_dialog->get_ok()->set_text(TTR("Overwrite"));
	add_child(overwrite_dialog);
	overwrite_dialog->connect("confirmed", this, "_move_with_overwrite");

	duplicate_dialog = memnew(ConfirmationDialog);
	VBoxContainer *duplicate_dialog_vb = memnew(VBoxContainer);
	duplicate_dialog->add_child(duplicate_dialog_vb);

	duplicate_dialog_text = memnew(LineEdit);
	duplicate_dialog_vb->add_margin_child(TTR("Name:"), duplicate_dialog_text);
	duplicate_dialog->get_ok()->set_text(TTR("Duplicate"));
	add_child(duplicate_dialog);
	duplicate_dialog->register_text_enter(duplicate_dialog_text);
	duplicate_dialog->connect("confirmed", this, "_duplicate_operation_confirm");

	make_dir_dialog = memnew(ConfirmationDialog);
	make_dir_dialog->set_title(TTR("Create Folder"));
	VBoxContainer *make_folder_dialog_vb = memnew(VBoxContainer);
	make_dir_dialog->add_child(make_folder_dialog_vb);

	make_dir_dialog_text = memnew(LineEdit);
	make_folder_dialog_vb->add_margin_child(TTR("Name:"), make_dir_dialog_text);
	add_child(make_dir_dialog);
	make_dir_dialog->register_text_enter(make_dir_dialog_text);
	make_dir_dialog->connect("confirmed", this, "_make_dir_confirm");

	make_script_dialog_text = memnew(ScriptCreateDialog);
	make_script_dialog_text->set_title(TTR("Create Script"));
	add_child(make_script_dialog_text);

	new_resource_dialog = memnew(CreateDialog);
	add_child(new_resource_dialog);
	new_resource_dialog->set_base_type("Resource");
	new_resource_dialog->connect("create", this, "_resource_created");

	searched_string = String();
	uncollapsed_paths_before_search = Vector<String>();

	updating_tree = false;
	tree_update_id = 0;
	initialized = false;
	import_dock_needs_update = false;

	history_pos = 0;
	history_max_size = 20;
	history.push_back("res://");

	display_mode = DISPLAY_MODE_SPLIT;
	file_list_display_mode = FILE_LIST_DISPLAY_THUMBNAILS;

	display_mode_setting = DISPLAY_MODE_SETTING_TREE_ONLY;
	old_display_mode_setting = DISPLAY_MODE_SETTING_TREE_ONLY;
	always_show_folders = false;
}

FileSystemDock::~FileSystemDock() {
}
