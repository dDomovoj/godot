/*************************************************************************/
/*  input_map.cpp                                                        */
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
#include "input_map.h"
#include "globals.h"

InputMap *InputMap::singleton=NULL;

void InputMap::_bind_methods() {

	ObjectTypeDB::bind_method(_MD("has_action","action"),&InputMap::has_action);
	ObjectTypeDB::bind_method(_MD("get_action_id","action"),&InputMap::get_action_id);
	ObjectTypeDB::bind_method(_MD("get_action_from_id","id"),&InputMap::get_action_from_id);
	ObjectTypeDB::bind_method(_MD("add_action","action"),&InputMap::add_action);
	ObjectTypeDB::bind_method(_MD("erase_action","action"),&InputMap::erase_action);
	
	ObjectTypeDB::bind_method(_MD("action_add_event","action","event"),&InputMap::action_add_event);
	ObjectTypeDB::bind_method(_MD("action_has_event","action","event"),&InputMap::action_has_event);
	ObjectTypeDB::bind_method(_MD("action_erase_event","action","event"),&InputMap::action_erase_event);
	ObjectTypeDB::bind_method(_MD("get_action_list","action"),&InputMap::_get_action_list);
	ObjectTypeDB::bind_method(_MD("event_is_action","event","action"),&InputMap::event_is_action);
	ObjectTypeDB::bind_method(_MD("load_from_globals"),&InputMap::load_from_globals);

}


void InputMap::add_action(const StringName& p_action) {

	ERR_FAIL_COND( input_map.has(p_action) );
	input_map[p_action]=Action();
	static int last_id=1;
	input_map[p_action].id=last_id;
	input_id_map[last_id]=p_action;
	last_id++;


}

void InputMap::erase_action(const StringName& p_action) {

	ERR_FAIL_COND( !input_map.has(p_action) );
	input_id_map.erase(input_map[p_action].id);
	input_map.erase(p_action);

}

StringName InputMap::get_action_from_id(int p_id) const {

	ERR_FAIL_COND_V(!input_id_map.has(p_id),StringName());
	return input_id_map[p_id];
}

List<InputEvent>::Element *InputMap::_find_event(List<InputEvent> &p_list,const InputEvent& p_event) const {

	for (List<InputEvent>::Element *E=p_list.front();E;E=E->next()) {

		const InputEvent& e=E->get();
		if(e.type!=p_event.type)
			continue;
		if (e.type!=InputEvent::KEY && e.device!=p_event.device)
			continue;

		bool same=false;

		switch(p_event.type) {

			case InputEvent::KEY: {

				same=(e.key.scancode==p_event.key.scancode && e.key.mod == p_event.key.mod);

			} break;
			case InputEvent::JOYSTICK_BUTTON: {

				same=(e.joy_button.button_index==p_event.joy_button.button_index);

			} break;
			case InputEvent::MOUSE_BUTTON: {

				same=(e.mouse_button.button_index==p_event.mouse_button.button_index);

			} break;
			case InputEvent::JOYSTICK_MOTION: {

				same=(e.joy_motion.axis==p_event.joy_motion.axis);

			} break;
		}

		if (same)
			return E;
	}


	return NULL;
}


bool InputMap::has_action(const StringName& p_action) const {

	return input_map.has(p_action);
}

void InputMap::action_add_event(const StringName& p_action,const InputEvent& p_event) {

	ERR_FAIL_COND(p_event.type==InputEvent::ACTION);
	ERR_FAIL_COND( !input_map.has(p_action) );
	if (_find_event(input_map[p_action].inputs,p_event))
		return; //already gots

	input_map[p_action].inputs.push_back(p_event);

}


int InputMap::get_action_id(const StringName& p_action) const {

	ERR_FAIL_COND_V(!input_map.has(p_action),-1);
	return input_map[p_action].id;
}

bool InputMap::action_has_event(const StringName& p_action,const InputEvent& p_event) {

	ERR_FAIL_COND_V( !input_map.has(p_action), false );
	return (_find_event(input_map[p_action].inputs,p_event)!=NULL);

}

void InputMap::action_erase_event(const StringName& p_action,const InputEvent& p_event) {

	ERR_FAIL_COND( !input_map.has(p_action) );

	List<InputEvent>::Element *E=_find_event(input_map[p_action].inputs,p_event);
	if (E)
		input_map[p_action].inputs.erase(E);
}


Array InputMap::_get_action_list(const StringName& p_action) {

	Array ret;
	const List<InputEvent> *al = get_action_list(p_action);
	if (al) {
		for(const List<InputEvent>::Element *E=al->front();E;E=E->next()) {

			ret.push_back(E->get());;
		}
	}

	return ret;

}

const List<InputEvent> *InputMap::get_action_list(const StringName& p_action) {

	const Map<StringName, Action>::Element *E=input_map.find(p_action);
	if (!E)
		return NULL;

	return &E->get().inputs;
}

bool InputMap::event_is_action(const InputEvent& p_event, const StringName& p_action) const {


	Map<StringName,Action >::Element *E=input_map.find(p_action);
	if(!E) {
		ERR_EXPLAIN("Request for nonexistent InputMap action: "+String(p_action));
		ERR_FAIL_COND_V(!E,false);
	}

	if (p_event.type==InputEvent::ACTION) {

		return p_event.action.action==E->get().id;
	}

	return _find_event(E->get().inputs,p_event)!=NULL;
}

bool InputMap::event_is_joy_motion_action_pressed(const InputEvent& p_event) const {

	ERR_FAIL_COND_V(p_event.type!=InputEvent::JOYSTICK_MOTION,false);
	bool pressed=false;

	//this could be optimized by having a separate list of joymotions?

	for (Map<StringName, Action>::Element *A=input_map.front();A;A=A->next()) {

		for (List<InputEvent>::Element *E=A->get().inputs.front();E;E=E->next()) {

			const InputEvent& e=E->get();
			if(e.type!=p_event.type)
				continue;
			if (e.type!=InputEvent::KEY && e.device!=p_event.device)
				continue;

			switch(p_event.type) {

				case InputEvent::KEY: {

					if (e.key.scancode==p_event.key.scancode && e.key.mod == p_event.key.mod)
						return e.key.pressed;

				} break;
				case InputEvent::JOYSTICK_BUTTON: {

					if (e.joy_button.button_index==p_event.joy_button.button_index) {
						return e.joy_button.pressed;
					}

				} break;
				case InputEvent::MOUSE_BUTTON: {

					if (e.mouse_button.button_index==p_event.mouse_button.button_index) {
						return e.mouse_button.pressed;
					}

				} break;
				case InputEvent::JOYSTICK_MOTION: {

					if (e.joy_motion.axis==p_event.joy_motion.axis) {
						if (
								(e.joy_motion.axis_value * p_event.joy_motion.axis_value >0) && //same axis
								ABS(e.joy_motion.axis_value)>0.5 && ABS(p_event.joy_motion.axis_value)>0.5 )
							pressed=true;
					}

				} break;
			}

		}
	}

	return pressed;

}

void InputMap::load_from_globals() {

	input_map.clear();;

	List<PropertyInfo> pinfo;
	Globals::get_singleton()->get_property_list(&pinfo);

	for(List<PropertyInfo>::Element *E=pinfo.front();E;E=E->next()) {
		const PropertyInfo &pi=E->get();

		if (!pi.name.begins_with("input/"))
			continue;

		String name = pi.name.substr(pi.name.find("/")+1,pi.name.length());

		add_action(name);

		Array va = Globals::get_singleton()->get(pi.name);;

		for(int i=0;i<va.size();i++) {

			InputEvent ie=va[i];
			if (ie.type==InputEvent::NONE)
				continue;
			action_add_event(name,ie);

		}

	}

}

InputMap::InputMap() {

	ERR_FAIL_COND(singleton);
	singleton=this;
}
