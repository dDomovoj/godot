#include "visual_script_yield_nodes.h"
#include "scene/main/scene_main_loop.h"
#include "os/os.h"
#include "scene/main/node.h"
#include "visual_script_nodes.h"

//////////////////////////////////////////
////////////////YIELD///////////
//////////////////////////////////////////

int VisualScriptYield::get_output_sequence_port_count() const {

	return 1;
}

bool VisualScriptYield::has_input_sequence_port() const{

	return true;
}

int VisualScriptYield::get_input_value_port_count() const{

	return 0;
}
int VisualScriptYield::get_output_value_port_count() const{

	return 0;
}

String VisualScriptYield::get_output_sequence_port_text(int p_port) const {

	return String();
}

PropertyInfo VisualScriptYield::get_input_value_port_info(int p_idx) const{

	return PropertyInfo();
}

PropertyInfo VisualScriptYield::get_output_value_port_info(int p_idx) const{

	return PropertyInfo();
}


String VisualScriptYield::get_caption() const {

	return "Wait";
}

String VisualScriptYield::get_text() const {

	switch (yield_mode) {
		case YIELD_FRAME: return "Next Frame"; break;
		case YIELD_FIXED_FRAME:  return "Next Fixed Frame"; break;
		case YIELD_WAIT:  return rtos(wait_time)+" sec(s)"; break;
	}

	return String();
}


class VisualScriptNodeInstanceYield : public VisualScriptNodeInstance {
public:

	VisualScriptYield::YieldMode mode;
	float wait_time;

	virtual int get_working_memory_size() const { return 1; } //yield needs at least 1
	//virtual bool is_output_port_unsequenced(int p_idx) const { return false; }
	//virtual bool get_output_port_unsequenced(int p_idx,Variant* r_value,Variant* p_working_mem,String &r_error) const { return false; }

	virtual int step(const Variant** p_inputs,Variant** p_outputs,StartMode p_start_mode,Variant* p_working_mem,Variant::CallError& r_error,String& r_error_str) {

		if (p_start_mode==START_MODE_RESUME_YIELD) {
			return 0; //resuming yield
		} else {
			//yield


			SceneTree *tree = OS::get_singleton()->get_main_loop()->cast_to<SceneTree>();
			if (!tree) {
				r_error_str="Main Loop is not SceneTree";
				r_error.error=Variant::CallError::CALL_ERROR_INVALID_METHOD;
				return 0;
			}

			Ref<VisualScriptFunctionState> state;
			state.instance();

			switch(mode) {

				case VisualScriptYield::YIELD_FRAME: state->connect_to_signal(tree,"idle_frame",Array()); break;
				case VisualScriptYield::YIELD_FIXED_FRAME:  state->connect_to_signal(tree,"fixed_frame",Array()); break;
				case VisualScriptYield::YIELD_WAIT:  state->connect_to_signal(tree->create_timer(wait_time).ptr(),"timeout",Array()); break;

			}

			*p_working_mem=state;

			return STEP_YIELD_BIT;
		}
	}

};

VisualScriptNodeInstance* VisualScriptYield::instance(VisualScriptInstance* p_instance) {

	VisualScriptNodeInstanceYield * instance = memnew(VisualScriptNodeInstanceYield );
	//instance->instance=p_instance;
	instance->mode=yield_mode;
	instance->wait_time=wait_time;
	return instance;
}

void VisualScriptYield::set_yield_mode(YieldMode p_mode) {

	if (yield_mode==p_mode)
		return;
	yield_mode=p_mode;
	ports_changed_notify();
	_change_notify();
}

VisualScriptYield::YieldMode VisualScriptYield::get_yield_mode(){

	return yield_mode;
}

void VisualScriptYield::set_wait_time(float p_time) {

	if (wait_time==p_time)
		return;
	wait_time=p_time;
	ports_changed_notify();

}

float VisualScriptYield::get_wait_time(){

	return wait_time;
}


void VisualScriptYield::_validate_property(PropertyInfo& property) const {


	if (property.name=="wait_time") {
		if (yield_mode!=YIELD_WAIT) {
			property.usage=0;
		}
	}
}

void VisualScriptYield::_bind_methods() {

	ObjectTypeDB::bind_method(_MD("set_yield_mode","mode"),&VisualScriptYield::set_yield_mode);
	ObjectTypeDB::bind_method(_MD("get_yield_mode"),&VisualScriptYield::get_yield_mode);

	ObjectTypeDB::bind_method(_MD("set_wait_time","sec"),&VisualScriptYield::set_wait_time);
	ObjectTypeDB::bind_method(_MD("get_wait_time"),&VisualScriptYield::get_wait_time);

	ADD_PROPERTY(PropertyInfo(Variant::INT,"mode",PROPERTY_HINT_ENUM,"Frame,FixedFrame,Time",PROPERTY_USAGE_NOEDITOR),_SCS("set_yield_mode"),_SCS("get_yield_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::REAL,"wait_time"),_SCS("set_wait_time"),_SCS("get_wait_time"));


	BIND_CONSTANT( YIELD_FRAME );
	BIND_CONSTANT( YIELD_FIXED_FRAME );
	BIND_CONSTANT( YIELD_WAIT );

}

VisualScriptYield::VisualScriptYield() {

	yield_mode=YIELD_FRAME;
	wait_time=1;

}


template<VisualScriptYield::YieldMode MODE>
static Ref<VisualScriptNode> create_yield_node(const String& p_name) {

	Ref<VisualScriptYield> node;
	node.instance();
	node->set_yield_mode(MODE);
	return node;
}

///////////////////////////////////////////////////
////////////////YIELD SIGNAL//////////////////////
//////////////////////////////////////////////////

int VisualScriptYieldSignal::get_output_sequence_port_count() const {

	return 1;
}

bool VisualScriptYieldSignal::has_input_sequence_port() const{

	return true;
}
#ifdef TOOLS_ENABLED

static Node* _find_script_node(Node* p_edited_scene,Node* p_current_node,const Ref<Script> &script) {

	if (p_edited_scene!=p_current_node && p_current_node->get_owner()!=p_edited_scene)
		return NULL;

	Ref<Script> scr = p_current_node->get_script();

	if (scr.is_valid() && scr==script)
		return p_current_node;

	for(int i=0;i<p_current_node->get_child_count();i++) {
		Node *n = _find_script_node(p_edited_scene,p_current_node->get_child(i),script);
		if (n)
			return n;
	}

	return NULL;
}

#endif
Node *VisualScriptYieldSignal::_get_base_node() const {

#ifdef TOOLS_ENABLED
	Ref<Script> script = get_visual_script();
	if (!script.is_valid())
		return NULL;

	MainLoop * main_loop = OS::get_singleton()->get_main_loop();
	if (!main_loop)
		return NULL;

	SceneTree *scene_tree = main_loop->cast_to<SceneTree>();

	if (!scene_tree)
		return NULL;

	Node *edited_scene = scene_tree->get_edited_scene_root();

	if (!edited_scene)
		return NULL;

	Node* script_node = _find_script_node(edited_scene,edited_scene,script);

	if (!script_node)
		return NULL;

	if (!script_node->has_node(base_path))
		return NULL;

	Node *path_to = script_node->get_node(base_path);

	return path_to;
#else

	return NULL;
#endif
}

StringName VisualScriptYieldSignal::_get_base_type() const {

	if (call_mode==CALL_MODE_SELF && get_visual_script().is_valid())
		return get_visual_script()->get_instance_base_type();
	else if (call_mode==CALL_MODE_NODE_PATH && get_visual_script().is_valid()) {
		Node *path = _get_base_node();
		if (path)
			return path->get_type();

	}

	return base_type;
}

int VisualScriptYieldSignal::get_input_value_port_count() const{

	if (call_mode==CALL_MODE_INSTANCE)
		return 1;
	else
		return 0;

}
int VisualScriptYieldSignal::get_output_value_port_count() const{


	MethodInfo sr;

	if (!ObjectTypeDB::get_signal(_get_base_type(),signal,&sr))
		return 0;

	return sr.arguments.size();

}

String VisualScriptYieldSignal::get_output_sequence_port_text(int p_port) const {

	return String();
}

PropertyInfo VisualScriptYieldSignal::get_input_value_port_info(int p_idx) const{

	if (call_mode==CALL_MODE_INSTANCE)
		return PropertyInfo(Variant::OBJECT,"instance");
	else
		return PropertyInfo();

}

PropertyInfo VisualScriptYieldSignal::get_output_value_port_info(int p_idx) const{

	MethodInfo sr;

	if (!ObjectTypeDB::get_signal(_get_base_type(),signal,&sr))
		return PropertyInfo(); //no signal
	ERR_FAIL_INDEX_V(p_idx,sr.arguments.size(),PropertyInfo());
	return sr.arguments[p_idx];

}


String VisualScriptYieldSignal::get_caption() const {

	static const char*cname[3]= {
		"WaitSignal",
		"WaitNodeSignal",
		"WaitInstanceSigna;",
	};

	return cname[call_mode];
}

String VisualScriptYieldSignal::get_text() const {

	if (call_mode==CALL_MODE_SELF)
		return "  "+String(signal)+"()";
	else
		return "  "+_get_base_type()+"."+String(signal)+"()";

}


void VisualScriptYieldSignal::set_base_type(const StringName& p_type) {

	if (base_type==p_type)
		return;

	base_type=p_type;

	_change_notify();
	ports_changed_notify();
}

StringName VisualScriptYieldSignal::get_base_type() const{

	return base_type;
}

void VisualScriptYieldSignal::set_signal(const StringName& p_type){

	if (signal==p_type)
		return;

	signal=p_type;

	_change_notify();
	ports_changed_notify();
}
StringName VisualScriptYieldSignal::get_signal() const {


	return signal;
}

void VisualScriptYieldSignal::set_base_path(const NodePath& p_type) {

	if (base_path==p_type)
		return;

	base_path=p_type;

	_change_notify();
	ports_changed_notify();
}

NodePath VisualScriptYieldSignal::get_base_path() const {

	return base_path;
}


void VisualScriptYieldSignal::set_call_mode(CallMode p_mode) {

	if (call_mode==p_mode)
		return;

	call_mode=p_mode;

	_change_notify();
	ports_changed_notify();

}

VisualScriptYieldSignal::CallMode VisualScriptYieldSignal::get_call_mode() const {

	return call_mode;
}


void VisualScriptYieldSignal::_validate_property(PropertyInfo& property) const {

	if (property.name=="signal/base_type") {
		if (call_mode!=CALL_MODE_INSTANCE) {
			property.usage=PROPERTY_USAGE_NOEDITOR;
		}
	}


	if (property.name=="signal/node_path") {
		if (call_mode!=CALL_MODE_NODE_PATH) {
			property.usage=0;
		} else {

			Node *bnode = _get_base_node();
			if (bnode) {
				property.hint_string=bnode->get_path(); //convert to loong string
			} else {

			}
		}
	}

	if (property.name=="signal/signal") {
		property.hint=PROPERTY_HINT_ENUM;


		List<MethodInfo> methods;

		ObjectTypeDB::get_signal_list(_get_base_type(),&methods);

		List<String> mstring;
		for (List<MethodInfo>::Element *E=methods.front();E;E=E->next()) {
			if (E->get().name.begins_with("_"))
				continue;
			mstring.push_back(E->get().name.get_slice(":",0));
		}

		mstring.sort();

		String ml;
		for (List<String>::Element *E=mstring.front();E;E=E->next()) {

			if (ml!=String())
				ml+=",";
			ml+=E->get();
		}

		property.hint_string=ml;
	}


}


void VisualScriptYieldSignal::_bind_methods() {

	ObjectTypeDB::bind_method(_MD("set_base_type","base_type"),&VisualScriptYieldSignal::set_base_type);
	ObjectTypeDB::bind_method(_MD("get_base_type"),&VisualScriptYieldSignal::get_base_type);

	ObjectTypeDB::bind_method(_MD("set_signal","signal"),&VisualScriptYieldSignal::set_signal);
	ObjectTypeDB::bind_method(_MD("get_signal"),&VisualScriptYieldSignal::get_signal);

	ObjectTypeDB::bind_method(_MD("set_call_mode","mode"),&VisualScriptYieldSignal::set_call_mode);
	ObjectTypeDB::bind_method(_MD("get_call_mode"),&VisualScriptYieldSignal::get_call_mode);

	ObjectTypeDB::bind_method(_MD("set_base_path","base_path"),&VisualScriptYieldSignal::set_base_path);
	ObjectTypeDB::bind_method(_MD("get_base_path"),&VisualScriptYieldSignal::get_base_path);



	String bt;
	for(int i=0;i<Variant::VARIANT_MAX;i++) {
		if (i>0)
			bt+=",";

		bt+=Variant::get_type_name(Variant::Type(i));
	}

	ADD_PROPERTY(PropertyInfo(Variant::INT,"signal/call_mode",PROPERTY_HINT_ENUM,"Self,Node Path,Instance",PROPERTY_USAGE_NOEDITOR),_SCS("set_call_mode"),_SCS("get_call_mode"));
	ADD_PROPERTY(PropertyInfo(Variant::STRING,"signal/base_type",PROPERTY_HINT_TYPE_STRING,"Object"),_SCS("set_base_type"),_SCS("get_base_type"));
	ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH,"signal/node_path",PROPERTY_HINT_NODE_PATH_TO_EDITED_NODE),_SCS("set_base_path"),_SCS("get_base_path"));
	ADD_PROPERTY(PropertyInfo(Variant::STRING,"signal/signal"),_SCS("set_signal"),_SCS("get_signal"));


	BIND_CONSTANT( CALL_MODE_SELF );
	BIND_CONSTANT( CALL_MODE_NODE_PATH);
	BIND_CONSTANT( CALL_MODE_INSTANCE);

}

class VisualScriptNodeInstanceYieldSignal : public VisualScriptNodeInstance {
public:


	VisualScriptYieldSignal::CallMode call_mode;
	NodePath node_path;
	int output_args;
	StringName signal;

	VisualScriptYieldSignal *node;
	VisualScriptInstance *instance;



	virtual int get_working_memory_size() const { return 1; }
	//virtual bool is_output_port_unsequenced(int p_idx) const { return false; }
	//virtual bool get_output_port_unsequenced(int p_idx,Variant* r_value,Variant* p_working_mem,String &r_error) const { return true; }

	virtual int step(const Variant** p_inputs,Variant** p_outputs,StartMode p_start_mode,Variant* p_working_mem,Variant::CallError& r_error,String& r_error_str) {

		if (p_start_mode==START_MODE_RESUME_YIELD) {
			return 0; //resuming yield
		} else {
			//yield

			Object * object;

			switch(call_mode) {

				case VisualScriptYieldSignal::CALL_MODE_SELF: {

					object=instance->get_owner_ptr();

				} break;
				case VisualScriptYieldSignal::CALL_MODE_NODE_PATH: {

					Node* node = instance->get_owner_ptr()->cast_to<Node>();
					if (!node) {
						r_error.error=Variant::CallError::CALL_ERROR_INVALID_METHOD;
						r_error_str="Base object is not a Node!";
						return 0;
					}

					Node* another = node->get_node(node_path);
					if (!node) {
						r_error.error=Variant::CallError::CALL_ERROR_INVALID_METHOD;
						r_error_str="Path does not lead Node!";
						return 0;
					}

					object=another;

				} break;
				case VisualScriptYieldSignal::CALL_MODE_INSTANCE: {

					object = *p_inputs[0];
					if (!object) {
						r_error.error=Variant::CallError::CALL_ERROR_INVALID_METHOD;
						r_error_str="Supplied instance input is null.";
						return 0;

					}

				} break;

			}

			Ref<VisualScriptFunctionState> state;
			state.instance();

			state->connect_to_signal(object,signal,Array());

			*p_working_mem=state;

			return STEP_YIELD_BIT;
		}


	}


};

VisualScriptNodeInstance* VisualScriptYieldSignal::instance(VisualScriptInstance* p_instance) {

	VisualScriptNodeInstanceYieldSignal * instance = memnew(VisualScriptNodeInstanceYieldSignal );
	instance->node=this;
	instance->instance=p_instance;
	instance->signal=signal;
	instance->call_mode=call_mode;
	instance->node_path=base_path;
	instance->output_args = get_output_value_port_count();
	return instance;
}
VisualScriptYieldSignal::VisualScriptYieldSignal() {

	call_mode=CALL_MODE_INSTANCE;
	base_type="Object";

}

template<VisualScriptYieldSignal::CallMode cmode>
static Ref<VisualScriptNode> create_yield_signal_node(const String& p_name) {

	Ref<VisualScriptYieldSignal> node;
	node.instance();
	node->set_call_mode(cmode);
	return node;
}

void register_visual_script_yield_nodes() {

	VisualScriptLanguage::singleton->add_register_func("functions/wait/wait_frame",create_yield_node<VisualScriptYield::YIELD_FRAME>);
	VisualScriptLanguage::singleton->add_register_func("functions/wait/wait_fixed_frame",create_yield_node<VisualScriptYield::YIELD_FIXED_FRAME>);
	VisualScriptLanguage::singleton->add_register_func("functions/wait/wait_time",create_yield_node<VisualScriptYield::YIELD_WAIT>);

	VisualScriptLanguage::singleton->add_register_func("functions/yield/instance_signal",create_yield_signal_node<VisualScriptYieldSignal::CALL_MODE_INSTANCE>);
	VisualScriptLanguage::singleton->add_register_func("functions/yield/self_signal",create_yield_signal_node<VisualScriptYieldSignal::CALL_MODE_SELF>);
	VisualScriptLanguage::singleton->add_register_func("functions/yield/node_signal",create_yield_signal_node<VisualScriptYieldSignal::CALL_MODE_NODE_PATH>);

}
