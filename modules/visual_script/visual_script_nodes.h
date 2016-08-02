#ifndef VISUAL_SCRIPT_NODES_H
#define VISUAL_SCRIPT_NODES_H

#include "visual_script.h"

class VisualScriptFunction : public VisualScriptNode {

	OBJ_TYPE(VisualScriptFunction,VisualScriptNode)


	struct Argument {
		String name;
		Variant::Type type;
	};

	Vector<Argument> arguments;
protected:

	bool _set(const StringName& p_name, const Variant& p_value);
	bool _get(const StringName& p_name,Variant &r_ret) const;
	void _get_property_list( List<PropertyInfo> *p_list) const;

public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void add_argument(Variant::Type p_type,const String& p_name,int p_index=-1);
	void set_argument_type(int p_argidx,Variant::Type p_type);
	Variant::Type get_argument_type(int p_argidx) const;
	void set_argument_name(int p_argidx,const String& p_name);
	String get_argument_name(int p_argidx) const;
	void remove_argument(int p_argidx);
	int get_argument_count() const;

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptFunction();
};


class VisualScriptOperator : public VisualScriptNode {

	OBJ_TYPE(VisualScriptOperator,VisualScriptNode)


	Variant::Operator op;
protected:

	static void _bind_methods();
public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void set_operator(Variant::Operator p_op);
	Variant::Operator get_operator() const;

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptOperator();
};


class VisualScriptVariable : public VisualScriptNode {

	OBJ_TYPE(VisualScriptVariable,VisualScriptNode)


	StringName variable;
protected:

	virtual void _validate_property(PropertyInfo& property) const;
	static void _bind_methods();
public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void set_variable(StringName p_var);
	StringName get_variable() const;

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptVariable();
};

class VisualScriptConstant : public VisualScriptNode {

	OBJ_TYPE(VisualScriptConstant,VisualScriptNode)


	Variant::Type type;
	Variant value;
protected:
	virtual void _validate_property(PropertyInfo& property) const;
	static void _bind_methods();

public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void set_constant_type(Variant::Type p_type);
	Variant::Type get_constant_type() const;

	void set_constant_value(Variant p_value);
	Variant get_constant_value() const;

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptConstant();
};


class VisualScriptIndexGet : public VisualScriptNode {

	OBJ_TYPE(VisualScriptIndexGet,VisualScriptNode)


public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptIndexGet();
};


class VisualScriptIndexSet : public VisualScriptNode {

	OBJ_TYPE(VisualScriptIndexSet,VisualScriptNode)


public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptIndexSet();
};



class VisualScriptGlobalConstant : public VisualScriptNode {

	OBJ_TYPE(VisualScriptGlobalConstant,VisualScriptNode)

	int index;

	static void _bind_methods();
public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void set_global_constant(int p_which);
	int get_global_constant();

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptGlobalConstant();
};



class VisualScriptMathConstant : public VisualScriptNode {

	OBJ_TYPE(VisualScriptMathConstant,VisualScriptNode)
public:

	enum MathConstant {
		MATH_CONSTANT_ONE,
		MATH_CONSTANT_PI,
		MATH_CONSTANT_2PI,
		MATH_CONSTANT_HALF_PI,
		MATH_CONSTANT_E,
		MATH_CONSTANT_SQRT2,
		MATH_CONSTANT_MAX,
	};

private:
	static const char* const_name[MATH_CONSTANT_MAX];
	MathConstant constant;
protected:
	static void _bind_methods();
public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void set_math_constant(MathConstant p_which);
	MathConstant get_math_constant();

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptMathConstant();
};

VARIANT_ENUM_CAST( VisualScriptMathConstant::MathConstant )

class VisualScriptSingleton : public VisualScriptNode {

	OBJ_TYPE(VisualScriptSingleton,VisualScriptNode)

	String singleton;

	static void _bind_methods();
public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void set_singleton(const String &p_string);
	String get_singleton();

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptSingleton();
};




class VisualScriptSceneNode : public VisualScriptNode {

	OBJ_TYPE(VisualScriptSceneNode,VisualScriptNode)

	NodePath path;
protected:
	virtual void _validate_property(PropertyInfo& property) const;
	static void _bind_methods();
public:

	virtual int get_output_sequence_port_count() const;
	virtual bool has_input_sequence_port() const;


	virtual String get_output_sequence_port_text(int p_port) const;


	virtual int get_input_value_port_count() const;
	virtual int get_output_value_port_count() const;


	virtual PropertyInfo get_input_value_port_info(int p_idx) const;
	virtual PropertyInfo get_output_value_port_info(int p_idx) const;

	virtual String get_caption() const;
	virtual String get_text() const;

	void set_node_path(const NodePath &p_path);
	NodePath get_node_path();

	virtual VisualScriptNodeInstance* instance(VScriptInstance* p_instance);

	VisualScriptSceneNode();
};



void register_visual_script_nodes();

#endif // VISUAL_SCRIPT_NODES_H
