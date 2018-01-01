/*************************************************************************/
/*  gdscript_compiler.cpp                                                */
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
#include "gdscript_compiler.h"

#include "gdscript.h"

bool GDScriptCompiler::_is_class_member_property(CodeGen &codegen, const StringName &p_name) {

	if (!codegen.function_node || codegen.function_node->_static)
		return false;

	return _is_class_member_property(codegen.script, p_name);
}

bool GDScriptCompiler::_is_class_member_property(GDScript *owner, const StringName &p_name) {

	GDScript *scr = owner;
	GDScriptNativeClass *nc = NULL;
	while (scr) {

		if (scr->native.is_valid())
			nc = scr->native.ptr();
		scr = scr->_base;
	}

	ERR_FAIL_COND_V(!nc, false);

	return ClassDB::has_property(nc->get_name(), p_name);
}

void GDScriptCompiler::_set_error(const String &p_error, const GDScriptParser::Node *p_node) {

	if (error != "")
		return;

	error = p_error;
	if (p_node) {
		err_line = p_node->line;
		err_column = p_node->column;
	} else {
		err_line = 0;
		err_column = 0;
	}
}

bool GDScriptCompiler::_create_unary_operator(CodeGen &codegen, const GDScriptParser::OperatorNode *on, Variant::Operator op, int p_stack_level) {

	ERR_FAIL_COND_V(on->arguments.size() != 1, false);

	int src_address_a = _parse_expression(codegen, on->arguments[0], p_stack_level);
	if (src_address_a < 0)
		return false;

	codegen.opcodes.push_back(GDScriptFunction::OPCODE_OPERATOR); // perform operator
	codegen.opcodes.push_back(op); //which operator
	codegen.opcodes.push_back(src_address_a); // argument 1
	codegen.opcodes.push_back(src_address_a); // argument 2 (repeated)
	//codegen.opcodes.push_back(GDScriptFunction::ADDR_TYPE_NIL); // argument 2 (unary only takes one parameter)
	return true;
}

bool GDScriptCompiler::_create_binary_operator(CodeGen &codegen, const GDScriptParser::OperatorNode *on, Variant::Operator op, int p_stack_level, bool p_initializer) {

	ERR_FAIL_COND_V(on->arguments.size() != 2, false);

	int src_address_a = _parse_expression(codegen, on->arguments[0], p_stack_level, false, p_initializer);
	if (src_address_a < 0)
		return false;
	if (src_address_a & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS)
		p_stack_level++; //uses stack for return, increase stack

	int src_address_b = _parse_expression(codegen, on->arguments[1], p_stack_level, false, p_initializer);
	if (src_address_b < 0)
		return false;

	codegen.opcodes.push_back(GDScriptFunction::OPCODE_OPERATOR); // perform operator
	codegen.opcodes.push_back(op); //which operator
	codegen.opcodes.push_back(src_address_a); // argument 1
	codegen.opcodes.push_back(src_address_b); // argument 2 (unary only takes one parameter)
	return true;
}

/*
int GDScriptCompiler::_parse_subexpression(CodeGen& codegen,const GDScriptParser::Node *p_expression) {


	int ret = _parse_expression(codegen,p_expression);
	if (ret<0)
		return ret;

	if (ret&(GDScriptFunction::ADDR_TYPE_STACK<<GDScriptFunction::ADDR_BITS)) {
		codegen.stack_level++;
		codegen.check_max_stack_level();
		//stack was used, keep value
	}

	return ret;
}
*/

int GDScriptCompiler::_parse_assign_right_expression(CodeGen &codegen, const GDScriptParser::OperatorNode *p_expression, int p_stack_level) {

	Variant::Operator var_op = Variant::OP_MAX;

	switch (p_expression->op) {

		case GDScriptParser::OperatorNode::OP_ASSIGN_ADD: var_op = Variant::OP_ADD; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_SUB: var_op = Variant::OP_SUBTRACT; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_MUL: var_op = Variant::OP_MULTIPLY; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_DIV: var_op = Variant::OP_DIVIDE; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_MOD: var_op = Variant::OP_MODULE; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_SHIFT_LEFT: var_op = Variant::OP_SHIFT_LEFT; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_SHIFT_RIGHT: var_op = Variant::OP_SHIFT_RIGHT; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_BIT_AND: var_op = Variant::OP_BIT_AND; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_BIT_OR: var_op = Variant::OP_BIT_OR; break;
		case GDScriptParser::OperatorNode::OP_ASSIGN_BIT_XOR: var_op = Variant::OP_BIT_XOR; break;
		case GDScriptParser::OperatorNode::OP_INIT_ASSIGN:
		case GDScriptParser::OperatorNode::OP_ASSIGN: {

			//none
		} break;
		default: {

			ERR_FAIL_V(-1);
		}
	}

	bool initializer = p_expression->op == GDScriptParser::OperatorNode::OP_INIT_ASSIGN;

	if (var_op == Variant::OP_MAX) {

		return _parse_expression(codegen, p_expression->arguments[1], p_stack_level, false, initializer);
	}

	if (!_create_binary_operator(codegen, p_expression, var_op, p_stack_level, initializer))
		return -1;

	int dst_addr = (p_stack_level) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
	codegen.opcodes.push_back(dst_addr); // append the stack level as destination address of the opcode
	codegen.alloc_stack(p_stack_level);
	return dst_addr;
}

int GDScriptCompiler::_parse_expression(CodeGen &codegen, const GDScriptParser::Node *p_expression, int p_stack_level, bool p_root, bool p_initializer) {

	switch (p_expression->type) {
		//should parse variable declaration and adjust stack accordingly...
		case GDScriptParser::Node::TYPE_IDENTIFIER: {
			//return identifier
			//wait, identifier could be a local variable or something else... careful here, must reference properly
			//as stack may be more interesting to work with

			//This could be made much simpler by just indexing "self", but done this way (with custom self-addressing modes) increases peformance a lot.

			const GDScriptParser::IdentifierNode *in = static_cast<const GDScriptParser::IdentifierNode *>(p_expression);

			StringName identifier = in->name;

			if (_is_class_member_property(codegen, identifier)) {
				//get property
				codegen.opcodes.push_back(GDScriptFunction::OPCODE_GET_MEMBER); // perform operator
				codegen.opcodes.push_back(codegen.get_name_map_pos(identifier)); // argument 2 (unary only takes one parameter)
				int dst_addr = (p_stack_level) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
				codegen.opcodes.push_back(dst_addr); // append the stack level as destination address of the opcode
				codegen.alloc_stack(p_stack_level);
				return dst_addr;
			}

			// TRY STACK!
			if (!p_initializer && codegen.stack_identifiers.has(identifier)) {

				int pos = codegen.stack_identifiers[identifier];
				return pos | (GDScriptFunction::ADDR_TYPE_STACK_VARIABLE << GDScriptFunction::ADDR_BITS);
			}
			//TRY MEMBERS!
			if (!codegen.function_node || !codegen.function_node->_static) {

				// TRY MEMBER VARIABLES!
				//static function
				if (codegen.script->member_indices.has(identifier)) {

					int idx = codegen.script->member_indices[identifier].index;
					return idx | (GDScriptFunction::ADDR_TYPE_MEMBER << GDScriptFunction::ADDR_BITS); //argument (stack root)
				}
			}

			//TRY CLASS CONSTANTS

			GDScript *owner = codegen.script;
			while (owner) {

				GDScript *scr = owner;
				GDScriptNativeClass *nc = NULL;
				while (scr) {

					if (scr->constants.has(identifier)) {

						//int idx=scr->constants[identifier];
						int idx = codegen.get_name_map_pos(identifier);
						return idx | (GDScriptFunction::ADDR_TYPE_CLASS_CONSTANT << GDScriptFunction::ADDR_BITS); //argument (stack root)
					}
					if (scr->native.is_valid())
						nc = scr->native.ptr();
					scr = scr->_base;
				}

				// CLASS C++ Integer Constant

				if (nc) {

					bool success = false;
					int constant = ClassDB::get_integer_constant(nc->get_name(), identifier, &success);
					if (success) {
						Variant key = constant;
						int idx;

						if (!codegen.constant_map.has(key)) {

							idx = codegen.constant_map.size();
							codegen.constant_map[key] = idx;

						} else {
							idx = codegen.constant_map[key];
						}

						return idx | (GDScriptFunction::ADDR_TYPE_LOCAL_CONSTANT << GDScriptFunction::ADDR_BITS); //make it a local constant (faster access)
					}
				}

				owner = owner->_owner;
			}

			/*
			 handled in constants now
			 if (codegen.script->subclasses.has(identifier)) {
				//same with a subclass, make it a local constant.
				int idx = codegen.get_constant_pos(codegen.script->subclasses[identifier]);
				return idx|(GDScriptFunction::ADDR_TYPE_LOCAL_CONSTANT<<GDScriptFunction::ADDR_BITS); //make it a local constant (faster access)

			}*/

			if (GDScriptLanguage::get_singleton()->get_global_map().has(identifier)) {

				int idx = GDScriptLanguage::get_singleton()->get_global_map()[identifier];
				return idx | (GDScriptFunction::ADDR_TYPE_GLOBAL << GDScriptFunction::ADDR_BITS); //argument (stack root)
			}

			//not found, error

			_set_error("Identifier not found: " + String(identifier), p_expression);

			return -1;

		} break;
		case GDScriptParser::Node::TYPE_CONSTANT: {
			//return constant
			const GDScriptParser::ConstantNode *cn = static_cast<const GDScriptParser::ConstantNode *>(p_expression);

			int idx;

			if (!codegen.constant_map.has(cn->value)) {

				idx = codegen.constant_map.size();
				codegen.constant_map[cn->value] = idx;

			} else {
				idx = codegen.constant_map[cn->value];
			}

			return idx | (GDScriptFunction::ADDR_TYPE_LOCAL_CONSTANT << GDScriptFunction::ADDR_BITS); //argument (stack root)

		} break;
		case GDScriptParser::Node::TYPE_SELF: {
			//return constant
			if (codegen.function_node && codegen.function_node->_static) {
				_set_error("'self' not present in static function!", p_expression);
				return -1;
			}
			return (GDScriptFunction::ADDR_TYPE_SELF << GDScriptFunction::ADDR_BITS);
		} break;
		case GDScriptParser::Node::TYPE_ARRAY: {

			const GDScriptParser::ArrayNode *an = static_cast<const GDScriptParser::ArrayNode *>(p_expression);
			Vector<int> values;

			int slevel = p_stack_level;

			for (int i = 0; i < an->elements.size(); i++) {

				int ret = _parse_expression(codegen, an->elements[i], slevel);
				if (ret < 0)
					return ret;
				if (ret & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
					slevel++;
					codegen.alloc_stack(slevel);
				}

				values.push_back(ret);
			}

			codegen.opcodes.push_back(GDScriptFunction::OPCODE_CONSTRUCT_ARRAY);
			codegen.opcodes.push_back(values.size());
			for (int i = 0; i < values.size(); i++)
				codegen.opcodes.push_back(values[i]);

			int dst_addr = (p_stack_level) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
			codegen.opcodes.push_back(dst_addr); // append the stack level as destination address of the opcode
			codegen.alloc_stack(p_stack_level);
			return dst_addr;

		} break;
		case GDScriptParser::Node::TYPE_DICTIONARY: {

			const GDScriptParser::DictionaryNode *dn = static_cast<const GDScriptParser::DictionaryNode *>(p_expression);
			Vector<int> values;

			int slevel = p_stack_level;

			for (int i = 0; i < dn->elements.size(); i++) {

				int ret = _parse_expression(codegen, dn->elements[i].key, slevel);
				if (ret < 0)
					return ret;
				if (ret & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
					slevel++;
					codegen.alloc_stack(slevel);
				}

				values.push_back(ret);

				ret = _parse_expression(codegen, dn->elements[i].value, slevel);
				if (ret < 0)
					return ret;
				if (ret & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
					slevel++;
					codegen.alloc_stack(slevel);
				}

				values.push_back(ret);
			}

			codegen.opcodes.push_back(GDScriptFunction::OPCODE_CONSTRUCT_DICTIONARY);
			codegen.opcodes.push_back(dn->elements.size());
			for (int i = 0; i < values.size(); i++)
				codegen.opcodes.push_back(values[i]);

			int dst_addr = (p_stack_level) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
			codegen.opcodes.push_back(dst_addr); // append the stack level as destination address of the opcode
			codegen.alloc_stack(p_stack_level);
			return dst_addr;

		} break;
		case GDScriptParser::Node::TYPE_OPERATOR: {
			//hell breaks loose

			const GDScriptParser::OperatorNode *on = static_cast<const GDScriptParser::OperatorNode *>(p_expression);
			switch (on->op) {

				//call/constructor operator
				case GDScriptParser::OperatorNode::OP_PARENT_CALL: {

					ERR_FAIL_COND_V(on->arguments.size() < 1, -1);

					const GDScriptParser::IdentifierNode *in = (const GDScriptParser::IdentifierNode *)on->arguments[0];

					Vector<int> arguments;
					int slevel = p_stack_level;
					for (int i = 1; i < on->arguments.size(); i++) {

						int ret = _parse_expression(codegen, on->arguments[i], slevel);
						if (ret < 0)
							return ret;
						if (ret & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
							slevel++;
							codegen.alloc_stack(slevel);
						}
						arguments.push_back(ret);
					}

					//push call bytecode
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_CALL_SELF_BASE); // basic type constructor

					codegen.opcodes.push_back(codegen.get_name_map_pos(in->name)); //instance
					codegen.opcodes.push_back(arguments.size()); //argument count
					codegen.alloc_call(arguments.size());
					for (int i = 0; i < arguments.size(); i++)
						codegen.opcodes.push_back(arguments[i]); //arguments

				} break;
				case GDScriptParser::OperatorNode::OP_CALL: {

					if (on->arguments[0]->type == GDScriptParser::Node::TYPE_TYPE) {
						//construct a basic type
						ERR_FAIL_COND_V(on->arguments.size() < 1, -1);

						const GDScriptParser::TypeNode *tn = (const GDScriptParser::TypeNode *)on->arguments[0];
						int vtype = tn->vtype;

						Vector<int> arguments;
						int slevel = p_stack_level;
						for (int i = 1; i < on->arguments.size(); i++) {

							int ret = _parse_expression(codegen, on->arguments[i], slevel);
							if (ret < 0)
								return ret;
							if (ret & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
								slevel++;
								codegen.alloc_stack(slevel);
							}
							arguments.push_back(ret);
						}

						//push call bytecode
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_CONSTRUCT); // basic type constructor
						codegen.opcodes.push_back(vtype); //instance
						codegen.opcodes.push_back(arguments.size()); //argument count
						codegen.alloc_call(arguments.size());
						for (int i = 0; i < arguments.size(); i++)
							codegen.opcodes.push_back(arguments[i]); //arguments

					} else if (on->arguments[0]->type == GDScriptParser::Node::TYPE_BUILT_IN_FUNCTION) {
						//built in function

						ERR_FAIL_COND_V(on->arguments.size() < 1, -1);

						Vector<int> arguments;
						int slevel = p_stack_level;
						for (int i = 1; i < on->arguments.size(); i++) {

							int ret = _parse_expression(codegen, on->arguments[i], slevel);
							if (ret < 0)
								return ret;

							if (ret & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
								slevel++;
								codegen.alloc_stack(slevel);
							}

							arguments.push_back(ret);
						}

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_CALL_BUILT_IN);
						codegen.opcodes.push_back(static_cast<const GDScriptParser::BuiltInFunctionNode *>(on->arguments[0])->function);
						codegen.opcodes.push_back(on->arguments.size() - 1);
						codegen.alloc_call(on->arguments.size() - 1);
						for (int i = 0; i < arguments.size(); i++)
							codegen.opcodes.push_back(arguments[i]);

					} else {
						//regular function
						ERR_FAIL_COND_V(on->arguments.size() < 2, -1);

						const GDScriptParser::Node *instance = on->arguments[0];

						if (instance->type == GDScriptParser::Node::TYPE_SELF) {
							//room for optimization
						}

						Vector<int> arguments;
						int slevel = p_stack_level;

						for (int i = 0; i < on->arguments.size(); i++) {

							int ret;

							if (i == 0 && on->arguments[i]->type == GDScriptParser::Node::TYPE_SELF && codegen.function_node && codegen.function_node->_static) {
								//static call to self
								ret = (GDScriptFunction::ADDR_TYPE_CLASS << GDScriptFunction::ADDR_BITS);
							} else if (i == 1) {

								if (on->arguments[i]->type != GDScriptParser::Node::TYPE_IDENTIFIER) {
									_set_error("Attempt to call a non-identifier.", on);
									return -1;
								}
								GDScriptParser::IdentifierNode *id = static_cast<GDScriptParser::IdentifierNode *>(on->arguments[i]);
								ret = codegen.get_name_map_pos(id->name);

							} else {

								ret = _parse_expression(codegen, on->arguments[i], slevel);
								if (ret < 0)
									return ret;
								if (ret & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
									slevel++;
									codegen.alloc_stack(slevel);
								}
							}
							arguments.push_back(ret);
						}

						codegen.opcodes.push_back(p_root ? GDScriptFunction::OPCODE_CALL : GDScriptFunction::OPCODE_CALL_RETURN); // perform operator
						codegen.opcodes.push_back(on->arguments.size() - 2);
						codegen.alloc_call(on->arguments.size() - 2);
						for (int i = 0; i < arguments.size(); i++)
							codegen.opcodes.push_back(arguments[i]);
					}
				} break;
				case GDScriptParser::OperatorNode::OP_YIELD: {

					ERR_FAIL_COND_V(on->arguments.size() && on->arguments.size() != 2, -1);

					Vector<int> arguments;
					int slevel = p_stack_level;
					for (int i = 0; i < on->arguments.size(); i++) {

						int ret = _parse_expression(codegen, on->arguments[i], slevel);
						if (ret < 0)
							return ret;
						if (ret & (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS)) {
							slevel++;
							codegen.alloc_stack(slevel);
						}
						arguments.push_back(ret);
					}

					//push call bytecode
					codegen.opcodes.push_back(arguments.size() == 0 ? GDScriptFunction::OPCODE_YIELD : GDScriptFunction::OPCODE_YIELD_SIGNAL); // basic type constructor
					for (int i = 0; i < arguments.size(); i++)
						codegen.opcodes.push_back(arguments[i]); //arguments
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_YIELD_RESUME);
					//next will be where to place the result :)

				} break;

				//indexing operator
				case GDScriptParser::OperatorNode::OP_INDEX:
				case GDScriptParser::OperatorNode::OP_INDEX_NAMED: {

					ERR_FAIL_COND_V(on->arguments.size() != 2, -1);

					int slevel = p_stack_level;
					bool named = (on->op == GDScriptParser::OperatorNode::OP_INDEX_NAMED);

					int from = _parse_expression(codegen, on->arguments[0], slevel);
					if (from < 0)
						return from;

					int index;
					if (named) {
						if (on->arguments[0]->type == GDScriptParser::Node::TYPE_SELF && codegen.script && codegen.function_node && !codegen.function_node->_static) {

							GDScriptParser::IdentifierNode *identifier = static_cast<GDScriptParser::IdentifierNode *>(on->arguments[1]);
							const Map<StringName, GDScript::MemberInfo>::Element *MI = codegen.script->member_indices.find(identifier->name);

#ifdef DEBUG_ENABLED
							if (MI && MI->get().getter == codegen.function_node->name) {
								String n = static_cast<GDScriptParser::IdentifierNode *>(on->arguments[1])->name;
								_set_error("Must use '" + n + "' instead of 'self." + n + "' in getter.", on);
								return -1;
							}
#endif

							if (MI && MI->get().getter == "") {
								// Faster than indexing self (as if no self. had been used)
								return (MI->get().index) | (GDScriptFunction::ADDR_TYPE_MEMBER << GDScriptFunction::ADDR_BITS);
							}
						}

						index = codegen.get_name_map_pos(static_cast<GDScriptParser::IdentifierNode *>(on->arguments[1])->name);

					} else {

						if (on->arguments[1]->type == GDScriptParser::Node::TYPE_CONSTANT && static_cast<const GDScriptParser::ConstantNode *>(on->arguments[1])->value.get_type() == Variant::STRING) {
							//also, somehow, named (speed up anyway)
							StringName name = static_cast<const GDScriptParser::ConstantNode *>(on->arguments[1])->value;
							index = codegen.get_name_map_pos(name);
							named = true;

						} else {
							//regular indexing
							if (from & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
								slevel++;
								codegen.alloc_stack(slevel);
							}

							index = _parse_expression(codegen, on->arguments[1], slevel);
							if (index < 0)
								return index;
						}
					}

					codegen.opcodes.push_back(named ? GDScriptFunction::OPCODE_GET_NAMED : GDScriptFunction::OPCODE_GET); // perform operator
					codegen.opcodes.push_back(from); // argument 1
					codegen.opcodes.push_back(index); // argument 2 (unary only takes one parameter)

				} break;
				case GDScriptParser::OperatorNode::OP_AND: {

					// AND operator with early out on failure

					int res = _parse_expression(codegen, on->arguments[0], p_stack_level);
					if (res < 0)
						return res;
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF_NOT);
					codegen.opcodes.push_back(res);
					int jump_fail_pos = codegen.opcodes.size();
					codegen.opcodes.push_back(0);

					res = _parse_expression(codegen, on->arguments[1], p_stack_level);
					if (res < 0)
						return res;

					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF_NOT);
					codegen.opcodes.push_back(res);
					int jump_fail_pos2 = codegen.opcodes.size();
					codegen.opcodes.push_back(0);

					codegen.alloc_stack(p_stack_level); //it will be used..
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN_TRUE);
					codegen.opcodes.push_back(p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
					codegen.opcodes.push_back(codegen.opcodes.size() + 3);
					codegen.opcodes[jump_fail_pos] = codegen.opcodes.size();
					codegen.opcodes[jump_fail_pos2] = codegen.opcodes.size();
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN_FALSE);
					codegen.opcodes.push_back(p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
					return p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS;

				} break;
				case GDScriptParser::OperatorNode::OP_OR: {

					// OR operator with early out on success

					int res = _parse_expression(codegen, on->arguments[0], p_stack_level);
					if (res < 0)
						return res;
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF);
					codegen.opcodes.push_back(res);
					int jump_success_pos = codegen.opcodes.size();
					codegen.opcodes.push_back(0);

					res = _parse_expression(codegen, on->arguments[1], p_stack_level);
					if (res < 0)
						return res;

					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF);
					codegen.opcodes.push_back(res);
					int jump_success_pos2 = codegen.opcodes.size();
					codegen.opcodes.push_back(0);

					codegen.alloc_stack(p_stack_level); //it will be used..
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN_FALSE);
					codegen.opcodes.push_back(p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
					codegen.opcodes.push_back(codegen.opcodes.size() + 3);
					codegen.opcodes[jump_success_pos] = codegen.opcodes.size();
					codegen.opcodes[jump_success_pos2] = codegen.opcodes.size();
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN_TRUE);
					codegen.opcodes.push_back(p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
					return p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS;

				} break;
				// ternary operators
				case GDScriptParser::OperatorNode::OP_TERNARY_IF: {

					// x IF a ELSE y operator with early out on failure

					int res = _parse_expression(codegen, on->arguments[0], p_stack_level);
					if (res < 0)
						return res;
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF_NOT);
					codegen.opcodes.push_back(res);
					int jump_fail_pos = codegen.opcodes.size();
					codegen.opcodes.push_back(0);

					res = _parse_expression(codegen, on->arguments[1], p_stack_level);
					if (res < 0)
						return res;

					codegen.alloc_stack(p_stack_level); //it will be used..
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN);
					codegen.opcodes.push_back(p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
					codegen.opcodes.push_back(res);
					codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
					int jump_past_pos = codegen.opcodes.size();
					codegen.opcodes.push_back(0);

					codegen.opcodes[jump_fail_pos] = codegen.opcodes.size();
					res = _parse_expression(codegen, on->arguments[2], p_stack_level);
					if (res < 0)
						return res;

					codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN);
					codegen.opcodes.push_back(p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
					codegen.opcodes.push_back(res);

					codegen.opcodes[jump_past_pos] = codegen.opcodes.size();

					return p_stack_level | GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS;

				} break;
				//unary operators
				case GDScriptParser::OperatorNode::OP_NEG: {
					if (!_create_unary_operator(codegen, on, Variant::OP_NEGATE, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_NOT: {
					if (!_create_unary_operator(codegen, on, Variant::OP_NOT, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_BIT_INVERT: {
					if (!_create_unary_operator(codegen, on, Variant::OP_BIT_NEGATE, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_PREINC: {
				} break; //?
				case GDScriptParser::OperatorNode::OP_PREDEC: {
				} break;
				case GDScriptParser::OperatorNode::OP_INC: {
				} break;
				case GDScriptParser::OperatorNode::OP_DEC: {
				} break;
				//binary operators (in precedence order)
				case GDScriptParser::OperatorNode::OP_IN: {
					if (!_create_binary_operator(codegen, on, Variant::OP_IN, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_EQUAL: {
					if (!_create_binary_operator(codegen, on, Variant::OP_EQUAL, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_NOT_EQUAL: {
					if (!_create_binary_operator(codegen, on, Variant::OP_NOT_EQUAL, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_LESS: {
					if (!_create_binary_operator(codegen, on, Variant::OP_LESS, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_LESS_EQUAL: {
					if (!_create_binary_operator(codegen, on, Variant::OP_LESS_EQUAL, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_GREATER: {
					if (!_create_binary_operator(codegen, on, Variant::OP_GREATER, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_GREATER_EQUAL: {
					if (!_create_binary_operator(codegen, on, Variant::OP_GREATER_EQUAL, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_ADD: {
					if (!_create_binary_operator(codegen, on, Variant::OP_ADD, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_SUB: {
					if (!_create_binary_operator(codegen, on, Variant::OP_SUBTRACT, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_MUL: {
					if (!_create_binary_operator(codegen, on, Variant::OP_MULTIPLY, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_DIV: {
					if (!_create_binary_operator(codegen, on, Variant::OP_DIVIDE, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_MOD: {
					if (!_create_binary_operator(codegen, on, Variant::OP_MODULE, p_stack_level)) return -1;
				} break;
				//case GDScriptParser::OperatorNode::OP_SHIFT_LEFT: { if (!_create_binary_operator(codegen,on,Variant::OP_SHIFT_LEFT,p_stack_level)) return -1;} break;
				//case GDScriptParser::OperatorNode::OP_SHIFT_RIGHT: { if (!_create_binary_operator(codegen,on,Variant::OP_SHIFT_RIGHT,p_stack_level)) return -1;} break;
				case GDScriptParser::OperatorNode::OP_BIT_AND: {
					if (!_create_binary_operator(codegen, on, Variant::OP_BIT_AND, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_BIT_OR: {
					if (!_create_binary_operator(codegen, on, Variant::OP_BIT_OR, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_BIT_XOR: {
					if (!_create_binary_operator(codegen, on, Variant::OP_BIT_XOR, p_stack_level)) return -1;
				} break;
				//shift
				case GDScriptParser::OperatorNode::OP_SHIFT_LEFT: {
					if (!_create_binary_operator(codegen, on, Variant::OP_SHIFT_LEFT, p_stack_level)) return -1;
				} break;
				case GDScriptParser::OperatorNode::OP_SHIFT_RIGHT: {
					if (!_create_binary_operator(codegen, on, Variant::OP_SHIFT_RIGHT, p_stack_level)) return -1;
				} break;
				//assignment operators
				case GDScriptParser::OperatorNode::OP_ASSIGN_ADD:
				case GDScriptParser::OperatorNode::OP_ASSIGN_SUB:
				case GDScriptParser::OperatorNode::OP_ASSIGN_MUL:
				case GDScriptParser::OperatorNode::OP_ASSIGN_DIV:
				case GDScriptParser::OperatorNode::OP_ASSIGN_MOD:
				case GDScriptParser::OperatorNode::OP_ASSIGN_SHIFT_LEFT:
				case GDScriptParser::OperatorNode::OP_ASSIGN_SHIFT_RIGHT:
				case GDScriptParser::OperatorNode::OP_ASSIGN_BIT_AND:
				case GDScriptParser::OperatorNode::OP_ASSIGN_BIT_OR:
				case GDScriptParser::OperatorNode::OP_ASSIGN_BIT_XOR:
				case GDScriptParser::OperatorNode::OP_INIT_ASSIGN:
				case GDScriptParser::OperatorNode::OP_ASSIGN: {

					ERR_FAIL_COND_V(on->arguments.size() != 2, -1);

					if (on->arguments[0]->type == GDScriptParser::Node::TYPE_OPERATOR && (static_cast<GDScriptParser::OperatorNode *>(on->arguments[0])->op == GDScriptParser::OperatorNode::OP_INDEX || static_cast<GDScriptParser::OperatorNode *>(on->arguments[0])->op == GDScriptParser::OperatorNode::OP_INDEX_NAMED)) {

					// SET (chained) MODE!
#ifdef DEBUG_ENABLED
						if (static_cast<GDScriptParser::OperatorNode *>(on->arguments[0])->op == GDScriptParser::OperatorNode::OP_INDEX_NAMED) {
							const GDScriptParser::OperatorNode *inon = static_cast<GDScriptParser::OperatorNode *>(on->arguments[0]);

							if (inon->arguments[0]->type == GDScriptParser::Node::TYPE_SELF && codegen.script && codegen.function_node && !codegen.function_node->_static) {

								const Map<StringName, GDScript::MemberInfo>::Element *MI = codegen.script->member_indices.find(static_cast<GDScriptParser::IdentifierNode *>(inon->arguments[1])->name);
								if (MI && MI->get().setter == codegen.function_node->name) {
									String n = static_cast<GDScriptParser::IdentifierNode *>(inon->arguments[1])->name;
									_set_error("Must use '" + n + "' instead of 'self." + n + "' in setter.", inon);
									return -1;
								}
							}
						}
#endif

						int slevel = p_stack_level;

						GDScriptParser::OperatorNode *op = static_cast<GDScriptParser::OperatorNode *>(on->arguments[0]);

						/* Find chain of sets */

						StringName assign_property;

						List<GDScriptParser::OperatorNode *> chain;

						{
							//create get/set chain
							GDScriptParser::OperatorNode *n = op;
							while (true) {

								chain.push_back(n);
								if (n->arguments[0]->type != GDScriptParser::Node::TYPE_OPERATOR) {

									//check for a built-in property
									if (n->arguments[0]->type == GDScriptParser::Node::TYPE_IDENTIFIER) {

										GDScriptParser::IdentifierNode *identifier = static_cast<GDScriptParser::IdentifierNode *>(n->arguments[0]);
										if (_is_class_member_property(codegen, identifier->name)) {
											assign_property = identifier->name;
										}
									}
									break;
								}
								n = static_cast<GDScriptParser::OperatorNode *>(n->arguments[0]);
								if (n->op != GDScriptParser::OperatorNode::OP_INDEX && n->op != GDScriptParser::OperatorNode::OP_INDEX_NAMED)
									break;
							}
						}

						/* Chain of gets */

						//get at (potential) root stack pos, so it can be returned
						int prev_pos = _parse_expression(codegen, chain.back()->get()->arguments[0], slevel);
						if (prev_pos < 0)
							return prev_pos;
						int retval = prev_pos;

						//print_line("retval: "+itos(retval));

						if (retval & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
							slevel++;
							codegen.alloc_stack(slevel);
						}

						Vector<int> setchain;

						if (assign_property != StringName()) {

							// recover and assign at the end, this allows stuff like
							// position.x+=2.0
							// in Node2D
							setchain.push_back(prev_pos);
							setchain.push_back(codegen.get_name_map_pos(assign_property));
							setchain.push_back(GDScriptFunction::OPCODE_SET_MEMBER);
						}

						for (List<GDScriptParser::OperatorNode *>::Element *E = chain.back(); E; E = E->prev()) {

							if (E == chain.front()) //ignore first
								break;

							bool named = E->get()->op == GDScriptParser::OperatorNode::OP_INDEX_NAMED;
							int key_idx;

							if (named) {

								key_idx = codegen.get_name_map_pos(static_cast<const GDScriptParser::IdentifierNode *>(E->get()->arguments[1])->name);
								//printf("named key %x\n",key_idx);

							} else {

								if (prev_pos & (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS)) {
									slevel++;
									codegen.alloc_stack(slevel);
								}

								GDScriptParser::Node *key = E->get()->arguments[1];
								key_idx = _parse_expression(codegen, key, slevel);
								//printf("expr key %x\n",key_idx);

								//stack was raised here if retval was stack but..
							}

							if (key_idx < 0) //error
								return key_idx;

							codegen.opcodes.push_back(named ? GDScriptFunction::OPCODE_GET_NAMED : GDScriptFunction::OPCODE_GET);
							codegen.opcodes.push_back(prev_pos);
							codegen.opcodes.push_back(key_idx);
							slevel++;
							codegen.alloc_stack(slevel);
							int dst_pos = (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) | slevel;

							codegen.opcodes.push_back(dst_pos);

							//add in reverse order, since it will be reverted

							setchain.push_back(dst_pos);
							setchain.push_back(key_idx);
							setchain.push_back(prev_pos);
							setchain.push_back(named ? GDScriptFunction::OPCODE_SET_NAMED : GDScriptFunction::OPCODE_SET);

							prev_pos = dst_pos;
						}

						setchain.invert();

						int set_index;
						bool named = false;

						if (static_cast<const GDScriptParser::OperatorNode *>(op)->op == GDScriptParser::OperatorNode::OP_INDEX_NAMED) {

							set_index = codegen.get_name_map_pos(static_cast<const GDScriptParser::IdentifierNode *>(op->arguments[1])->name);
							named = true;
						} else {

							set_index = _parse_expression(codegen, op->arguments[1], slevel + 1);
							named = false;
						}

						if (set_index < 0) //error
							return set_index;

						if (set_index & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
							slevel++;
							codegen.alloc_stack(slevel);
						}

						int set_value = _parse_assign_right_expression(codegen, on, slevel + 1);
						if (set_value < 0) //error
							return set_value;

						codegen.opcodes.push_back(named ? GDScriptFunction::OPCODE_SET_NAMED : GDScriptFunction::OPCODE_SET);
						codegen.opcodes.push_back(prev_pos);
						codegen.opcodes.push_back(set_index);
						codegen.opcodes.push_back(set_value);

						for (int i = 0; i < setchain.size(); i++) {

							codegen.opcodes.push_back(setchain[i]);
						}

						return retval;

					} else if (on->arguments[0]->type == GDScriptParser::Node::TYPE_IDENTIFIER && _is_class_member_property(codegen, static_cast<GDScriptParser::IdentifierNode *>(on->arguments[0])->name)) {
						//assignment to member property

						int slevel = p_stack_level;

						int src_address = _parse_assign_right_expression(codegen, on, slevel);
						if (src_address < 0)
							return -1;

						StringName name = static_cast<GDScriptParser::IdentifierNode *>(on->arguments[0])->name;

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_SET_MEMBER);
						codegen.opcodes.push_back(codegen.get_name_map_pos(name));
						codegen.opcodes.push_back(src_address);

						return GDScriptFunction::ADDR_TYPE_NIL << GDScriptFunction::ADDR_BITS;
					} else {

						//REGULAR ASSIGNMENT MODE!!

						int slevel = p_stack_level;

						int dst_address_a = _parse_expression(codegen, on->arguments[0], slevel, false, on->op == GDScriptParser::OperatorNode::OP_INIT_ASSIGN);
						if (dst_address_a < 0)
							return -1;

						if (dst_address_a & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) {
							slevel++;
							codegen.alloc_stack(slevel);
						}

						int src_address_b = _parse_assign_right_expression(codegen, on, slevel);
						if (src_address_b < 0)
							return -1;

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN); // perform operator
						codegen.opcodes.push_back(dst_address_a); // argument 1
						codegen.opcodes.push_back(src_address_b); // argument 2 (unary only takes one parameter)
						return dst_address_a; //if anything, returns wathever was assigned or correct stack position
					}

				} break;
				case GDScriptParser::OperatorNode::OP_IS: {

					ERR_FAIL_COND_V(on->arguments.size() != 2, false);

					int slevel = p_stack_level;

					int src_address_a = _parse_expression(codegen, on->arguments[0], slevel);
					if (src_address_a < 0)
						return -1;

					if (src_address_a & GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS)
						slevel++; //uses stack for return, increase stack

					int src_address_b = _parse_expression(codegen, on->arguments[1], slevel);
					if (src_address_b < 0)
						return -1;

					codegen.opcodes.push_back(GDScriptFunction::OPCODE_EXTENDS_TEST); // perform operator
					codegen.opcodes.push_back(src_address_a); // argument 1
					codegen.opcodes.push_back(src_address_b); // argument 2 (unary only takes one parameter)

				} break;
				default: {

					ERR_EXPLAIN("Bug in bytecode compiler, unexpected operator #" + itos(on->op) + " in parse tree while parsing expression.");
					ERR_FAIL_V(0); //unreachable code

				} break;
			}

			int dst_addr = (p_stack_level) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
			codegen.opcodes.push_back(dst_addr); // append the stack level as destination address of the opcode
			codegen.alloc_stack(p_stack_level);
			return dst_addr;
		} break;
		//TYPE_TYPE,
		default: {

			ERR_EXPLAIN("Bug in bytecode compiler, unexpected node in parse tree while parsing expression.");
			ERR_FAIL_V(-1); //unreachable code
		} break;
	}

	ERR_FAIL_V(-1); //unreachable code
}

Error GDScriptCompiler::_parse_block(CodeGen &codegen, const GDScriptParser::BlockNode *p_block, int p_stack_level, int p_break_addr, int p_continue_addr) {

	codegen.push_stack_identifiers();
	int new_identifiers = 0;
	codegen.current_line = p_block->line;

	for (int i = 0; i < p_block->statements.size(); i++) {

		const GDScriptParser::Node *s = p_block->statements[i];

		switch (s->type) {
			case GDScriptParser::Node::TYPE_NEWLINE: {
#ifdef DEBUG_ENABLED
				const GDScriptParser::NewLineNode *nl = static_cast<const GDScriptParser::NewLineNode *>(s);
				codegen.opcodes.push_back(GDScriptFunction::OPCODE_LINE);
				codegen.opcodes.push_back(nl->line);
				codegen.current_line = nl->line;
#endif
			} break;
			case GDScriptParser::Node::TYPE_CONTROL_FLOW: {
				// try subblocks

				const GDScriptParser::ControlFlowNode *cf = static_cast<const GDScriptParser::ControlFlowNode *>(s);

				switch (cf->cf_type) {

					case GDScriptParser::ControlFlowNode::CF_MATCH: {
						GDScriptParser::MatchNode *match = cf->match;

						GDScriptParser::IdentifierNode *id = memnew(GDScriptParser::IdentifierNode);
						id->name = "#match_value";

						// var #match_value
						// copied because there is no _parse_statement :(
						codegen.add_stack_identifier(id->name, p_stack_level++);
						codegen.alloc_stack(p_stack_level);
						new_identifiers++;

						GDScriptParser::OperatorNode *op = memnew(GDScriptParser::OperatorNode);
						op->op = GDScriptParser::OperatorNode::OP_ASSIGN;
						op->arguments.push_back(id);
						op->arguments.push_back(match->val_to_match);

						int ret = _parse_expression(codegen, op, p_stack_level);
						if (ret < 0) {
							return ERR_PARSE_ERROR;
						}

						// break address
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(codegen.opcodes.size() + 3);
						int break_addr = codegen.opcodes.size();
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(0); // break addr

						for (int j = 0; j < match->compiled_pattern_branches.size(); j++) {
							GDScriptParser::MatchNode::CompiledPatternBranch branch = match->compiled_pattern_branches[j];

							// jump over continue
							// jump unconditionally
							// continue address
							// compile the condition
							int ret = _parse_expression(codegen, branch.compiled_pattern, p_stack_level);
							if (ret < 0) {
								return ERR_PARSE_ERROR;
							}

							codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF);
							codegen.opcodes.push_back(ret);
							codegen.opcodes.push_back(codegen.opcodes.size() + 3);
							int continue_addr = codegen.opcodes.size();
							codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
							codegen.opcodes.push_back(0);

							Error err = _parse_block(codegen, branch.body, p_stack_level, p_break_addr, continue_addr);
							if (err) {
								return ERR_PARSE_ERROR;
							}

							codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
							codegen.opcodes.push_back(break_addr);

							codegen.opcodes[continue_addr + 1] = codegen.opcodes.size();
						}

						codegen.opcodes[break_addr + 1] = codegen.opcodes.size();

					} break;

					case GDScriptParser::ControlFlowNode::CF_IF: {

#ifdef DEBUG_ENABLED
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_LINE);
						codegen.opcodes.push_back(cf->line);
						codegen.current_line = cf->line;
#endif
						int ret = _parse_expression(codegen, cf->arguments[0], p_stack_level, false);
						if (ret < 0)
							return ERR_PARSE_ERROR;

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF_NOT);
						codegen.opcodes.push_back(ret);
						int else_addr = codegen.opcodes.size();
						codegen.opcodes.push_back(0); //temporary

						Error err = _parse_block(codegen, cf->body, p_stack_level, p_break_addr, p_continue_addr);
						if (err)
							return err;

						if (cf->body_else) {

							codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
							int end_addr = codegen.opcodes.size();
							codegen.opcodes.push_back(0);
							codegen.opcodes[else_addr] = codegen.opcodes.size();

							Error err = _parse_block(codegen, cf->body_else, p_stack_level, p_break_addr, p_continue_addr);
							if (err)
								return err;

							codegen.opcodes[end_addr] = codegen.opcodes.size();
						} else {
							//end without else
							codegen.opcodes[else_addr] = codegen.opcodes.size();
						}

					} break;
					case GDScriptParser::ControlFlowNode::CF_FOR: {

						int slevel = p_stack_level;
						int iter_stack_pos = slevel;
						int iterator_pos = (slevel++) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
						int counter_pos = (slevel++) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
						int container_pos = (slevel++) | (GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS);
						codegen.alloc_stack(slevel);

						codegen.push_stack_identifiers();
						codegen.add_stack_identifier(static_cast<const GDScriptParser::IdentifierNode *>(cf->arguments[0])->name, iter_stack_pos);

						int ret = _parse_expression(codegen, cf->arguments[1], slevel, false);
						if (ret < 0)
							return ERR_COMPILATION_FAILED;

						//assign container
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSIGN);
						codegen.opcodes.push_back(container_pos);
						codegen.opcodes.push_back(ret);

						//begin loop
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_ITERATE_BEGIN);
						codegen.opcodes.push_back(counter_pos);
						codegen.opcodes.push_back(container_pos);
						codegen.opcodes.push_back(codegen.opcodes.size() + 4);
						codegen.opcodes.push_back(iterator_pos);
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP); //skip code for next
						codegen.opcodes.push_back(codegen.opcodes.size() + 8);
						//break loop
						int break_pos = codegen.opcodes.size();
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP); //skip code for next
						codegen.opcodes.push_back(0); //skip code for next
						//next loop
						int continue_pos = codegen.opcodes.size();
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_ITERATE);
						codegen.opcodes.push_back(counter_pos);
						codegen.opcodes.push_back(container_pos);
						codegen.opcodes.push_back(break_pos);
						codegen.opcodes.push_back(iterator_pos);

						Error err = _parse_block(codegen, cf->body, slevel, break_pos, continue_pos);
						if (err)
							return err;

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(continue_pos);
						codegen.opcodes[break_pos + 1] = codegen.opcodes.size();

						codegen.pop_stack_identifiers();

					} break;
					case GDScriptParser::ControlFlowNode::CF_WHILE: {

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(codegen.opcodes.size() + 3);
						int break_addr = codegen.opcodes.size();
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(0);
						int continue_addr = codegen.opcodes.size();

						int ret = _parse_expression(codegen, cf->arguments[0], p_stack_level, false);
						if (ret < 0)
							return ERR_PARSE_ERROR;
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_IF_NOT);
						codegen.opcodes.push_back(ret);
						codegen.opcodes.push_back(break_addr);
						Error err = _parse_block(codegen, cf->body, p_stack_level, break_addr, continue_addr);
						if (err)
							return err;
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(continue_addr);

						codegen.opcodes[break_addr + 1] = codegen.opcodes.size();

					} break;
					case GDScriptParser::ControlFlowNode::CF_SWITCH: {

					} break;
					case GDScriptParser::ControlFlowNode::CF_BREAK: {

						if (p_break_addr < 0) {

							_set_error("'break'' not within loop", cf);
							return ERR_COMPILATION_FAILED;
						}
						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(p_break_addr);

					} break;
					case GDScriptParser::ControlFlowNode::CF_CONTINUE: {

						if (p_continue_addr < 0) {

							_set_error("'continue' not within loop", cf);
							return ERR_COMPILATION_FAILED;
						}

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP);
						codegen.opcodes.push_back(p_continue_addr);

					} break;
					case GDScriptParser::ControlFlowNode::CF_RETURN: {

						int ret;

						if (cf->arguments.size()) {

							ret = _parse_expression(codegen, cf->arguments[0], p_stack_level, false);
							if (ret < 0)
								return ERR_PARSE_ERROR;

						} else {

							ret = GDScriptFunction::ADDR_TYPE_NIL << GDScriptFunction::ADDR_BITS;
						}

						codegen.opcodes.push_back(GDScriptFunction::OPCODE_RETURN);
						codegen.opcodes.push_back(ret);

					} break;
				}
			} break;
			case GDScriptParser::Node::TYPE_ASSERT: {
#ifdef DEBUG_ENABLED
				// try subblocks

				const GDScriptParser::AssertNode *as = static_cast<const GDScriptParser::AssertNode *>(s);

				int ret = _parse_expression(codegen, as->condition, p_stack_level, false);
				if (ret < 0)
					return ERR_PARSE_ERROR;

				codegen.opcodes.push_back(GDScriptFunction::OPCODE_ASSERT);
				codegen.opcodes.push_back(ret);
#endif
			} break;
			case GDScriptParser::Node::TYPE_BREAKPOINT: {
#ifdef DEBUG_ENABLED
				// try subblocks
				codegen.opcodes.push_back(GDScriptFunction::OPCODE_BREAKPOINT);
#endif
			} break;
			case GDScriptParser::Node::TYPE_LOCAL_VAR: {

				const GDScriptParser::LocalVarNode *lv = static_cast<const GDScriptParser::LocalVarNode *>(s);

				if (_is_class_member_property(codegen, lv->name)) {
					_set_error("Name for local variable '" + String(lv->name) + "' can't shadow class property of the same name.", lv);
					return ERR_ALREADY_EXISTS;
				}

				codegen.add_stack_identifier(lv->name, p_stack_level++);
				codegen.alloc_stack(p_stack_level);
				new_identifiers++;

			} break;
			default: {
				//expression
				int ret = _parse_expression(codegen, s, p_stack_level, true);
				if (ret < 0)
					return ERR_PARSE_ERROR;
			} break;
		}
	}
	codegen.pop_stack_identifiers();
	return OK;
}

Error GDScriptCompiler::_parse_function(GDScript *p_script, const GDScriptParser::ClassNode *p_class, const GDScriptParser::FunctionNode *p_func, bool p_for_ready) {

	Vector<int> bytecode;
	CodeGen codegen;

	codegen.class_node = p_class;
	codegen.script = p_script;
	codegen.function_node = p_func;
	codegen.stack_max = 0;
	codegen.current_line = 0;
	codegen.call_max = 0;
	codegen.debug_stack = ScriptDebugger::get_singleton() != NULL;
	Vector<StringName> argnames;

	int stack_level = 0;

	if (p_func) {
		for (int i = 0; i < p_func->arguments.size(); i++) {
			if (_is_class_member_property(p_script, p_func->arguments[i])) {
				_set_error("Name for argument '" + String(p_func->arguments[i]) + "' can't shadow class property of the same name.", p_func);
				return ERR_ALREADY_EXISTS;
			}
			codegen.add_stack_identifier(p_func->arguments[i], i);
#ifdef TOOLS_ENABLED
			argnames.push_back(p_func->arguments[i]);
#endif
		}
		stack_level = p_func->arguments.size();
	}

	codegen.alloc_stack(stack_level);

	/* Parse initializer -if applies- */

	bool is_initializer = !p_for_ready && !p_func;

	if (is_initializer || (p_func && String(p_func->name) == "_init")) {
		//parse initializer for class members
		if (!p_func && p_class->extends_used && p_script->native.is_null()) {

			//call implicit parent constructor
			codegen.opcodes.push_back(GDScriptFunction::OPCODE_CALL_SELF_BASE);
			codegen.opcodes.push_back(codegen.get_name_map_pos("_init"));
			codegen.opcodes.push_back(0);
			codegen.opcodes.push_back((GDScriptFunction::ADDR_TYPE_STACK << GDScriptFunction::ADDR_BITS) | 0);
		}
		Error err = _parse_block(codegen, p_class->initializer, stack_level);
		if (err)
			return err;
		is_initializer = true;
	}

	if (p_for_ready || (p_func && String(p_func->name) == "_ready")) {
		//parse initializer for class members
		if (p_class->ready->statements.size()) {
			Error err = _parse_block(codegen, p_class->ready, stack_level);
			if (err)
				return err;
		}
	}

	/* Parse default argument code -if applies- */

	Vector<int> defarg_addr;
	StringName func_name;

	if (p_func) {

		if (p_func->default_values.size()) {

			codegen.opcodes.push_back(GDScriptFunction::OPCODE_JUMP_TO_DEF_ARGUMENT);
			defarg_addr.push_back(codegen.opcodes.size());
			for (int i = 0; i < p_func->default_values.size(); i++) {

				_parse_expression(codegen, p_func->default_values[i], stack_level, true);
				defarg_addr.push_back(codegen.opcodes.size());
			}

			defarg_addr.invert();
		}

		Error err = _parse_block(codegen, p_func->body, stack_level);
		if (err)
			return err;

		func_name = p_func->name;
	} else {
		if (p_for_ready)
			func_name = "_ready";
		else
			func_name = "_init";
	}

	codegen.opcodes.push_back(GDScriptFunction::OPCODE_END);

	/*
	if (String(p_func->name)=="") { //initializer func
		gdfunc = &p_script->initializer;
	*/
	//} else { //regular func
	p_script->member_functions[func_name] = memnew(GDScriptFunction);
	GDScriptFunction *gdfunc = p_script->member_functions[func_name];
	//}

	if (p_func) {
		gdfunc->_static = p_func->_static;
		gdfunc->rpc_mode = p_func->rpc_mode;
	}

#ifdef TOOLS_ENABLED
	gdfunc->arg_names = argnames;
#endif
	//constants
	if (codegen.constant_map.size()) {
		gdfunc->_constant_count = codegen.constant_map.size();
		gdfunc->constants.resize(codegen.constant_map.size());
		gdfunc->_constants_ptr = &gdfunc->constants[0];
		const Variant *K = NULL;
		while ((K = codegen.constant_map.next(K))) {
			int idx = codegen.constant_map[*K];
			gdfunc->constants[idx] = *K;
		}
	} else {

		gdfunc->_constants_ptr = NULL;
		gdfunc->_constant_count = 0;
	}
	//global names
	if (codegen.name_map.size()) {

		gdfunc->global_names.resize(codegen.name_map.size());
		gdfunc->_global_names_ptr = &gdfunc->global_names[0];
		for (Map<StringName, int>::Element *E = codegen.name_map.front(); E; E = E->next()) {

			gdfunc->global_names[E->get()] = E->key();
		}
		gdfunc->_global_names_count = gdfunc->global_names.size();

	} else {
		gdfunc->_global_names_ptr = NULL;
		gdfunc->_global_names_count = 0;
	}

	if (codegen.opcodes.size()) {

		gdfunc->code = codegen.opcodes;
		gdfunc->_code_ptr = &gdfunc->code[0];
		gdfunc->_code_size = codegen.opcodes.size();

	} else {

		gdfunc->_code_ptr = NULL;
		gdfunc->_code_size = 0;
	}

	if (defarg_addr.size()) {

		gdfunc->default_arguments = defarg_addr;
		gdfunc->_default_arg_count = defarg_addr.size() - 1;
		gdfunc->_default_arg_ptr = &gdfunc->default_arguments[0];
	} else {
		gdfunc->_default_arg_count = 0;
		gdfunc->_default_arg_ptr = NULL;
	}

	gdfunc->_argument_count = p_func ? p_func->arguments.size() : 0;
	gdfunc->_stack_size = codegen.stack_max;
	gdfunc->_call_size = codegen.call_max;
	gdfunc->name = func_name;
#ifdef DEBUG_ENABLED
	if (ScriptDebugger::get_singleton()) {
		String signature;
		//path
		if (p_script->get_path() != String())
			signature += p_script->get_path();
		//loc
		if (p_func) {
			signature += "::" + itos(p_func->body->line);
		} else {
			signature += "::0";
		}

		//function and class

		if (p_class->name) {
			signature += "::" + String(p_class->name) + "." + String(func_name);
		} else {
			signature += "::" + String(func_name);
		}

		gdfunc->profile.signature = signature;
	}
#endif
	gdfunc->_script = p_script;
	gdfunc->source = source;

#ifdef DEBUG_ENABLED

	{
		gdfunc->func_cname = (String(source) + " - " + String(func_name)).utf8();
		gdfunc->_func_cname = gdfunc->func_cname.get_data();
	}

#endif
	if (p_func) {
		gdfunc->_initial_line = p_func->line;
#ifdef TOOLS_ENABLED

		p_script->member_lines[func_name] = p_func->line;
#endif
	} else {
		gdfunc->_initial_line = 0;
	}

	if (codegen.debug_stack)
		gdfunc->stack_debug = codegen.stack_debug;

	if (is_initializer)
		p_script->initializer = gdfunc;

	return OK;
}

Error GDScriptCompiler::_parse_class(GDScript *p_script, GDScript *p_owner, const GDScriptParser::ClassNode *p_class, bool p_keep_state) {

	Map<StringName, Ref<GDScript> > old_subclasses;

	if (p_keep_state) {
		old_subclasses = p_script->subclasses;
	}

	p_script->native = Ref<GDScriptNativeClass>();
	p_script->base = Ref<GDScript>();
	p_script->_base = NULL;
	p_script->members.clear();
	p_script->constants.clear();
	for (Map<StringName, GDScriptFunction *>::Element *E = p_script->member_functions.front(); E; E = E->next()) {
		memdelete(E->get());
	}
	p_script->member_functions.clear();
	p_script->member_indices.clear();
	p_script->member_info.clear();
	p_script->_signals.clear();
	p_script->initializer = NULL;

	p_script->subclasses.clear();
	p_script->_owner = p_owner;
	p_script->tool = p_class->tool;
	p_script->name = p_class->name;

	Ref<GDScriptNativeClass> native;

	if (p_class->extends_used) {
		//do inheritance
		String path = p_class->extends_file;

		Ref<GDScript> script;

		if (path != "") {
			//path (and optionally subclasses)

			if (path.is_rel_path()) {

				String base;

				if (p_owner) {
					GDScript *current_class = p_owner;
					while (current_class != NULL) {
						base = current_class->get_path();
						if (base == "")
							current_class = current_class->_owner;
						else
							break;
					}
				} else {
					base = p_script->get_path();
				}

				if (base == "" || base.is_rel_path()) {
					_set_error("Could not resolve relative path for parent class: " + path, p_class);
					return ERR_FILE_NOT_FOUND;
				}
				path = base.get_base_dir().plus_file(path).simplify_path();
			}
			script = ResourceLoader::load(path);
			if (script.is_null()) {
				_set_error("Could not load base class: " + path, p_class);
				return ERR_FILE_NOT_FOUND;
			}
			if (!script->valid) {

				_set_error("Script not fully loaded (cyclic preload?): " + path, p_class);
				return ERR_BUSY;
			}
			//print_line("EXTENDS PATH: "+path+" script is "+itos(script.is_valid())+" indices is "+itos(script->member_indices.size())+" valid? "+itos(script->valid));

			if (p_class->extends_class.size()) {

				for (int i = 0; i < p_class->extends_class.size(); i++) {

					String sub = p_class->extends_class[i];
					if (script->subclasses.has(sub)) {

						Ref<Script> subclass = script->subclasses[sub]; //avoid reference from disappearing
						script = subclass;
					} else {

						_set_error("Could not find subclass: " + sub, p_class);
						return ERR_FILE_NOT_FOUND;
					}
				}
			}

		} else {

			ERR_FAIL_COND_V(p_class->extends_class.size() == 0, ERR_BUG);
			//look around for the subclasses

			String base = p_class->extends_class[0];
			GDScript *p = p_owner;
			Ref<GDScript> base_class;

			while (p) {

				if (p->subclasses.has(base)) {

					base_class = p->subclasses[base];
					break;
				}

				if (p->constants.has(base)) {

					base_class = p->constants[base];
					if (base_class.is_null()) {
						_set_error("Constant is not a class: " + base, p_class);
						return ERR_SCRIPT_FAILED;
					}
					break;
				}

				p = p->_owner;
			}

			if (base_class.is_valid()) {

				String ident = base;

				for (int i = 1; i < p_class->extends_class.size(); i++) {

					String subclass = p_class->extends_class[i];

					ident += ("." + subclass);

					if (base_class->subclasses.has(subclass)) {

						base_class = base_class->subclasses[subclass];
					} else if (base_class->constants.has(subclass)) {

						Ref<GDScript> new_base_class = base_class->constants[subclass];
						if (new_base_class.is_null()) {
							_set_error("Constant is not a class: " + ident, p_class);
							return ERR_SCRIPT_FAILED;
						}
						base_class = new_base_class;
					} else {

						_set_error("Could not find subclass: " + ident, p_class);
						return ERR_FILE_NOT_FOUND;
					}
				}

				script = base_class;

			} else {

				if (p_class->extends_class.size() > 1) {

					_set_error("Invalid inheritance (unknown class+subclasses)", p_class);
					return ERR_FILE_NOT_FOUND;
				}
				//if not found, try engine classes
				if (!GDScriptLanguage::get_singleton()->get_global_map().has(base)) {

					_set_error("Unknown class: '" + base + "'", p_class);
					return ERR_FILE_NOT_FOUND;
				}

				int base_idx = GDScriptLanguage::get_singleton()->get_global_map()[base];
				native = GDScriptLanguage::get_singleton()->get_global_array()[base_idx];
				if (!native.is_valid()) {

					_set_error("Global not a class: '" + base + "'", p_class);

					return ERR_FILE_NOT_FOUND;
				}
			}
		}

		if (script.is_valid()) {

			p_script->base = script;
			p_script->_base = p_script->base.ptr();
			p_script->member_indices = script->member_indices;

		} else if (native.is_valid()) {

			p_script->native = native;
		} else {

			_set_error("Could not determine inheritance", p_class);
			return ERR_FILE_NOT_FOUND;
		}

	} else {
		// without extends, implicitly extend Reference
		int native_idx = GDScriptLanguage::get_singleton()->get_global_map()["Reference"];
		native = GDScriptLanguage::get_singleton()->get_global_array()[native_idx];
		ERR_FAIL_COND_V(native.is_null(), ERR_BUG);
		p_script->native = native;
	}

	//print_line("Script: "+p_script->get_path()+" indices: "+itos(p_script->member_indices.size()));

	for (int i = 0; i < p_class->variables.size(); i++) {

		StringName name = p_class->variables[i].identifier;
		if (p_script->member_indices.has(name)) {
			_set_error("Member '" + name + "' already exists (in current or parent class)", p_class);
			return ERR_ALREADY_EXISTS;
		}
		if (_is_class_member_property(p_script, name)) {
			_set_error("Member '" + name + "' already exists as a class property.", p_class);
			return ERR_ALREADY_EXISTS;
		}

		if (p_class->variables[i]._export.type != Variant::NIL) {

			p_script->member_info[name] = p_class->variables[i]._export;
#ifdef TOOLS_ENABLED
			if (p_class->variables[i].default_value.get_type() != Variant::NIL) {

				p_script->member_default_values[name] = p_class->variables[i].default_value;
			}
#endif
		} else {

			p_script->member_info[name] = PropertyInfo(Variant::NIL, name, PROPERTY_HINT_NONE, "", PROPERTY_USAGE_SCRIPT_VARIABLE);
		}

		//int new_idx = p_script->member_indices.size();
		GDScript::MemberInfo minfo;
		minfo.index = p_script->member_indices.size();
		minfo.setter = p_class->variables[i].setter;
		minfo.getter = p_class->variables[i].getter;
		minfo.rpc_mode = p_class->variables[i].rpc_mode;

		p_script->member_indices[name] = minfo;
		p_script->members.insert(name);

#ifdef TOOLS_ENABLED

		p_script->member_lines[name] = p_class->variables[i].line;
#endif
	}

	for (int i = 0; i < p_class->constant_expressions.size(); i++) {

		StringName name = p_class->constant_expressions[i].identifier;
		ERR_CONTINUE(p_class->constant_expressions[i].expression->type != GDScriptParser::Node::TYPE_CONSTANT);

		if (_is_class_member_property(p_script, name)) {
			_set_error("Member '" + name + "' already exists as a class property.", p_class);
			return ERR_ALREADY_EXISTS;
		}

		GDScriptParser::ConstantNode *constant = static_cast<GDScriptParser::ConstantNode *>(p_class->constant_expressions[i].expression);

		p_script->constants.insert(name, constant->value);
//p_script->constants[constant->value].make_const();
#ifdef TOOLS_ENABLED

		p_script->member_lines[name] = p_class->constant_expressions[i].expression->line;
#endif
	}

	for (int i = 0; i < p_class->_signals.size(); i++) {

		StringName name = p_class->_signals[i].name;

		GDScript *c = p_script;

		while (c) {

			if (c->_signals.has(name)) {
				_set_error("Signal '" + name + "' redefined (in current or parent class)", p_class);
				return ERR_ALREADY_EXISTS;
			}

			if (c->base.is_valid()) {
				c = c->base.ptr();
			} else {
				c = NULL;
			}
		}

		if (native.is_valid()) {
			if (ClassDB::has_signal(native->get_name(), name)) {
				_set_error("Signal '" + name + "' redefined (original in native class '" + String(native->get_name()) + "')", p_class);
				return ERR_ALREADY_EXISTS;
			}
		}

		p_script->_signals[name] = p_class->_signals[i].arguments;
	}
	//parse sub-classes

	for (int i = 0; i < p_class->subclasses.size(); i++) {
		StringName name = p_class->subclasses[i]->name;

		Ref<GDScript> subclass;

		if (old_subclasses.has(name)) {
			subclass = old_subclasses[name];
		} else {
			subclass.instance();
		}

		Error err = _parse_class(subclass.ptr(), p_script, p_class->subclasses[i], p_keep_state);
		if (err)
			return err;

#ifdef TOOLS_ENABLED

		p_script->member_lines[name] = p_class->subclasses[i]->line;
#endif

		p_script->constants.insert(name, subclass); //once parsed, goes to the list of constants
		p_script->subclasses.insert(name, subclass);
	}

	//parse methods

	bool has_initializer = false;
	bool has_ready = false;

	for (int i = 0; i < p_class->functions.size(); i++) {

		if (!has_initializer && p_class->functions[i]->name == "_init")
			has_initializer = true;
		if (!has_ready && p_class->functions[i]->name == "_ready")
			has_ready = true;
		Error err = _parse_function(p_script, p_class, p_class->functions[i]);
		if (err)
			return err;
	}

	//parse static methods

	for (int i = 0; i < p_class->static_functions.size(); i++) {

		Error err = _parse_function(p_script, p_class, p_class->static_functions[i]);
		if (err)
			return err;
	}

	if (!has_initializer) {
		//create a constructor
		Error err = _parse_function(p_script, p_class, NULL);
		if (err)
			return err;
	}

	if (!has_ready && p_class->ready->statements.size()) {
		//create a constructor
		Error err = _parse_function(p_script, p_class, NULL, true);
		if (err)
			return err;
	}

#ifdef DEBUG_ENABLED
	//validate setters/getters if debug is enabled
	for (int i = 0; i < p_class->variables.size(); i++) {

		if (p_class->variables[i].setter) {
			const Map<StringName, GDScriptFunction *>::Element *E = p_script->get_member_functions().find(p_class->variables[i].setter);
			if (!E) {
				_set_error("Setter function '" + String(p_class->variables[i].setter) + "' not found in class.", NULL);
				err_line = p_class->variables[i].line;
				err_column = 0;
				return ERR_PARSE_ERROR;
			}

			if (E->get()->is_static()) {

				_set_error("Setter function '" + String(p_class->variables[i].setter) + "' is static.", NULL);
				err_line = p_class->variables[i].line;
				err_column = 0;
				return ERR_PARSE_ERROR;
			}
		}
		if (p_class->variables[i].getter) {
			const Map<StringName, GDScriptFunction *>::Element *E = p_script->get_member_functions().find(p_class->variables[i].getter);
			if (!E) {
				_set_error("Getter function '" + String(p_class->variables[i].getter) + "' not found in class.", NULL);
				err_line = p_class->variables[i].line;
				err_column = 0;
				return ERR_PARSE_ERROR;
			}

			if (E->get()->is_static()) {

				_set_error("Getter function '" + String(p_class->variables[i].getter) + "' is static.", NULL);
				err_line = p_class->variables[i].line;
				err_column = 0;
				return ERR_PARSE_ERROR;
			}
		}
	}

	//validate instances if keeping state

	if (p_keep_state) {

		print_line("RELOAD KEEP " + p_script->path);
		for (Set<Object *>::Element *E = p_script->instances.front(); E;) {

			Set<Object *>::Element *N = E->next();

			ScriptInstance *si = E->get()->get_script_instance();
			if (si->is_placeholder()) {
#ifdef TOOLS_ENABLED
				PlaceHolderScriptInstance *psi = static_cast<PlaceHolderScriptInstance *>(si);

				if (p_script->is_tool()) {
					//re-create as an instance
					p_script->placeholders.erase(psi); //remove placeholder

					GDScriptInstance *instance = memnew(GDScriptInstance);
					instance->base_ref = Object::cast_to<Reference>(E->get());
					instance->members.resize(p_script->member_indices.size());
					instance->script = Ref<GDScript>(p_script);
					instance->owner = E->get();

					//needed for hot reloading
					for (Map<StringName, GDScript::MemberInfo>::Element *E = p_script->member_indices.front(); E; E = E->next()) {
						instance->member_indices_cache[E->key()] = E->get().index;
					}
					instance->owner->set_script_instance(instance);

					/* STEP 2, INITIALIZE AND CONSRTUCT */

					Variant::CallError ce;
					p_script->initializer->call(instance, NULL, 0, ce);

					if (ce.error != Variant::CallError::CALL_OK) {
						//well, tough luck, not goinna do anything here
					}
				}
#endif
			} else {

				GDScriptInstance *gi = static_cast<GDScriptInstance *>(si);
				gi->reload_members();
			}

			E = N;
		}
	}
#endif

	p_script->valid = true;
	return OK;
}

Error GDScriptCompiler::compile(const GDScriptParser *p_parser, GDScript *p_script, bool p_keep_state) {

	err_line = -1;
	err_column = -1;
	error = "";
	parser = p_parser;
	const GDScriptParser::Node *root = parser->get_parse_tree();
	ERR_FAIL_COND_V(root->type != GDScriptParser::Node::TYPE_CLASS, ERR_INVALID_DATA);

	source = p_script->get_path();

	Error err = _parse_class(p_script, NULL, static_cast<const GDScriptParser::ClassNode *>(root), p_keep_state);

	if (err)
		return err;

	return OK;
}

String GDScriptCompiler::get_error() const {

	return error;
}
int GDScriptCompiler::get_error_line() const {

	return err_line;
}
int GDScriptCompiler::get_error_column() const {

	return err_column;
}

GDScriptCompiler::GDScriptCompiler() {
}
