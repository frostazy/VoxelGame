#include "HyVoxelPrivatePCH.h"

// suppress _itoa warnings
#define _CRT_SECURE_NO_WARNINGS

#include "ModCompiler.h"
#include "Reduction.h"
#include "LSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>

using namespace std;
using namespace HyVoxel;

namespace GoldCPP {

enum ActionType
{
	ACTION_VAR,
	ACTION_SNAP,
	ACTION_MOVE,
	ACTION_ROTATE,
	ACTION_SCALE,
	ACTION_CENTER_PIVOT,
	ACTION_RESET_PIVOT,
	ACTION_ALIGN,
	ACTION_ROTATEX,
	ACTION_ROTATEY,
	ACTION_ROTATEZ,
	ACTION_OCCLUDE,
	ACTION_LOCAL,
	ACTION_BOX,
	ACTION_NGON,
	ACTION_LOFT_BOX,
	ACTION_LOFT_NGON,
	ACTION_PRISM,
	ACTION_ROOF,
	ACTION_CALL,
	ACTION_DEFER,
	ACTION_SELECT,
	ACTION_INSTANCE,
	ACTION_DIVIDE,
	ACTION_REPEAT,
	ACTION_PUSH,
	ACTION_POP,
	ACTION_MATERIAL
};

struct Action
{
	ActionType type;
	map<string, string> attributes;
	string filename;
	int line;
	int column;

	Action()
		: line(-1)
		, column(-1)
	{}
};

class Rule
{
public:
	list<Action>		actions;
	list<string>		conditionals;
	list<string>		occlusionScopes;
	list<string>		occlusionVolumeIds;
	list<string>		occlusionMargins;

	Rule() {}

	Action* snap(const string& dir, const string& id, const string& radius)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_SNAP;
		action.attributes["dir"] = dir;
		action.attributes["id"] = id;
		action.attributes["radius"] = radius;
		actions.push_back(action);
		return &action;
	}

	Action* var(const string& identifier, const string& value, bool constant)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_VAR;
		action.attributes["id"] = identifier;
		action.attributes["value"] = value;
		if (constant)
			action.attributes["const"] = "true";
		else
			action.attributes["const"] = "false";
		return &action;
	}

	Action* loft_box(const string& positions, const string& widths, const string& depths)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_LOFT_BOX;
		action.attributes["positions"] = positions;
		action.attributes["widths"] = widths;
		action.attributes["depths"] = depths;
		return &action;
	}

	Action* loft_ngon(const string& sides, const string& positions, const string& radii)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_LOFT_NGON;
		action.attributes["sides"] = sides;
		action.attributes["positions"] = positions;
		action.attributes["radii"] = radii;
		return &action;
	}

	Action* select(const string& subscope, const string& moduleName)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_SELECT;
		action.attributes["dimension"] = "face";
		action.attributes["type"] = subscope;
		action.attributes["module"] = moduleName;
		return &action;
	}

	Action* instance(const string& instance)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_INSTANCE;
		action.attributes["geometry"] = instance;
		return &action;
	}

	Action* repeat_x(const string& length, const string& module, const string& flex, const string& snap)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_REPEAT;
		action.attributes["dir"] = "x";
		action.attributes["length"] = length;
		action.attributes["snap"] = snap;
		action.attributes["module"] = module;
		action.attributes["flex"] = flex;
		return &action;
	}

	Action* repeat_y(const string& length, const string& module, const string& flex, const string& snap)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_REPEAT;
		action.attributes["dir"] = "y";
		action.attributes["length"] = length;
		action.attributes["snap"] = snap;
		action.attributes["module"] = module;
		action.attributes["flex"] = flex;
		return &action;
	}

	Action* repeat_z(const string& length, const string& module, const string& flex, const string& snap)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_REPEAT;
		action.attributes["dir"] = "z";
		action.attributes["length"] = length;
		action.attributes["snap"] = snap;
		action.attributes["module"] = module;
		action.attributes["flex"] = flex;
		return &action;
	}

	Action* subdivide_x(const string& pattern, const string& modules)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_DIVIDE;
		action.attributes["type"] = "x";
		action.attributes["pattern"] = pattern;
		action.attributes["modules"] = modules;
		return &action;
	}

	Action* subdivide_y(const string& pattern, const string& modules)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_DIVIDE;
		action.attributes["type"] = "y";
		action.attributes["pattern"] = pattern;
		action.attributes["modules"] = modules;
		return &action;
	}

	Action* subdivide_z(const string& pattern, const string& modules)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_DIVIDE;
		action.attributes["type"] = "z";
		action.attributes["pattern"] = pattern;
		action.attributes["modules"] = modules;
		return &action;
	}

	Action* subdivide_prism(const string& position, const string& center, const string& top, const string& left, const string& right)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_DIVIDE;
		action.attributes["type"] = "prism";
		action.attributes["position"] = position;
		action.attributes["center"] = center;
		action.attributes["top"] = top;
		action.attributes["left"] = left;
		action.attributes["right"] = right;
		return &action;
	}

	Action* module(const string& name)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CALL;
		action.attributes["name"] = name;
		action.attributes["generations"] = "-1";
		return &action;
	}

	Action* push()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_PUSH;
		return &action;
	}

	Action* pop()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_POP;
		return &action;
	}

	Action* defer(const string& name)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_DEFER;
		action.attributes["name"] = name;
		return &action;
	}

	Action* box()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_BOX;
		return &action;
	}

	Action* prism(const string& type)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_PRISM;
		action.attributes["type"] = type;
		return &action;
	}

	Action* ngon(const string& sides)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_NGON;
		action.attributes["sides"] = sides;
		return &action;
	}

	Action* center()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CENTER_PIVOT;
		action.attributes["dir"] = "xyz";
		return &action;
	}

	Action* center_x()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CENTER_PIVOT;
		action.attributes["dir"] = "x";
		return &action;
	}

	Action* center_y()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CENTER_PIVOT;
		action.attributes["dir"] = "y";
		return &action;
	}

	Action* center_z()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CENTER_PIVOT;
		action.attributes["dir"] = "z";
		return &action;
	}

	Action* center_xy()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CENTER_PIVOT;
		action.attributes["dir"] = "xy";
		return &action;
	}

	Action* center_yz()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CENTER_PIVOT;
		action.attributes["dir"] = "yz";
		return &action;
	}

	Action* center_xz()
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_CENTER_PIVOT;
		action.attributes["dir"] = "xz";
		return &action;
	}

	Action* move(const string& x, const string& y, const string& z)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_MOVE;
		action.attributes["x"] = x;
		action.attributes["y"] = y;
		action.attributes["z"] = z;
		return &action;
	}

	Action* rotate_x(const string& angle)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_ROTATE;
		action.attributes["x"] = angle;
		action.attributes["y"] = "0";
		action.attributes["z"] = "0";
		return &action;
	}

	Action* rotate_y(const string& angle)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_ROTATE;
		action.attributes["x"] = "0";
		action.attributes["y"] = angle;
		action.attributes["z"] = "0";
		return &action;
	}

	Action* rotate_z(const string& angle)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_ROTATE;
		action.attributes["x"] = "0";
		action.attributes["y"] = "0";
		action.attributes["z"] = angle;
		return &action;
	}

	Action* scale(const string& x, const string& y, const string& z)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_SCALE;
		action.attributes["x"] = x;
		action.attributes["y"] = y;
		action.attributes["z"] = z;
		return &action;
	}

	Action* occludes(const string& volumeid)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_OCCLUDE;
		action.attributes["volumeid"] = volumeid;
		return &action;
	}

	Action* material(const string& id)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_MATERIAL;
		action.attributes["id"] = id;
		return &action;
	}

	Action* local(const string& volumeid)
	{
		actions.push_back(Action());
		Action& action = actions.back();

		action.type = ActionType::ACTION_LOCAL;
		action.attributes["volumeid"] = volumeid;
		return &action;
	}

	void condition(const string& expression)
	{
		conditionals.push_back(expression);
	}

	void if_occluded_rule(const string& scope, const string& volumeId, const string& margin)
	{
		occlusionScopes.push_back(scope);
		occlusionVolumeIds.push_back(volumeId);
		occlusionMargins.push_back(margin);
	}
};

class Module
{
public:
	map<string, list<Rule>> rules;

	Module() {}

	Rule& AddNewRule(const string& name)
	{
		map<string, list<Rule>>::iterator it = rules.find(name);

		list<Rule>& subRules = it != rules.end() ? it->second : rules[name];
		subRules.push_back(Rule());
		return subRules.back();
	}

	void getCode(string& code, string& debuginfo)
	{
		code = "";
		debuginfo = "";
		
		string endl = "\n";

		char itoaBuf[64];

		for (auto const& kv : rules)
		{
			const string& ruleId = kv.first;
			const list<Rule>& subRules = kv.second;

			for (auto const& rule : subRules)
			{
				// Save rule ID
				code += ruleId + endl;

				// Save rule conditionals
				code += _itoa((int)rule.conditionals.size(), itoaBuf, 10) + endl;
				for (auto const& icond : rule.conditionals)
				{
					code += icond + endl;
				}

				// Save occlusion
				code += _itoa((int)rule.occlusionScopes.size(), itoaBuf, 10) + endl;
				for (auto const& s : rule.occlusionScopes)
				{
					code += s + endl;
				}
				code += _itoa((int)rule.occlusionVolumeIds.size(), itoaBuf, 10) + endl;
				for (auto const& s : rule.occlusionVolumeIds)
				{
					code += s + endl;
				}
				code += _itoa((int)rule.occlusionMargins.size(), itoaBuf, 10) + endl;
				for (auto const& s : rule.occlusionMargins)
				{
					code += s + endl;
				}

				/*
				// Save scope modes for occlusion tests
				f << rule.occlusionScopes.size() << endl;
				for (vector<string>::iterator icond = rule.occlusionScopes.begin();
				icond != rule.occlusionScopes.end(); ++icond)
				{
				f << *icond << endl;
				}

				// Save IDs for occlusion tests
				f << rule.occlusionVolumeIds.size() << endl;
				for (vector<string>::iterator icond = rule.occlusionVolumeIds.begin();
				icond != rule.occlusionVolumeIds.end(); ++icond)
				{
				f << *icond << endl;
				}

				// Save margins for occlusion tests
				f << rule.occlusionMargins.size() << endl;
				for (vector<string>::iterator icond = rule.occlusionMargins.begin();
				icond != rule.occlusionMargins.end(); ++icond)
				{
				f << *icond << endl;
				}
				*/

				// Save actions
				code += _itoa((int)rule.actions.size(), itoaBuf, 10) + endl;
				for (auto const& action : rule.actions)
				{
					// Get reference to action
					code += _itoa((int)(action.type), itoaBuf, 10) + endl;

					// Save action attributes
					code += _itoa((int)action.attributes.size(), itoaBuf, 10) + endl;
					for (auto const& kvAttr : action.attributes)
					{
						const string& attrId = kvAttr.first;
						code += attrId + endl;
						code += kvAttr.second + endl;
					}

					// Save debug info
					debuginfo += action.filename + endl;
					debuginfo += _itoa(action.line, itoaBuf, 10) + endl;
					debuginfo += _itoa(action.column, itoaBuf, 10) + endl;
				}
			}
		}
	}
};

static void ReplaceAll(string& str, const string& from, const string& to)
{
	if (from.empty())
		return;

	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
	}
}

// forward declaration
static void buildExpression(Reduction& reduction, /*not const*/ string& expression);
static void buildRule(
	GoldCPP::Reduction& reduction,
	/*not const*/ Rule& rule,
	Module& module,
	const string& filename);

static void buildParameterList(Reduction& reduction, string& list)
{
	Reduction& op = reduction;
	if (op.Parent->TableIndex == ProductionConstants::PROD_PARAMETERLIST_COMMA)
	{
		string parameter = "";
		buildExpression(*(op.get_Data(0).ReductionData), parameter);
		list += parameter + ",";
		buildParameterList(*(op.get_Data(2).ReductionData), list);
	}
	if (op.Parent->TableIndex == ProductionConstants::PROD_PARAMETERLIST)
	{
		string parameter = "";
		buildExpression(*(op.get_Data(0).ReductionData), parameter);
		list += parameter;
	}
}

static void buildExpression(Reduction& reduction, /*not const*/ string& expression)
{
	Reduction& op = reduction;
	if (op.Count() == 3 && op.get_Data(1).ReductionData == NULL)
	{
		string lop = "", rop = "", oper;
		buildExpression(*(op.get_Data(0).ReductionData), lop);
		buildExpression(*(op.get_Data(2).ReductionData), rop);
		oper = op.get_Data(1).StringData;
		expression += lop + "," + rop + "," + oper;
	}
	else if (op.Count() == 2 && op.get_Data(0).ReductionData == NULL)
	{
		string rop = "", oper;
		buildExpression(*(op.get_Data(1).ReductionData), rop);
		oper = op.get_Data(0).StringData;
		expression += rop + ",~" + oper;
	}
	else if (op.Parent->TableIndex == ProductionConstants::PROD_PERCENTAGE_PERCENT)
	{
		string rop = "";
		buildExpression(*(op.get_Data(0).ReductionData), rop);
		expression += rop + ",~%";
	}
	else if (op.Parent->TableIndex == ProductionConstants::PROD_FUNCTIONCALL_IDENTIFIER_LPAREN_RPAREN)
	{
		string functionName = op.get_Data(0).StringData;
		string parameters = "";
		buildParameterList(*(op.get_Data(2).ReductionData), parameters);
		if (parameters.length() > 0)
			expression += "?," + parameters + ",#" + functionName;
		else
			expression += "?,#" + functionName;
	}
	else if (op.Parent->TableIndex == ProductionConstants::PROD_VALUE_LPAREN_RPAREN)
	{
		if (op.get_Data(1).ReductionData != NULL)
			buildExpression(*(op.get_Data(1).ReductionData), expression);
		else
			expression += op.get_Data(1).StringData;
	}
	else
	{
		if (op.get_Data(0).ReductionData != NULL)
			buildExpression(*(op.get_Data(0).ReductionData), expression);
		else
			expression += op.get_Data(0).StringData;
	}
}

static void buildVector(Reduction& reduction, vector<string>& vect, short type)
{
	if (reduction.Parent->TableIndex == type)
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL &&
				reduction.get_Data(i).ReductionData->Parent->TableIndex == ProductionConstants::PROD_EXPRESSION)
			{
				string expression = "";
				buildExpression(*(reduction.get_Data(i).ReductionData), expression);
				vect.push_back(expression);
			}
		}
	}
	else
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL)
			{
				buildVector(*(reduction.get_Data(i).ReductionData), vect, type);
			}
		}
	}
}

static void buildVector2DList(Reduction& reduction, string& positions, string& widths)
{
	if (reduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR2D_LBRACKET_COMMA_RBRACKET)
	{
		Reduction& vectReduction = reduction;
		vector<string> vect;
		buildVector(vectReduction, vect, ProductionConstants::PROD_VECTOR2D_LBRACKET_COMMA_RBRACKET);
		positions += vect[0] + " ";
		widths += vect[1] + " ";
	}
	else if (reduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR2DLIST ||
		reduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR2DLIST2)
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL)
			{
				buildVector2DList(*(reduction.get_Data(i).ReductionData), positions, widths);
			}
		}
	}
}

static void buildVector3DList(Reduction& reduction, string& positions, string& widths, string& depths)
{
	if (reduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR3D_LBRACKET_COMMA_COMMA_RBRACKET)
	{
		Reduction& vectReduction = reduction;
		vector<string> vect;
		buildVector(vectReduction, vect, ProductionConstants::PROD_VECTOR3D_LBRACKET_COMMA_COMMA_RBRACKET);
		positions += vect[0] + " ";
		widths += vect[1] + " ";
		depths += vect[2] + " ";
	}
	else if (reduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR3DLIST ||
		reduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR3DLIST2)
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL)
			{
				buildVector3DList(*(reduction.get_Data(i).ReductionData), positions, widths, depths);
			}
		}
	}
}

static void getStringValue(Reduction& reduction, string& value)
{
	if (reduction.get_Data(0).ReductionData == NULL)
	{
		value = reduction.get_Data(0).StringData;
	}
	else
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL)
			{
				getStringValue(*(reduction.get_Data(i).ReductionData), value);
			}
		}
	}
}

static string NewGuidString()
{
	static int sFakeGuid = 10;
	char buf[64];
	return _itoa(sFakeGuid++, buf, 10);
}

static void buildModuleRef(Reduction& reduction, string& moduleList, Module& module, const string& filename)
{
	Reduction& moduleReduction = reduction;
	if (moduleReduction.Parent->TableIndex == ProductionConstants::PROD_MODULEREF_IDENTIFIER)
	{
		const string& name = moduleReduction.Branches[0]->StringData;
		moduleList = name;
	}
	else if (moduleReduction.Parent->TableIndex == ProductionConstants::PROD_MODULEREF)
	{
		string rulename = NewGuidString();
		rulename = "$" + rulename;
		moduleList = rulename;
		Rule& anon = module.AddNewRule(rulename);
		buildRule(moduleReduction, anon, module, filename);
	}
}

static void buildModuleSequence(Reduction& reduction, string& dim, string& snap, string& moduleList, Module& module, const string& filename)
{
	if (reduction.Parent->TableIndex == ProductionConstants::PROD_MODULEREF1D)
	{
		Reduction& vectReduction = *(reduction.get_Data(0).ReductionData);
		if (vectReduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR1D_LBRACKET_RBRACKET)
		{
			vector<string> vect;
			buildVector(vectReduction, vect, ProductionConstants::PROD_VECTOR1D_LBRACKET_RBRACKET);
			dim += vect[0] + " ";
			snap += "* ";
		}
		else if (vectReduction.Parent->TableIndex == ProductionConstants::PROD_VECTOR1D_LBRACKET_PIPE_IDENTIFIER_RBRACKET)
		{
			vector<string> vect;
			buildVector(vectReduction, vect, ProductionConstants::PROD_VECTOR1D_LBRACKET_PIPE_IDENTIFIER_RBRACKET);
			dim += vect[0] + " ";
			string snapLine = vectReduction.get_Data(3).StringData;
			snap += snapLine + " ";
		}
		else
			dim += "* ";

		Reduction& moduleReduction = *(reduction.get_Data(1).ReductionData);
		if (moduleReduction.Parent->TableIndex == ProductionConstants::PROD_MODULEREF_IDENTIFIER)
		{
			string name = moduleReduction.get_Data(0).StringData;
			moduleList += name + " ";
		}
		else if (moduleReduction.Parent->TableIndex == ProductionConstants::PROD_MODULEREF)
		{
			string rulename = NewGuidString();
			rulename = "$" + rulename;
			moduleList += rulename + " ";
			Rule& anon = module.AddNewRule(rulename);
			buildRule(moduleReduction, anon, module, filename);
		}
	}
	else if (reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESEQUENCE1D ||
		reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESEQUENCE1D2)
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL)
			{
				buildModuleSequence(*(reduction.get_Data(i).ReductionData), dim, snap, moduleList, module, filename);
			}
		}
	}
}

static void buildRule(
	GoldCPP::Reduction& reduction,
	/*not const*/ Rule& rule,
	Module& module,
	const string& filename)
{
	int type = reduction.Parent->TableIndex;
	Action* action = NULL;
	switch (type)
	{
	case ProductionConstants::PROD_RULESNAP_SNAP_IDENTIFIER:
	{
		string direction, identifier, expression = "";
		identifier = reduction.Branches[1]->StringData;
		ReplaceAll(identifier, "\"", "");

		int dir = (reduction.Branches[2]->ReductionData)->Parent->TableIndex;
		if (dir == ProductionConstants::PROD_DIRECTION_X)
			direction = "x";
		else if (dir == ProductionConstants::PROD_DIRECTION_Y)
			direction = "y";
		else // if (dir == ProductionConstants::PROD_DIRECTION_Z)
			direction = "z";
		buildExpression(*(reduction.Branches[3]->ReductionData), expression);
		action = rule.snap(direction, identifier, expression);
	}
	break;
	case ProductionConstants::PROD_ASSIGMENT_DEFINE_IDENTIFIER:
	{
		string identifier, expression = "";
		identifier = reduction.Branches[1]->StringData;
		ReplaceAll(identifier, "\"", "");
		buildExpression(*(reduction.get_Data(2).ReductionData), expression);
		action = rule.var(identifier, expression, false);
	}
	break;
	case ProductionConstants::PROD_ASSIGMENTCONST_CONST_IDENTIFIER:
	{
		string identifier, expression = "";
		identifier = reduction.Branches[1]->StringData;
		ReplaceAll(identifier, "\"", "");
		buildExpression(*(reduction.get_Data(2).ReductionData), expression);
		action = rule.var(identifier, expression, true);
	}
	break;
	case ProductionConstants::PROD_RULELOFT_LOFT:
	{
		string positions = "", widths = "", depths = "";
		buildVector3DList(*(reduction.get_Data(1).ReductionData), positions, widths, depths);
		action = rule.loft_box(positions, widths, depths);
	}
	break;
	case ProductionConstants::PROD_RULELOFTNGON_LOFT_NGON:
	{
		string sides = "";
		buildExpression(*(reduction.get_Data(2).ReductionData), sides);
		string positions = "", widths = "";
		buildVector2DList(*(reduction.get_Data(3).ReductionData), positions, widths);
		action = rule.loft_ngon(sides, positions, widths);
	}
	break;
	case ProductionConstants::PROD_RULESELECT_SELECT:
	{
		string subscope = "";
		getStringValue(*(reduction.get_Data(1).ReductionData), subscope);
		if (subscope.compare("xy") == 0)
			subscope = "face_x";
		else if (subscope.compare("yz") == 0)
			subscope = "face_z";
		else if (subscope.compare("prismsides") == 0)
			subscope = "prism_sides";
		string moduleName = "";
		buildModuleRef(*(reduction.get_Data(2).ReductionData), moduleName, module, filename);
		action = rule.select(subscope, moduleName);
	}
	break;
	case ProductionConstants::PROD_RULEINSTANCE_INSTANCE_STRINGLITERAL:
	{
		string instance;
		instance = reduction.Branches[1]->StringData;
		ReplaceAll(instance, "\"", "");
		action = rule.instance(instance);
	}
	break;
	case ProductionConstants::PROD_RULEREPEAT_REPEAT:
	{
		string dim = "", snap = "", modules = "";
		buildModuleSequence(*(reduction.get_Data(2).ReductionData), dim, snap, modules, module, filename);
		int dir = reduction.Branches[1]->ReductionData->Parent->TableIndex;
		if (dir == ProductionConstants::PROD_DIRECTION_X)
			action = rule.repeat_x(dim, modules, "flex", snap);
		else if (dir == ProductionConstants::PROD_DIRECTION_Y)
			action = rule.repeat_y(dim, modules, "flex", snap);
		else if (dir == ProductionConstants::PROD_DIRECTION_Z)
			action = rule.repeat_z(dim, modules, "flex", snap);
	}
	break;
	case ProductionConstants::PROD_RULEDIVIDE_DIVIDE:
	{
		string dim = "", snap = "", modules = "";
		buildModuleSequence(*(reduction.get_Data(2).ReductionData), dim, snap, modules, module, filename);
		int dir = reduction.Branches[1]->ReductionData->Parent->TableIndex;
		if (dir == ProductionConstants::PROD_DIRECTION_X)
			action = rule.subdivide_x(dim, modules);
		else if (dir == ProductionConstants::PROD_DIRECTION_Y)
			action = rule.subdivide_y(dim, modules);
		else if (dir == ProductionConstants::PROD_DIRECTION_Z)
			action = rule.subdivide_z(dim, modules);
	}
	break;
	case ProductionConstants::PROD_RULEDIVIDE_DIVIDE_PRISM:
	{
		string position = "", center = "", top = "", left = "", right = "";

		buildExpression(*(reduction.get_Data(2).ReductionData), position);
		buildModuleRef(*(reduction.get_Data(3).ReductionData), center, module, filename);
		buildModuleRef(*(reduction.get_Data(4).ReductionData), top, module, filename);
		buildModuleRef(*(reduction.get_Data(5).ReductionData), left, module, filename);
		buildModuleRef(*(reduction.get_Data(6).ReductionData), right, module, filename);

		action = rule.subdivide_prism(position, center, top, left, right);
	}
	break;
	case ProductionConstants::PROD_RULECALL_MODULE_IDENTIFIER:
	{
		string instance;
		instance = reduction.Branches[1]->StringData;
		action = rule.module(instance);
	}
	break;
	case ProductionConstants::PROD_RULECALL_MODULE:
	{
		rule.push();
		buildRule(*(reduction.Branches[1]->ReductionData), rule, module, filename);
		rule.pop();
	}
	break;
	case ProductionConstants::PROD_RULEDEFER_DEFER_IDENTIFIER:
	{
		string instance;
		instance = reduction.Branches[1]->StringData;
		action = rule.defer(instance);
	}
	break;
	case ProductionConstants::PROD_RULEBOX_BOX:
		action = rule.box();
		break;
	case ProductionConstants::PROD_RULEPRISM_PRISM:
		action = rule.prism("");
		break;
	case ProductionConstants::PROD_RULENGON_NGON:
	{
		string expression = "";
		buildExpression(*(reduction.get_Data(1).ReductionData), expression);
		action = rule.ngon(expression);
	}
	break;
	case ProductionConstants::PROD_RULECENTER_CENTER:
		action = rule.center();
		break;
	case ProductionConstants::PROD_RULECENTER_CENTER2:
	case ProductionConstants::PROD_RULECENTER_CENTER3:
	{
		// (((GOLD::Production^)((((GOLD::Reduction^)(((System::Object^)(([1])->Data)))))->Parent)))->m_TableIndex
		int dir = reduction.Branches[1]->ReductionData->Parent->TableIndex;
		if (dir == ProductionConstants::PROD_DIRECTION_X)
			action = rule.center_x();
		else if (dir == ProductionConstants::PROD_DIRECTION_Y)
			action = rule.center_y();
		else if (dir == ProductionConstants::PROD_DIRECTION_Z)
			action = rule.center_z();
		else if (dir == ProductionConstants::PROD_DIRECTION2D_XY)
			action = rule.center_xy();
		else if (dir == ProductionConstants::PROD_DIRECTION2D_YZ)
			action = rule.center_yz();
		else if (dir == ProductionConstants::PROD_DIRECTION2D_XZ)
			action = rule.center_xz();
	}
	break;
	case ProductionConstants::PROD_RULEMOVE_MOVE:
	{
		vector<string> vect;
		buildVector(reduction, vect, ProductionConstants::PROD_VECTOR3D_LBRACKET_COMMA_COMMA_RBRACKET);
		action = rule.move(vect[0], vect[1], vect[2]);
	}
	break;
	case ProductionConstants::PROD_RULEROTATE_ROTATE:
	{
		string expression = "";
		buildExpression(*(reduction.get_Data(2).ReductionData), expression);
		int dir = reduction.Branches[1]->ReductionData->Parent->TableIndex;
		if (dir == ProductionConstants::PROD_DIRECTION_X)
			action = rule.rotate_x(expression);
		else if (dir == ProductionConstants::PROD_DIRECTION_Y)
			action = rule.rotate_y(expression);
		else if (dir == ProductionConstants::PROD_DIRECTION_Z)
			action = rule.rotate_z(expression);
	}
	break;
	case ProductionConstants::PROD_RULESCALE_SCALE3:
	{
		vector<string> vect;
		buildVector(reduction, vect, ProductionConstants::PROD_VECTOR3D_LBRACKET_COMMA_COMMA_RBRACKET);
		action = rule.scale(vect[0], vect[1], vect[2]);
	}
	break;
	case ProductionConstants::PROD_RULEOCCLUDE_OCCLUDES_IDENTIFIER:
	{
		action = rule.occludes(reduction.Branches[1]->StringData);
	}
	break;
	case ProductionConstants::PROD_RULEMATERIAL_MATERIAL_IDENTIFIER:
	{
		action = rule.material(reduction.Branches[1]->StringData);
	}
	break;
	case ProductionConstants::PROD_RULELOCAL_LOCAL_IDENTIFIER:
	{
		action = rule.local(reduction.Branches[1]->StringData);
	}
	break;
	default:
		for (int i = 0; i < reduction.Branches.Count(); i++)
		{
			if (reduction.Branches[i]->ReductionData != NULL)
			{
				buildRule(*(reduction.Branches[i]->ReductionData), rule, module, filename);
			}
		}
		break;
	}

	if (action != NULL)
	{
		if (reduction.Branches.Count() > 0)
		{
			// add debug information
			action->filename = filename;
			GoldCPP::Token& token = *(reduction.Branches[0]);
			action->line = token.Pos.Line;
			action->column = token.Pos.Column;
		}
	}
}

static void buildSelectors(Reduction& reduction, Rule& rule, Module& module)
{
	if (reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESELECTOR ||
		reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESELECTOR_COMMA ||
		reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESELECTOR2 ||
		reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESELECTOR_COMMA2)
	{
		Reduction& condition = *(reduction.get_Data(0).ReductionData);
		if (condition.Parent->TableIndex == ProductionConstants::PROD_OCCLUSIONEXPRESSION_OCCLUDED_IDENTIFIER)
		{
			string scope, margin = "", volumeId;
			scope = (condition.get_Data(1).ReductionData)->get_Data(0).StringData;
			volumeId = condition.get_Data(2).StringData;
			buildExpression(*(condition.get_Data(3).ReductionData), margin);
			rule.if_occluded_rule(scope, volumeId, margin);
		}
		else
		{
			string expression = "";
			buildExpression(condition, expression);
			rule.condition(expression);
		}
		if (reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESELECTOR_COMMA ||
			reduction.Parent->TableIndex == ProductionConstants::PROD_MODULESELECTOR_COMMA2)
			buildSelectors(*(reduction.get_Data(2).ReductionData), rule, module);
	}
	else
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL)
			{
				buildSelectors(*(reduction.get_Data(i).ReductionData), rule, module);
			}
		}
	}
}

static void buildModule(Reduction& reduction, Module& module, const string& filename)
{
	if (reduction.Parent->TableIndex == ProductionConstants::PROD_MODULE_IDENTIFIER)
	{
		Rule* rule = NULL;
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.Branches[i] != NULL)
			{
				GoldCPP::Token& token = *(reduction.Branches[i]);
				if (token.Parent->TableIndex == SymbolIndex::SYM_IDENTIFIER)
				{
					string name = token.StringData;

					rule = &module.AddNewRule(name);
				}
				if (rule != NULL && token.Parent->TableIndex == SymbolIndex::SYM_INLINEMODULE)
				{
					buildRule(reduction, *rule, module, filename);
				}
				if (token.Parent->TableIndex == SymbolIndex::SYM_SELECTOR)
				{
					buildSelectors(*(token.ReductionData), *rule, module);
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < reduction.Count(); i++)
		{
			if (reduction.get_Data(i).ReductionData != NULL)
			{
				buildModule(*(reduction.get_Data(i).ReductionData), module, filename);
			}
		}
	}
}


bool ModCompiler::Compile(std::vector<SourceFile>& files, std::string& outCode, std::string& outDbg, std::vector<GrammarError>& errors)
{
	LoadEgtIfNeeded();

	if (!isEgtLoaded)
		return false;

	errors.clear();
	Module main;
	bool error = false;
	bool hasCode = false;

	for (size_t i = 0; i < files.size(); i++)
	{
		SourceFile& file = files[i];

		parser.Open(file.modSrc);

		GoldCPP::ParseMessage status = GoldCPP::ParseMessage::NotLoadedError;
		bool done = false;
		while (!done)
		{
			status = parser.Parse();
			done =
				status == GoldCPP::ParseMessage::Accept ||
				status == GoldCPP::ParseMessage::GroupError ||
				status == GoldCPP::ParseMessage::InternalError ||
				status == GoldCPP::ParseMessage::LexicalError ||
				status == GoldCPP::ParseMessage::NotLoadedError ||
				status == GoldCPP::ParseMessage::SyntaxError;
		}
		if (status == GoldCPP::ParseMessage::LexicalError || status == GoldCPP::ParseMessage::SyntaxError)
		{
			// add error to report
			GrammarError e;
			e.filename = file.filename;
			e.line = parser.GetCurrentPosition().Line;
			e.column = parser.GetCurrentPosition().Column;
			e.message = "Expected " + parser.GetExpectedSymbols().GetText();
			errors.push_back(e);
			error = true;
		}
		else if (status == GoldCPP::ParseMessage::Accept)
		{
			buildModule(*parser.GetCurrentReduction(), main, file.filename);
			hasCode = true;
		}
	}

	if (hasCode && !error)
	{
		main.getCode(outCode, outDbg);
	}
	else
	{
		outCode = "";
		outDbg = "";
	}
	return !error;
}


void ModCompiler::LoadEgtIfNeeded()
{
	if (isEgtLoaded)
		return;

	ifstream egtInput(egtPath, ios::binary);
	if (egtInput.fail())
		return;

	vector<char> egtBuffer((istreambuf_iterator<char>(egtInput)), (istreambuf_iterator<char>()));
	egtInput.close();

	if (!parser.LoadTables((uint8_t*)egtBuffer.data(), egtBuffer.size()))
		return;

	isEgtLoaded = true;
}


}

