/************************************************************
* (C) Voxel Farm Inc. 2015
*/
#include "HyVoxelPrivatePCH.h"

#include "HyVoxelConfig.h"
#include "ProcModeling/Architecture/Grammar.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <math.h>

#include "Common/Matrix.h"

#include <string.h>
#include <stdio.h>

#include "Common/vfminmax.h"

using namespace HyVoxel;
using namespace HyVoxel::Algebra;
using namespace HyVoxel::Architecture;

#define is_operator(c)  (c == '+' || c == '-' || c == '/' || c == '*' || c == '!' || c == '%' || c == '=')
#define is_function(c)  (c == '#')
#define is_argument_chain(c)  (c == '?')
#define is_ident(c)     ((c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c == '.') || (c == '~'))

// CNode
CNode::~CNode()
{
	for (TVFVector<CNode*>::iterator i = children.begin(); i != children.end(); ++i)
	{
		CNode* node = *i;
		VF_DELETE node;
	}
}

void CNode::newGeneration(CNode* child, VFString module)
{
	recursePending.push_back(std::pair<CNode*, VFString>(child, module));
	child->generation = generation + 1;
}

void CNode::defer(CNode* child, VFString module)
{
	deferred.push_back(std::pair<CNode*, VFString>(child, module));
}


// CRule

CRule& CRule::var(VFString identifier, VFString value, bool constant)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_VAR;
	action.attributes["id"] = identifier;
	action.attributes["value"] = value;
	if (constant)
	{
		action.attributes["const"] = "true";
	}
	else
	{
		action.attributes["const"] = "false";
	}
	return *this;
}

CRule& CRule::snap(VFString dir, VFString id, VFString radius)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_SNAP;
	action.attributes["dir"] = dir;
	action.attributes["id"] = id;
	action.attributes["radius"] = radius;
	return *this;
}

CRule& CRule::condition(VFString exp)
{
	conditionals.push_back(exp);
	return *this;
}

CRule& CRule::scale(VFString x, VFString y, VFString z)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_SCALE;
	action.attributes["x"] = x;
	action.attributes["y"] = y;
	action.attributes["z"] = z;
	return *this;
}

CRule& CRule::move(VFString x, VFString y, VFString z)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_MOVE;
	action.attributes["x"] = x;
	action.attributes["y"] = y;
	action.attributes["z"] = z;
	return *this;
}

CRule& CRule::rotate_x(VFString angle)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_ROTATE;
	action.attributes["x"] = angle;
	action.attributes["y"] = "0";
	action.attributes["z"] = "0";
	return *this;
}

CRule& CRule::rotate_y(VFString angle)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_ROTATE;
	action.attributes["x"] = "0";
	action.attributes["y"] = angle;
	action.attributes["z"] = "0";
	return *this;
}

CRule& CRule::rotate_z(VFString angle)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_ROTATE;
	action.attributes["x"] = "0";
	action.attributes["y"] = "0";
	action.attributes["z"] = angle;
	return *this;
}

CRule& CRule::align(VFString dir)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_ALIGN;
	action.attributes["dir"] = dir;
	return *this;
}

CRule& CRule::center()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CENTER_PIVOT;
	action.attributes["dir"] = "xyz";
	return *this;
}

CRule& CRule::center_x()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CENTER_PIVOT;
	action.attributes["dir"] = "x";
	return *this;
}

CRule& CRule::center_y()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CENTER_PIVOT;
	action.attributes["dir"] = "y";
	return *this;
}

CRule& CRule::center_z()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CENTER_PIVOT;
	action.attributes["dir"] = "z";
	return *this;
}

CRule& CRule::center_xy()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CENTER_PIVOT;
	action.attributes["dir"] = "xy";
	return *this;
}

CRule& CRule::center_xz()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CENTER_PIVOT;
	action.attributes["dir"] = "xz";
	return *this;
}

CRule& CRule::center_yz()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CENTER_PIVOT;
	action.attributes["dir"] = "yz";
	return *this;
}

CRule& CRule::module(VFString name, VFString generations)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_CALL;
	action.attributes["name"] = name;
	action.attributes["generations"] = generations;
	return *this;
}

CRule& CRule::defer(VFString name)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_DEFER;
	action.attributes["name"] = name;
	return *this;
}

CRule& CRule::instance(VFString geometry)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_INSTANCE;
	action.attributes["geometry"] = geometry;
	return *this;
}

CRule& CRule::material(VFString id)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_MATERIAL;
	action.attributes["id"] = id;
	return *this;
}

CRule& CRule::box()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_BOX;
	return *this;
}

CRule& CRule::ngon(VFString sides)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_NGON;
	action.attributes["sides"] = sides;
	return *this;
}

CRule& CRule::prism(VFString type)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_PRISM;
	action.attributes["type"] = type;
	return *this;
}

CRule& CRule::loft_box(VFString positions, VFString widths, VFString depths)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_LOFT_BOX;
	action.attributes["positions"] = positions;
	action.attributes["widths"] = widths;
	action.attributes["depths"] = depths;
	return *this;
}

CRule& CRule::loft_ngon(VFString sides, VFString positions, VFString radii)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_LOFT_NGON;
	action.attributes["sides"] = sides;
	action.attributes["positions"] = positions;
	action.attributes["radii"] = radii;
	return *this;
}

CRule& CRule::select(VFString stype, VFString module)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_SELECT;
	action.attributes["dimension"] = "face";
	action.attributes["type"] = stype;
	action.attributes["module"] = module;
	return *this;
}

CRule& CRule::if_occluded_rule(VFString scope, VFString volumeId, VFString margin)
{
	occlusionVolumeIds.push_back(volumeId);
	occlusionMargins.push_back(margin);
	occlusionScopes.push_back(scope);
	return *this;
}

CRule& CRule::occludes(VFString volumeid)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_OCCLUDE;
	action.attributes["volumeid"] = volumeid;
	return *this;
}

CRule& CRule::local(VFString volumeid)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_LOCAL;
	action.attributes["volumeid"] = volumeid;
	return *this;
}

CRule& CRule::subdivide_x(VFString pattern, VFString modules)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_DIVIDE;
	action.attributes["type"] = "x";
	action.attributes["pattern"] = pattern;
	action.attributes["modules"] = modules;
	return *this;
}

CRule& CRule::subdivide_y(VFString pattern, VFString modules)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_DIVIDE;
	action.attributes["type"] = "y";
	action.attributes["pattern"] = pattern;
	action.attributes["modules"] = modules;
	return *this;
}

CRule& CRule::subdivide_z(VFString pattern, VFString modules)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_DIVIDE;
	action.attributes["type"] = "z";
	action.attributes["pattern"] = pattern;
	action.attributes["modules"] = modules;
	return *this;
}

CRule& CRule::subdivide_prism(VFString position, VFString center, VFString top, VFString left, VFString right)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_DIVIDE;
	action.attributes["type"] = "prism";
	action.attributes["position"] = position;
	action.attributes["center"] = center;
	action.attributes["top"] = top;
	action.attributes["left"] = left;
	action.attributes["right"] = right;
	return *this;
}

CRule& CRule::repeat_x(VFString length, VFString module, VFString flex, VFString snap)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_REPEAT;
	action.attributes["dir"] = "x";
	action.attributes["length"] = length;
	action.attributes["snap"] = snap;
	action.attributes["module"] = module;
	action.attributes["flex"] = flex;
	return *this;
}

CRule& CRule::repeat_y(VFString length, VFString module, VFString flex, VFString snap)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_REPEAT;
	action.attributes["dir"] = "y";
	action.attributes["length"] = length;
	action.attributes["snap"] = snap;
	action.attributes["module"] = module;
	action.attributes["flex"] = flex;
	return *this;
}

CRule& CRule::repeat_z(VFString length, VFString module, VFString flex, VFString snap)
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_REPEAT;
	action.attributes["dir"] = "z";
	action.attributes["length"] = length;
	action.attributes["snap"] = snap;
	action.attributes["module"] = module;
	action.attributes["flex"] = flex;
	return *this;
}

CRule& CRule::push()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_PUSH;
	return *this;
}

CRule& CRule::pop()
{
	actions.push_back(CAction());
	CAction& action = actions.back();
	action.type = ACTION_POP;
	return *this;
}

CModule& CRule::end()
{
	return *parent;
}


// CModule

CModule& CModule::axiom(VFString name)
{
	baseAxiom = name;
	return *this;
}

CModule& CModule::add(CRule& rule)
{
	return *this;
}

CRule& CModule::begin(VFString name)
{
	TVFVector<CRule>& rvec = rules[name];
	rvec.push_back(CRule());
	CRule& rule = rvec.back();
	rule.parent = this;
	return rule;
}

void CModule::saveToFile(const char* filename)
{
	// Create file
	std::ofstream f(filename, std::ios::out);

	// Save each set of rule
	for (TVFMap<VFString, TVFVector<CRule> >::iterator irule = rules.begin();
			irule != rules.end(); ++irule)
	{
		// Iterate over rule set
		TVFVector<CRule>& subrules = irule->second;
		for (TVFVector<CRule>::iterator isubrule = subrules.begin(); isubrule != subrules.end(); ++isubrule)
		{
			// Get reference to rule
			CRule& rule = *isubrule;

			// Save rule ID
			f << irule->first << std::endl;

			// Save rule conditionals
			f << rule.conditionals.size() << std::endl;
			for (TVFVector<VFString>::iterator icond = rule.conditionals.begin();
					icond != rule.conditionals.end(); ++icond)
			{
				f << *icond << std::endl;
			}

			// Save scope modes for occlusion tests
			f << rule.occlusionScopes.size() << std::endl;
			for (TVFVector<VFString>::iterator icond = rule.occlusionScopes.begin();
					icond != rule.occlusionScopes.end(); ++icond)
			{
				f << *icond << std::endl;
			}

			// Save IDs for occlusion tests
			f << rule.occlusionVolumeIds.size() << std::endl;
			for (TVFVector<VFString>::iterator icond = rule.occlusionVolumeIds.begin();
					icond != rule.occlusionVolumeIds.end(); ++icond)
			{
				f << *icond << std::endl;
			}

			// Save margins for occlusion tests
			f << rule.occlusionMargins.size() << std::endl;
			for (TVFVector<VFString>::iterator icond = rule.occlusionMargins.begin();
					icond != rule.occlusionMargins.end(); ++icond)
			{
				f << *icond << std::endl;
			}

			// Save actions
			f << rule.actions.size() << std::endl;
			for (TVFVector<CAction>::iterator iaction = rule.actions.begin();
					iaction != rule.actions.end(); ++iaction)
			{
				// Get reference to action
				CAction& action = *iaction;
				f << action.type << std::endl;

				// Save action attributes
				f << action.attributes.size() << std::endl;
				for (TVFMap<VFString, VFString>::iterator iattr = action.attributes.begin();
						iattr != action.attributes.end(); ++iattr)
				{
					f << iattr->first << std::endl;
					f << iattr->second << std::endl;
				}
			}
		}
	}
	f << std::endl;
	f.close();
}

template <typename TStream>
void SafeGetLine(TStream& f, std::string& line, int lineSize)
{
	std::getline(f, line);
	if (!line.empty() && *line.rbegin() == '\r') {
		line.erase(line.length() - 1, 1);
	}
}

template <typename TStream>
void LoadModuleFromStream(CModule& module, TStream& f)
{
	const int LINE_SIZE = 512;
	int actionCount = 0;
	bool eof = false;
	while (!eof)
	{
		eof = f.eof();

		// Read rule ID
		std::string line;
		SafeGetLine(f, line, LINE_SIZE);
		eof |= (line[0] == 0);
		CRule& rule = module.begin(line);

		// Read conditionals
		SafeGetLine(f, line, LINE_SIZE);
		eof |= (line[0] == 0);
		int condCount = atoi(line.c_str());
		for (int i = 0; i < condCount && !f.eof(); i++)
		{
			SafeGetLine(f, line, LINE_SIZE);
			rule.conditionals.push_back(line);
		}

		// Read occlusion scopes
		SafeGetLine(f, line, LINE_SIZE);
		eof |= (line[0] == 0);
		condCount = atoi(line.c_str());
		for (int i = 0; i < condCount && !f.eof(); i++)
		{
			SafeGetLine(f, line, LINE_SIZE);
			rule.occlusionScopes.push_back(line);
		}

		// Read occlusion ids
		SafeGetLine(f, line, LINE_SIZE);
		eof |= (line[0] == 0);
		condCount = atoi(line.c_str());
		for (int i = 0; i < condCount && !f.eof(); i++)
		{
			SafeGetLine(f, line, LINE_SIZE);
			rule.occlusionVolumeIds.push_back(line);
		}

		// Read occlusion margin
		SafeGetLine(f, line, LINE_SIZE);
		eof |= (line[0] == 0);
		condCount = atoi(line.c_str());
		for (int i = 0; i < condCount && !f.eof(); i++)
		{
			SafeGetLine(f, line, LINE_SIZE);
			rule.occlusionMargins.push_back(line);
		}

		// Read actions
		SafeGetLine(f, line, LINE_SIZE);
		eof |= (line[0] == 0);
		actionCount = atoi(line.c_str());
		for (int i = 0; i < actionCount && !f.eof(); i++)
		{
			rule.actions.push_back(CAction());
			CAction& action = rule.actions.back();
			SafeGetLine(f, line, LINE_SIZE);
			eof |= (line[0] == 0);
			action.type = (ActionType)atoi(line.c_str());

			// Read attributes
			SafeGetLine(f, line, LINE_SIZE);
			eof |= (line[0] == 0);
			int attrCount = atoi(line.c_str());
			for (int j = 0; j < attrCount && !f.eof(); j++)
			{
				std::string name;
				SafeGetLine(f, name, LINE_SIZE);
				std::string value;
				SafeGetLine(f, value, LINE_SIZE);
				action.attributes[name] = value;
			}
		}
	}
}

void CModule::loadFromFile(const char* filename, char* debugfile)
{
	// Open file
	std::ifstream f(filename, std::ios::in);
	if (!f.is_open())
		return;

	LoadModuleFromStream<>(*this, f);
}

void CModule::loadFromString(std::string& strBuf)
{
	std::istringstream iss(strBuf);

	LoadModuleFromStream<>(*this, iss);
}

// CGrammar

CModule& CGrammar::define(VFString name)
{
	return modules[name];
}


// CEvaluator

bool intersectSegPlane(Vector& N, Vector& p, Vector& sp1, Vector& sp2, Vector& pi)
{
	float d = Vector_dot(N, Vector_subtract(sp2, sp1));

	if (abs(d) < 0.000001)
	{
		return false;
	}

	float t = Vector_dot(N, Vector_subtract(p, sp1))/d;

	pi.x = sp1.x + t*(sp2.x - sp1.x);
	pi.y = sp1.y + t*(sp2.y - sp1.y);
	pi.z = sp1.z + t*(sp2.z - sp1.z);

	return (t >= 0.0f && t <= 1.0f);
}

bool intersectSegs(Vector& a1, Vector& a2, Vector& b1, Vector& b2, Vector& pi)
{
	Vector va = Vector_subtract(a2, a1);
	Vector vb = Vector_subtract(b2, b1);
	Vector vab = Vector_subtract(b1, a1);
	Vector n1 = Vector_cross(va, vb);
	Vector n2 = Vector_cross(vab, vb);

	float m1 = Vector_magnitude(n1);
	if (m1 < 0.000001f)
	{
		return false;
	}
	n1.x /= m1;
	n1.y /= m1;
	n1.z /= m1;

	float m2 = Vector_magnitude(n2);
	n2.x /= m2;
	n2.y /= m2;
	n2.z /= m2;

	float d = Vector_dot(n1, n2);
	if (abs(d) < 0.9999f)
	{
		return false;
	}

	float t = m2/m1;
	if (d < 0.0f)
	{
		t = -t;
	}
	if (t < 0.0f || t > 1.0f)
	{
		return false;
	}

	pi.x = a1.x + t*va.x;
	pi.y = a1.y + t*va.y;
	pi.z = a1.z + t*va.z;

	d = Vector_dot(
			Vector_subtract(b1, pi),
			Vector_subtract(b2, pi));

	return (d < 0.0f);
}

// Separating axis test
bool sat(float c1[4][2], float c2[4][2])
{
	//Debug
	bool ret = true;

	//For every face in c1
	for (int i = 0; i < 4; i++)
	{
		//Grab a face (face x, face y)
		float fx = c1[i][0] - c1[(i + 1) % 4][0];
		float fy = c1[i][1] - c1[(i + 1) % 4][1];

		//Create a perpendicular axis to project on (axis x, axis y)
		float ax = -fy, ay = fx;

		//Normalize the axis
		float len_v = sqrt(ax * ax + ay * ay);
		ax /= len_v;
		ay /= len_v;

		//Carve out the min and max values
		float c1_min = (std::numeric_limits<float>::max)();
		float c1_max = -c1_min;
		float c2_min = (std::numeric_limits<float>::max)();
		float c2_max = -c2_min;

		//Project every point in c1 on the axis and store min and max
		for (int j = 0; j < 4; j++)
		{
			float c1_proj = (ax * (c1[j][0]) + ay * (c1[j][1])) / (ax * ax + ay * ay);
			c1_min = std::min(c1_proj, c1_min);
			c1_max = std::max(c1_proj, c1_max);
		}

		//Project every point in c2 on the axis and store min and max
		for (int j = 0; j < 4; j++)
		{
			float c2_proj = (ax * (c2[j][0]) + ay * (c2[j][1])) / (ax * ax + ay * ay);
			c2_min = std::min(c2_proj, c2_min);
			c2_max = std::max(c2_proj, c2_max);
		}

		//Return if the projections do not overlap
		if (!(c1_max >= c2_min && c1_min <= c2_max))
		{
			ret = false;    //return false;
		}
	}
	return ret; //return true;
}

CEvaluator::CEvaluator() :rootOffsetX(0), rootOffsetY(0), rootOffsetZ(0), debugger(NULL), debugAction(DEBUG_STEP) {}

void CEvaluator::findSnapPoints(VFString& id, Vector segments[4][2], TVFVector<float>& points)
{
	// Find list of planes for the specified ID
	TVFMap<VFString, TVFVector<SnapPlane> >::iterator planes = snapPlanes.find(id);

	// If no planes are found, exit
	if (planes == snapPlanes.end())
	{
		return;
	}

	// For each plane:
	for (TVFVector<SnapPlane>::iterator i = planes->second.begin(); i != planes->second.end(); ++i)
	{
		// Get reference to plane
		SnapPlane& snapPlane = *i;

		// Compute normal
		Vector N = Vector_cross(
					   Vector_subtract(snapPlane.s1p2, snapPlane.s1p1),
					   Vector_subtract(snapPlane.s2p1, snapPlane.s1p1));

		// See if plane is parallel to the provided segments
		bool parallel = false;
		Vector segInt[4];
		for (int si = 0; si < 4 && !parallel; si++)
		{
			parallel = !intersectSegPlane(N, snapPlane.s1p2, segments[si][0], segments[si][1], segInt[si]);
		}

		// If it is parallel, there will be no intersection
		if (parallel)
		{
			continue;
		}

		// Project plane to the best aligned main coordenate plane
		float rectPlane[4][2];
		float rectSegs[4][2];
		float nx = abs(N.x);
		float ny = abs(N.y);
		float nz = abs(N.z);
		if (nx > ny && nx > nz)
		{
			rectPlane[0][0] = snapPlane.s1p1.y;
			rectPlane[0][1] = snapPlane.s1p1.z;
			rectPlane[1][0] = snapPlane.s1p2.y;
			rectPlane[1][1] = snapPlane.s1p2.z;
			rectPlane[2][0] = snapPlane.s2p2.y;
			rectPlane[2][1] = snapPlane.s2p2.z;
			rectPlane[3][0] = snapPlane.s2p1.y;
			rectPlane[3][1] = snapPlane.s2p1.z;
			rectSegs[0][0] = segInt[0].y;
			rectSegs[0][1] = segInt[0].z;
			rectSegs[1][0] = segInt[1].y;
			rectSegs[1][1] = segInt[1].z;
			rectSegs[2][0] = segInt[3].y;
			rectSegs[2][1] = segInt[3].z;
			rectSegs[3][0] = segInt[2].y;
			rectSegs[3][1] = segInt[2].z;
		}
		else if (ny > nx && ny > nz)
		{
			rectPlane[0][0] = snapPlane.s1p1.x;
			rectPlane[0][1] = snapPlane.s1p1.z;
			rectPlane[1][0] = snapPlane.s1p2.x;
			rectPlane[1][1] = snapPlane.s1p2.z;
			rectPlane[2][0] = snapPlane.s2p2.x;
			rectPlane[2][1] = snapPlane.s2p2.z;
			rectPlane[3][0] = snapPlane.s2p1.x;
			rectPlane[3][1] = snapPlane.s2p1.z;
			rectSegs[0][0] = segInt[0].x;
			rectSegs[0][1] = segInt[0].z;
			rectSegs[1][0] = segInt[1].x;
			rectSegs[1][1] = segInt[1].z;
			rectSegs[2][0] = segInt[3].x;
			rectSegs[2][1] = segInt[3].z;
			rectSegs[3][0] = segInt[2].x;
			rectSegs[3][1] = segInt[2].z;
		}
		else
		{
			rectPlane[0][0] = snapPlane.s1p1.y;
			rectPlane[0][1] = snapPlane.s1p1.x;
			rectPlane[1][0] = snapPlane.s1p2.y;
			rectPlane[1][1] = snapPlane.s1p2.x;
			rectPlane[2][0] = snapPlane.s2p2.y;
			rectPlane[2][1] = snapPlane.s2p2.x;
			rectPlane[3][0] = snapPlane.s2p1.y;
			rectPlane[3][1] = snapPlane.s2p1.x;
			rectSegs[0][0] = segInt[0].y;
			rectSegs[0][1] = segInt[0].x;
			rectSegs[1][0] = segInt[1].y;
			rectSegs[1][1] = segInt[1].x;
			rectSegs[2][0] = segInt[3].y;
			rectSegs[2][1] = segInt[3].x;
			rectSegs[3][0] = segInt[2].y;
			rectSegs[3][1] = segInt[2].x;
		}

		// Perform separation axis test
		if (!sat(rectPlane, rectSegs) && !sat(rectSegs, rectPlane))
		{
			continue;
		}

		// Compute intersections
		for (int si = 0; si < 4; si++)
		{
			float distToIntersection = Vector_magnitude(Vector_subtract(segInt[si], segments[si][0]));

			TVFVector<float>::iterator pos = points.begin();
			if (pos != points.end())
			{
				float distToSnap;
				do
				{
					float spi = *pos;
					distToSnap = spi;
					if (distToSnap < distToIntersection)
					{
						++pos;
					}
				}
				while (pos != points.end() && distToSnap < distToIntersection);
			}

			points.insert(pos, distToIntersection);
		}
	}
}

bool CEvaluator::evaluatePending(CNode* node, int maxgeneration)
{
	// See if the generation is unde the specified limit
	if (maxgeneration != -1 && node->generation >= maxgeneration)
	{
		return false;
	}

	// Collect visited nodes
	bool pending = false;
	HyVoxel::TVFSet<CNode*> visited;

	// Iterate over all the pending nodes
	for (TVFVector<std::pair<CNode*, VFString> >::iterator i = node->recursePending.begin();
			i != node->recursePending.end(); ++i)
	{
		// Get child node and process it
		CNode* child = i->first;
		VFString module = i->second;
		currentNode = child;
		recurse(module);
		pending = true;

		// Remember it as visited
		visited.insert(child);
	}
	// Clear list of pending nodes
	node->recursePending.clear();

	// Iterate over children nodes
	for (TVFVector<CNode*>::iterator i = node->children.begin();
			i != node->children.end(); ++i)
	{
		// If the node was not visited before, evaluate pending subnodes for it
		if (visited.find(*i) == visited.end())
			if (evaluatePending(*i, maxgeneration))
			{
				pending = true;
			}
	}

	currentNode = node;
	return pending;
}

bool CEvaluator::evaluateDeferred(CNode* node, int maxgeneration)
{
	bool pending = false;
	HyVoxel::TVFSet<CNode*> visited;

	// Iterate over all the deferred nodes
	for (TVFVector<std::pair<CNode*, VFString> >::iterator i = node->deferred.begin();
			i != node->deferred.end(); ++i)
	{
		// Get child node and process it
		CNode* child = i->first;
		VFString module = i->second;
		currentNode = child;
		recurse(module);
		pending = true;

		// Remember it as visited
		visited.insert(child);
	}
	node->deferred.clear();

	// Iterate over children nodes
	for (TVFVector<CNode*>::iterator i = node->children.begin();
			i != node->children.end(); ++i)
	{
		// If the node was not visited before, evaluate pending subnodes for it
		if (visited.find(*i) == visited.end())
			if (evaluateDeferred(*i, maxgeneration))
			{
				pending = true;
			}
	}

	currentNode = node;
	return pending;
}

void CEvaluator::run(CGrammar& grammar, VFString module, CreateInstance* instancer, int maxgeneration)
{
	this->instancer = instancer;
	currentGrammar = &grammar;
	currentModule = &grammar.modules[module];
	currentNode = &root;

	// Evaluate root node
	recurse(currentModule->baseAxiom);

	// Iterate until all nodes are evaluated
	bool done = false;
	while (!done)
	{
		// Evaluated directly linked nodes first
		while (evaluatePending(currentNode, maxgeneration));
		// Evaluated any remaining deferred nodes
		done = !evaluateDeferred(currentNode, maxgeneration);
	}

	// Delete occluder indices
	for (TVFMap<VFString, OccluderIndex*>::iterator i = occluders.begin();
			i != occluders.end(); ++i)
	{
		VF_DELETE i->second;
	}
	occluders.clear();
}

template <class T>
bool from_string(T& t,
				 const VFString& s,
				 std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

bool isOccluded(CNode* node, VFString volumeId);

void HyVoxel::Architecture::updateScope(CNode* node, Matrix* matrix, bool useCache)
{
	// See if matrix cache can be used
	if (useCache && node->matrixCacheValid)
	{
		*matrix = node->cachedMatrix;
		return;
	}

	// Update parent's scope
	if (node->parent != NULL)
	{
		updateScope(node->parent, matrix);
	}

	// Translate and rotate scope
	Matrix_translate(matrix, node->scope.x + node->scope.px, node->scope.y + node->scope.py, node->scope.z + node->scope.pz);
	Matrix_rotate(matrix, Vector_withValues(1.f, 0.f, 0.f), (float)(-node->scope.ax*piover180f));
	Matrix_rotate(matrix, Vector_withValues(0.f, 0.f, 1.f), (float)(-node->scope.az*piover180f));
	Matrix_rotate(matrix, Vector_withValues(0.f, 1.f, 0.f), (float)(-node->scope.ay*piover180f));
	Matrix_translate(matrix, -node->scope.px, -node->scope.py, -node->scope.pz);
	Matrix_translate(matrix, node->scope.dx, node->scope.dy, node->scope.dz);

	// Apply aligment rules if any
	if (node->scope.xaligned || node->scope.yaligned || node->scope.zaligned)
	{
		Vector origin = Vector_withValues(0.f, 0.f, 0.f);
		origin = Matrix_multiplyVector(*matrix, origin);

		if (node->scope.yaligned)
		{
			Vector xaxis0 = Vector_withValues(0.f, 0.f, 0.f);
			Vector xaxis1 = Vector_withValues(1.f, 0.f, 0.f);
			xaxis0 = Matrix_multiplyVector(*matrix, xaxis0);
			xaxis1 = Matrix_multiplyVector(*matrix, xaxis1);
			Vector tx = Vector_subtract(xaxis1, xaxis0);
			Vector pxz = Vector_withValues(tx.x, 0.f, tx.z);
			Vector_normalize(&pxz);

			Vector zaxis0 = Vector_withValues(0.f, 0.f, 0.f);
			Vector zaxis1 = Vector_withValues(0.f, 0.f, 1.f);
			zaxis0 = Matrix_multiplyVector(*matrix, zaxis0);
			zaxis1 = Matrix_multiplyVector(*matrix, zaxis1);
			Vector tz = Vector_subtract(zaxis1, zaxis0);
			pxz = Vector_withValues(tz.x, 0.f, tz.z);
			Vector_normalize(&pxz);

			Matrix_loadIdentity(matrix);
			Matrix_translate(matrix, origin.x, origin.y, origin.z);

			Vector_normalize(&pxz);

			float d = Vector_dot(Vector_withValues(1.f, 0.f, 0.f), pxz);
			float d2 = Vector_dot(Vector_withValues(0.f, 0.f, 1.f), pxz);
			float ay = acos(d);
			if (d < 0.f && d2 > 0.f)
			{
				Matrix_rotate(matrix, Vector_withValues(0.f, 1.f, 0.f), -ay);
			}
			else if (d < 0.f && d2 < 0.f)
			{
				Matrix_rotate(matrix, Vector_withValues(0.f, 1.f, 0.f), ay);
			}
			else if (d > 0.f && d2 < 0.f)
			{
				Matrix_rotate(matrix, Vector_withValues(0.f, 1.f, 0.f), ay);
			}
			else // (d > 0.f && d2 > 0.f)
			{
				Matrix_rotate(matrix, Vector_withValues(0.f, 1.f, 0.f), -ay);
			}
		}
	}
	// Cache matrix
	if (useCache)
	{
		node->cachedMatrix = *matrix;
		node->matrixCacheValid = true;
	}
}

float CEvaluator::eval(VFString exp, float total, CNode* startNode)
{
	// If no node was provide, use current evaluation node
	if (startNode == NULL)
	{
		startNode = currentNode;
	}

	// Get c-string for expression
	const char* rpnExp = exp.c_str();

	// Get pointer for scanning the string
	const char* pos = rpnExp;

	// Get pointer to the end of the string
	const char* strend = rpnExp + strlen(rpnExp);

	// Test whether the expression is just a literal number
	while (pos < strend && *pos != ',')
	{
		pos++;
	}
	if (pos == strend)
	{
		float result = (float)atof(rpnExp);
		if (result != 0.0f || result == 0.0f && rpnExp[0] == '0')
		{
			return result;
		}
	}

	// Expression comes in polish notation, evaluate using stack
	pos = rpnExp;
	TVFVector<float> stack;
	TVFVector<int> argumentStack;
	while (pos < strend)
	{
		char token[256];
		int p = 0;
		while (pos < strend && *pos != ',')
		{
			if (*pos != ' ')
			{
				token[p++] = *pos;
			}
			pos++;
		}
		token[p] = 0;
		pos++;

		char c = token[0];
		if (c == '+')
		{
			float op2 = stack.back();
			stack.pop_back();
			float op1 = stack.back();
			stack.pop_back();
			stack.push_back(op1 + op2);
		}
		else if (c == '-')
		{
			float op2 = stack.back();
			stack.pop_back();
			float op1 = stack.back();
			stack.pop_back();
			stack.push_back(op1 - op2);
		}
		else if (c == '*')
		{
			float op2 = stack.back();
			stack.pop_back();
			float op1 = stack.back();
			stack.pop_back();
			stack.push_back(op1*op2);
		}
		else if (c == '<')
		{
			float op2 = stack.back();
			stack.pop_back();
			float op1 = stack.back();
			stack.pop_back();
			if (token[1] == '=')
			{
				stack.push_back((op1 <= op2)?1.0f:0.0f);
			}
			else
			{
				stack.push_back((op1 < op2)?1.0f:0.0f);
			}
		}
		else if (c == '>')
		{
			float op2 = stack.back();
			stack.pop_back();
			float op1 = stack.back();
			stack.pop_back();
			if (token[1] == '=')
			{
				stack.push_back((op1 >= op2)?1.0f:0.0f);
			}
			else
			{
				stack.push_back((op1 > op2)?1.0f:0.0f);
			}
		}
		else if (c == '=')
		{
			float op2 = stack.back();
			stack.pop_back();
			float op1 = stack.back();
			stack.pop_back();
			stack.push_back((abs(op1 - op2) < 0.000001f)?1.0f:0.0f);
		}
		else if (c == '/')
		{
			float op2 = stack.back();
			stack.pop_back();
			float op1 = stack.back();
			stack.pop_back();
			if (abs(op2) > 0.00000001f)
			{
				stack.push_back(op1/op2);
			}
			else
			{
				stack.push_back(0.0f);
			}
		}
		else if (c == '~' && token[1] == '%')
		{
			float op = stack.back();
			stack.pop_back();
			if (total == -1.0f)
			{
				op = 0.0f;
			}
			stack.push_back(op*total/100.0f);
		}
		else if (c == '~' && token[1] == '-')
		{
			float op = stack.back();
			stack.pop_back();
			stack.push_back(-op);
		}
		else if (is_argument_chain(c))
		{
			argumentStack.push_back((int)stack.size());
		}
		else if (is_function(c))
		{
			VFString identifier = (char*)(&token[1]);
			unsigned int stackPos = argumentStack.back();
			argumentStack.pop_back();
			TVFVector<float> arguments;
			while (stack.size() > stackPos)
			{
				float arg = stack.back();
				stack.pop_back();
				arguments.push_back(arg);
			}
			float value = 0.0f;
			if (identifier.compare("sin") == 0 && arguments.size() == 1)
			{
				value = sin(arguments.at(0));
			}
			else if (identifier.compare("cos") == 0 && arguments.size() == 1)
			{
				value = cos(arguments.at(0));
			}
			else if (identifier.compare("tan") == 0 && arguments.size() == 1)
			{
				value = tan(arguments.at(0));
			}
			else if (identifier.compare("abs") == 0 && arguments.size() == 1)
			{
				value = abs(arguments.at(0));
			}
			else if (identifier.compare("floor") == 0 && arguments.size() == 1)
			{
				value = floorf(arguments.at(0));
			}
			else if (identifier.compare("sqrt") == 0 && arguments.size() == 1 && arguments.at(0) > 0.0f)
			{
				value = sqrt(arguments.at(0));
			}
			else if (identifier.compare("pow") == 0 && arguments.size() == 2)
			{
				value = pow(arguments.at(1), arguments.at(0));
			}
			else if (identifier.compare("mod") == 0 && arguments.size() == 2)
			{
				int o1 = (int)arguments.at(1);
				int o2 = (int)arguments.at(0);
				if (o2 != 0)
				{
					value = (float)(o1 % o2);
				}
				else
				{
					value = 0.0f;
				}
			}
			else if (identifier.compare("max") == 0 && arguments.size() == 2)
			{
				value = std::max(arguments.at(0), arguments.at(1));
			}
			else if (identifier.compare("min") == 0 && arguments.size() == 2)
			{
				value = std::min(arguments.at(0), arguments.at(1));
			}
			stack.push_back(value);
		}
		else if (is_ident(c))
		{
			VFString identifier = token;
			float value = 0.f;

			bool valueFound = from_string<float>(value, identifier, std::dec);

			if (!valueFound)
			{
				if (identifier.compare("rnd") == 0)
				{
					valueFound = true;
					Matrix matrix = Matrix_identity();
					updateScope(currentNode, &matrix);
					Vector origin = Matrix_multiplyVector(matrix, Vector_withValues(0.0, 0.0, 0.0));
					unsigned int x = (unsigned int)(1000*(origin.x + rootOffsetX));
					unsigned int y = (unsigned int)(1000*(origin.y + rootOffsetY));
					unsigned int z = (unsigned int)(1000*(origin.z + rootOffsetZ));
					int random = CWhiteNoise::getValue(x, y, z);
					value = ((float)random)/RAND_MAX;
				}
				else if (identifier.compare("width") == 0)
				{
					valueFound = true;
					value = currentNode->scope.sx;
				}
				else if (identifier.compare("height") == 0)
				{
					valueFound = true;
					value = currentNode->scope.sy;
				}
				else if (identifier.compare("depth") == 0)
				{
					valueFound = true;
					value = currentNode->scope.sz;
				}
				else if (identifier.compare("xpos") == 0)
				{
					valueFound = true;
					Matrix matrix = Matrix_identity();
					updateScope(currentNode, &matrix);
					Vector origin = Matrix_multiplyVector(matrix, Vector_withValues(0.0, 0.0, 0.0));
					value = origin.x;
				}
				else if (identifier.compare("ypos") == 0)
				{
					valueFound = true;
					Matrix matrix = Matrix_identity();
					updateScope(currentNode, &matrix);
					Vector origin = Matrix_multiplyVector(matrix, Vector_withValues(0.0, 0.0, 0.0));
					value = origin.y;
				}
				else if (identifier.compare("zpos") == 0)
				{
					valueFound = true;
					Matrix matrix = Matrix_identity();
					updateScope(currentNode, &matrix);
					Vector origin = Matrix_multiplyVector(matrix, Vector_withValues(0.0, 0.0, 0.0));
					value = origin.z;
				}
				else if (identifier.compare("vertical") == 0)
				{
					valueFound = true;

					Matrix matrix = Matrix_identity();;
					updateScope(currentNode, &matrix);
					Vector origin = Vector_withValues(0.0, 0.0, 0.0);
					Vector up = Vector_withValues(0.0, 1.0, 0.0);
					origin = Matrix_multiplyVector(matrix, origin);
					Vector scopeUp = Matrix_multiplyVector(matrix, up);
					scopeUp = Vector_subtract(scopeUp, origin);
					value = Vector_dot(up, scopeUp);
				}
				else if (identifier.compare("occluded") == 0)
				{
					valueFound = true;
					value = (isOccluded(currentNode, "full", false))?1.0f:0.0f;
				}
				else if (identifier.compare("occludedfront") == 0 ||
						 identifier.compare("occludedback") == 0)
				{
					valueFound = true;
					CNode testNode;
					testNode.parent = currentNode;
					testNode.scope.sx = currentNode->scope.sx;
					testNode.scope.sy = currentNode->scope.sy;
					testNode.scope.sz = 0.0001f;
					if (identifier.compare("occludedfront") == 0)
					{
						testNode.scope.z = -1.0f;
					}
					else if (identifier.compare("occludedback") == 0)
					{
						testNode.scope.z = currentNode->scope.sz + 1.0f;
					}
					value = (isOccluded(&testNode, "", false))?1.0f:0.0f;
				}
			}

			if (!valueFound)
			{
				TVFMap<VFString, float>::iterator iparm = parameters.find(identifier);
				valueFound = iparm != parameters.end();
				if (valueFound)
				{
					value = iparm->second;
				}
			}

			if (!valueFound)
			{
				CNode* node = startNode;
				while (!valueFound && node != NULL)
				{
					TVFMap<VFString, VFString>::iterator i = node->variables.find(identifier);
					valueFound = i != node->variables.end();
					node = node->parent;
					if (valueFound)
					{
						const VFString& rvalue = i->second;
						bool recursive = rvalue.size() > identifier.size() && rvalue.find(identifier + ",") != VFString::npos;
						if (!recursive)
						{
							value = eval(i->second, total, node);
						}
					}
				}
			}

			stack.push_back(value);
		}
	}

	// The top of the stack contains the final value
	return stack.back();

}

void scopeToSnapSegments(CNode* node, Vector segments[4][2], int dir, float rad)
{
	// Get transformation matrix for current scope
	Matrix m = Matrix_identity();
	updateScope(node, &m);

	// Create segments for scope depending on direction
	if (dir == 0)
	{
		segments[0][0] = Vector_withValues(-rad, 0.0f,	0.0f);
		segments[0][1] = Vector_withValues(node->scope.sx + rad, 0.0f, 0.0f);
		segments[1][0] = Vector_withValues(-rad, 0.0f,	node->scope.sz);
		segments[1][1] = Vector_withValues(node->scope.sx + rad, 0.0f, node->scope.sz);
		segments[2][0] = Vector_withValues(-rad, node->scope.sy, 0.0f);
		segments[2][1] = Vector_withValues(node->scope.sx + rad, node->scope.sy, 0.0f);
		segments[3][0] = Vector_withValues(-rad, node->scope.sy,	node->scope.sz);
		segments[3][1] = Vector_withValues(node->scope.sx + rad, node->scope.sy, node->scope.sz);
	}
	else if (dir == 1)
	{
		segments[0][0] = Vector_withValues(0, -rad,	0);
		segments[0][1] = Vector_withValues(0, node->scope.sy + rad,	0);
		segments[1][0] = Vector_withValues(node->scope.sx, -rad, 0);
		segments[1][1] = Vector_withValues(node->scope.sx, node->scope.sy + rad, 0);
		segments[2][0] = Vector_withValues(node->scope.sx, -rad, node->scope.sz);
		segments[2][1] = Vector_withValues(node->scope.sx, node->scope.sy + rad, node->scope.sz);
		segments[3][0] = Vector_withValues(0.0f, -rad, node->scope.sz);
		segments[3][1] = Vector_withValues(0.0f, node->scope.sy + rad, node->scope.sz);
	}
	else if (dir == 2)
	{
		segments[0][0] = Vector_withValues(0.0f, 0.0f, -rad);
		segments[0][1] = Vector_withValues(0.0f, 0.0f, node->scope.sz + rad);
		segments[1][0] = Vector_withValues(node->scope.sx, 0.0f, -rad);
		segments[1][1] = Vector_withValues(node->scope.sx, 0.0f, node->scope.sz + rad);
		segments[2][0] = Vector_withValues(0.0f, node->scope.sy, -rad);
		segments[2][1] = Vector_withValues(0.0f, node->scope.sy, node->scope.sz + rad);
		segments[3][0] = Vector_withValues(node->scope.sx, node->scope.sy, -rad);
		segments[3][1] = Vector_withValues(node->scope.sx, node->scope.sy, node->scope.sz + rad);
	}

	// Transform segments using matrix
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 2; j++)
		{
			segments[i][j] = Matrix_multiplyVector(m, segments[i][j]);
		}
}


void scopeToSnapPlane(CNode* node, SnapPlane& p, int dir, float rad)
{
	// Get transformation matrix for current scope
	Matrix m = Matrix_identity();
	updateScope(node, &m);

	// Create plane endpoints depending on direction
	if (dir == 0)
	{
		p.s1p1 = Vector_withValues(0, -rad,					-rad);
		p.s1p2 = Vector_withValues(0, -rad,					node->scope.sz + rad);
		p.s2p1 = Vector_withValues(0, node->scope.sy + rad, -rad);
		p.s2p2 = Vector_withValues(0, node->scope.sy + rad, node->scope.sz + rad);
	}
	else if (dir == 1)
	{
		p.s1p1 = Vector_withValues(-rad,					0, -rad);
		p.s1p2 = Vector_withValues(node->scope.sx + rad,	0, -rad);
		p.s2p1 = Vector_withValues(-rad,					0, node->scope.sz + rad);
		p.s2p2 = Vector_withValues(node->scope.sx + rad,	0, node->scope.sz + rad);
	}
	else if (dir == 2)
	{
		p.s1p1 = Vector_withValues(-rad,					-rad,				  0);
		p.s1p2 = Vector_withValues(node->scope.sx + rad,	-rad,				  0);
		p.s2p1 = Vector_withValues(-rad,					node->scope.sy + rad, 0);
		p.s2p2 = Vector_withValues(node->scope.sx + rad,	node->scope.sy + rad, 0);
	}

	// Transform points using matrix
	p.s1p1 = Matrix_multiplyVector(m, p.s1p1);
	p.s1p2 = Matrix_multiplyVector(m, p.s1p2);
	p.s2p1 = Matrix_multiplyVector(m, p.s2p1);
	p.s2p2 = Matrix_multiplyVector(m, p.s2p2);
}

bool intersectScopes(CNode* nodeA, CNode* nodeB, bool testPartial)
{
	// Get transformation matrices for nodes
	Matrix mA = Matrix_identity();
	Matrix mB = Matrix_identity();
	updateScope(nodeA, &mA, true);
	updateScope(nodeB, &mB, true);

	// Compute coordinates for the eight corners in nodeA's scope
	float asx = nodeA->scope.sx;
	float asy = nodeA->scope.sy;
	float asz = nodeA->scope.sz;
	if (nodeA->ngon.active || nodeA->ngonLoft.active)
	{
		float minD = std::min(asx, asz);
		asx = minD;
		asz = minD;
	}
	Vector vA[8];
	vA[0] = Vector_withValues(0, 0, 0);
	vA[1] = Vector_withValues(asx, 0, 0);
	vA[2] = Vector_withValues(asx, 0, asz);
	vA[3] = Vector_withValues(0, 0, asz);
	vA[4] = Vector_withValues(0, asy, 0);
	vA[5] = Vector_withValues(asx, asy, 0);
	vA[6] = Vector_withValues(asx, asy, asz);
	vA[7] = Vector_withValues(0, asy, asz);
	for (int i = 0; i < 8; i++)
	{
		vA[i] = Matrix_multiplyVector(mA, vA[i]);
	}

	// Compute coordinates for the eight corners in nodeB's scope
	float bsx = nodeB->scope.sx;
	float bsy = nodeB->scope.sy;
	float bsz = nodeB->scope.sz;
	if (nodeB->ngon.active || nodeB->ngonLoft.active)
	{
		float minD = std::min(bsx, bsz);
		bsx = minD;
		bsz = minD;
	}
	Vector vB[8];
	vB[0] = Vector_withValues(0, 0, 0);
	vB[1] = Vector_withValues(bsx, 0, 0);
	vB[2] = Vector_withValues(bsx, 0, bsz);
	vB[3] = Vector_withValues(0, 0, bsz);
	vB[4] = Vector_withValues(0, bsy, 0);
	vB[5] = Vector_withValues(bsx, bsy, 0);
	vB[6] = Vector_withValues(bsx, bsy, bsz);
	vB[7] = Vector_withValues(0, bsy, bsz);
	for (int i = 0; i < 8; i++)
	{
		vB[i] = Matrix_multiplyVector(mB, vB[i]);
	}


	// Peform separating axis test over the two boxes
	int boxFaces[6][3] =
	{
		{0, 3, 1}, // bottom
		{4, 5, 7}, // top
		{0, 1, 4}, // front
		{3, 7, 2}, // back
		{1, 2, 5}, // right
		{0, 4, 3}  // left
	};

	bool outside = false;
	int insideFacesA = 0;
	int insideFacesB = 0;

	// Test A vs B
	for (int f = 0; f < 6 && (!outside || testPartial); f++)
	{
		Vector& pA = vA[boxFaces[f][0]];
		Vector v10 = Vector_withValues(
						 vA[boxFaces[f][1]].x - pA.x,
						 vA[boxFaces[f][1]].y - pA.y,
						 vA[boxFaces[f][1]].z - pA.z);
		Vector v20 = Vector_withValues(
						 vA[boxFaces[f][2]].x - pA.x,
						 vA[boxFaces[f][2]].y - pA.y,
						 vA[boxFaces[f][2]].z - pA.z);
		Vector n = Vector_cross(v20, v10);

		bool alloutside = true;
		int insideCount = 0;
		int i = 0;
		for (; i < 8 && (alloutside || testPartial); i++)
		{
			Vector v = Vector_subtract(vB[i], pA);
			float d = Vector_dot(v, n);
			if (d < 0.0)
			{
				insideCount++;
			}
			alloutside &= d > 0.0;
		}
		if (insideCount == 8)
		{
			insideFacesA++;
		}
		outside |= alloutside;
	}

	// Test B vs A
	for (int f = 0; f < 6 && (!outside || testPartial); f++)
	{
		Vector& p = vB[boxFaces[f][0]];
		Vector v10 = Vector_withValues(
						 vB[boxFaces[f][1]].x - p.x,
						 vB[boxFaces[f][1]].y - p.y,
						 vB[boxFaces[f][1]].z - p.z);
		Vector v20 = Vector_withValues(
						 vB[boxFaces[f][2]].x - p.x,
						 vB[boxFaces[f][2]].y - p.y,
						 vB[boxFaces[f][2]].z - p.z);
		Vector n = Vector_cross(v20, v10);

		bool alloutside = true;
		int insideCount = 0;
		int i = 0;
		for (; i < 8 && (alloutside || testPartial); i++)
		{
			Vector v = Vector_subtract(vA[i], p);
			float d = Vector_dot(v, n);
			if (d < 0.0)
			{
				insideCount++;
			}
			alloutside &= d > 0.0;
		}
		if (insideCount == 8)
		{
			insideFacesB++;
		}
		outside |= alloutside;
	}

	if (outside)
	{
		return false;
	}
	else if (!testPartial)
	{
		return true;
	}
	else
	{
		return (insideFacesA != 6 && insideFacesB != 6);
	}
}

void projectRectPoint(Vector& v, int axis, float rect[2])
{
	if (axis == 0)
	{
		rect[0] = v.y;
		rect[1] = v.z;
	}
	else if (axis == 1)
	{
		rect[0] = v.x;
		rect[1] = v.z;
	}
	else if (axis == 2)
	{
		rect[0] = v.x;
		rect[1] = v.y;
	}
}

bool testForOcclusion(CNode* node, CNode* start, HyVoxel::TVFSet<CNode*>& exclude, VFString& volumeId)
{
	bool occludes = false;
	// Iterate over all children to test for occlusion
	for (TVFVector<CNode*>::iterator i = start->children.begin();
			i != start->children.end() && !occludes; ++i)
	{
		CNode* child = *i;

		// If child is not an occluder it can be skipped
		if (!child->occludes)
		{
			continue;
		}

		// If child is in the exclusion list, skip it
		if (exclude.find(child) != exclude.end())
		{
			continue;
		}

		// See if there is a geometric scope intersection
		if (intersectScopes(node, child, false), false)
		{
			occludes = true;
		}
	}

	// If there was no occlusion, test children
	for (TVFVector<CNode*>::iterator i = start->children.begin();
			i != start->children.end() && !occludes; ++i)
	{
		occludes = testForOcclusion(node, *i, exclude, volumeId);
	}
	return occludes;
}

bool isOccluded(CNode* node, VFString volumeId)
{
	// Collect nodes along the parent chain.
	// These nodes will be excluded from the occlusion tests.
	HyVoxel::TVFSet<CNode*> exclude;
	CNode* top = node->parent;
	while (top != NULL && top->parent != NULL)
	{
		if (top->occludes)
		{
			exclude.insert(top);
		}
		top = top->parent;
	}
	if (top != NULL)
	{
		return testForOcclusion(node, top, exclude, volumeId);
	}
	else
	{
		return false;
	}
}

namespace HyVoxel
{
	namespace Architecture
	{
		// This struct holds the search context for occluder RTree queries
		struct OcclusionIntersectionContext
		{
			CNode* node;
			CNode* occluder;
			bool partial;
		};
	}
}

#if TARGET_OS_IPHONE || defined(__ANDROID_API__)
bool occlusionIntersectionCallback(CNode* occluderNode, void* context)
#else
bool __cdecl occlusionIntersectionCallback(CNode* occluderNode, void* context)
#endif
{
	// The function will run if the AABB of node and occluder intersect
	OcclusionIntersectionContext* ctx = (OcclusionIntersectionContext*)context;

	// Now test for an accurate scope intersection
	if (intersectScopes(ctx->node, occluderNode, ctx->partial))
	{
		ctx->occluder = occluderNode;
		return false;
	}
	return true;
}

inline void translateName(CNode* node, VFString& name)
{
	CNode* localRoot = node;
	bool rootFound = false;
	// Navigate along the parent chain
	while (!rootFound && localRoot != NULL)
	{
		if (localRoot->hasLocalNames)
		{
			// See if there is an alias set up for the name and return it
			HyVoxel::TVFMap<VFString, VFString>::iterator il = localRoot->localNames.find(name);
			if (il != localRoot->localNames.end())
			{
				name = il->second;
				rootFound = true;
			}
		}
		localRoot = localRoot->parent;
	}
}

bool CEvaluator::isOccluded(CNode* node, VFString volumeId, bool partial)
{
	// Translate local aliases first
	translateName(node, volumeId);

	// See if there are occluders registered for that name
	TVFMap<VFString, OccluderIndex*>::iterator i = occluders.find(volumeId);
	if (i == occluders.end())
	{
		return false;
	}

	// Get index tree for occluders
	OccluderIndex* occluderTree = i->second;

	// Compute node's transform matrix
	Matrix matrix = Matrix_identity();
	updateScope(node, &matrix, true);

	// Compute node AABB from node's 8 corner points
	float asx = node->scope.sx;
	float asy = node->scope.sy;
	float asz = node->scope.sz;
	Vector v[8];
	v[0] = Vector_withValues(0, 0, 0);
	v[1] = Vector_withValues(asx, 0, 0);
	v[2] = Vector_withValues(asx, 0, asz);
	v[3] = Vector_withValues(0, 0, asz);
	v[4] = Vector_withValues(0, asy, 0);
	v[5] = Vector_withValues(asx, asy, 0);
	v[6] = Vector_withValues(asx, asy, asz);
	v[7] = Vector_withValues(0, asy, asz);
	float minp[3];
	float maxp[3];
	for (int i = 0; i < 8; i++)
	{
		Vector vp = Matrix_multiplyVector(matrix, v[i]);
		if (i == 0)
		{
			minp[0] = vp.x;
			minp[1] = vp.y;
			minp[2] = vp.z;
			maxp[0] = vp.x;
			maxp[1] = vp.y;
			maxp[2] = vp.z;
		}
		else
		{
			minp[0] = std::min(minp[0], vp.x);
			minp[1] = std::min(minp[1], vp.y);
			minp[2] = std::min(minp[2], vp.z);
			maxp[0] = std::max(maxp[0], vp.x);
			maxp[1] = std::max(maxp[1], vp.y);
			maxp[2] = std::max(maxp[2], vp.z);
		}
	}

	// Search the RTree
	OcclusionIntersectionContext ctx;
	ctx.node = node;
	ctx.occluder = NULL;
	ctx.partial = partial;
	occluderTree->Search(minp, maxp, occlusionIntersectionCallback, (void*)&ctx);

	if (ctx.occluder != NULL)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void CEvaluator::recurse(VFString rule)
{
	// Locate rule first, if it is not found, exit
	TVFMap<VFString, TVFVector<CRule> >::iterator irv = currentModule->rules.find(rule);
	if (irv == currentModule->rules.end())
	{
		return;
	}

	// Find list of definitions for the rule
	TVFVector<CRule>& rules = irv->second;

	// Iterate over definitions to find which one matches best
	CRule* match = NULL;
	for (TVFVector<CRule>::iterator i = rules.begin();
			i != rules.end() && match == NULL; ++i)
	{
		// Get reference to rule
		CRule& r = *i;

		// Evaluate conditionals for the rule
		bool fail = false;
		for (TVFVector<VFString>::iterator icond = r.conditionals.begin();
				icond != r.conditionals.end() && !fail; ++icond)
		{
			fail = abs(eval(*icond)) < 0.0000001f;
		}

		// See if conditionals not failed
		if (!fail)
		{
			// Perform occlusion tests
			for (unsigned int i = 0; i < r.occlusionScopes.size() && !fail; i++)
			{
				// Get occlusion test parameters
				VFString scope = r.occlusionScopes.at(i);
				VFString volumeId = r.occlusionVolumeIds.at(i);
				float margin = eval(r.occlusionMargins.at(i));

				// Create a test node and update its scope
				CNode testNode;
				testNode.parent = currentNode;
				memset(&testNode.scope, 0, sizeof(testNode.scope));
				bool partial = scope.compare("partial") == 0;
				if (partial || scope.compare("full") == 0)
				{
					testNode.scope.x = std::min(0.0f, -margin);
					testNode.scope.y = std::min(0.0f, -margin);
					testNode.scope.z = std::min(0.0f, -margin);
					testNode.scope.sx = std::max(0.0f, currentNode->scope.sx + margin);
					testNode.scope.sy = std::max(0.0f, currentNode->scope.sy + margin);
					testNode.scope.sz = std::max(0.0f, currentNode->scope.sz + margin);
				}
				else if (scope.compare("front") == 0)
				{
					testNode.scope.sx = currentNode->scope.sx;
					testNode.scope.sy = currentNode->scope.sy;
					testNode.scope.sz = 0.0001f;
					testNode.scope.z = -margin;
				}
				else if (scope.compare("back") == 0)
				{
					testNode.scope.sx = currentNode->scope.sx;
					testNode.scope.sy = currentNode->scope.sy;
					testNode.scope.sz = 0.0001f;
					testNode.scope.z = currentNode->scope.sz + margin;
				}

				// See if the test node is occluded
				fail = !isOccluded(&testNode, volumeId, partial);
			}
		}
		if (!fail)
		{
			match = &r;
		}
		else
		{
			match = NULL;
		}
	}

	// If no rule was matched, exit the function
	if (match == NULL)
	{
		return;
	}

	// Iterate over all actions in the rule and execute them
	for (TVFVector<CAction>::iterator i = match->actions.begin();
			i != match->actions.end(); ++i)
	{
		CAction& a = *i;

		if (debugger != NULL)
		{
			switch (debugAction)
			{
			case  DebugAction::DEBUG_RUN:
				break;
			case DebugAction::DEBUG_INSTANCE:
				if (a.type == ACTION_INSTANCE)
				{
					debugAction = debugger->debug(a.debugInfo.filename.c_str(), a.debugInfo.line, a.debugInfo.column, this);
				}
				break;
			default:
				debugAction = debugger->debug(a.debugInfo.filename.c_str(), a.debugInfo.line, a.debugInfo.column, this);
				break;
			}
		}

		switch (a.type)
		{
		case ACTION_PUSH:
		{
			// Create a new node and insert it in the execution tree.
			// This now becomes the current node
			CNode* newNode = VF_NEW CNode();
			currentNode->children.push_back(newNode);
			newNode->type = currentNode->type;
			newNode->parent = currentNode;
			newNode->scope.type = currentNode->scope.type;
			newNode->scope.sx = currentNode->scope.sx;
			newNode->scope.sy = currentNode->scope.sy;
			newNode->scope.sz = currentNode->scope.sz;
			currentNode = newNode;
		}
		break;
		case ACTION_POP:
		{
			// Execution node is moved to the parent
			currentNode = currentNode->parent;
		}
		break;
		case ACTION_MOVE:
		{
			// Move scope
			currentNode->scope.dx += eval(a.attributes["x"], currentNode->scope.sx);
			currentNode->scope.dy += eval(a.attributes["y"], currentNode->scope.sy);
			currentNode->scope.dz += eval(a.attributes["z"], currentNode->scope.sz);
		}
		break;
		case ACTION_ROTATE:
		{
			// Rotate scope
			currentNode->scope.ax += eval(a.attributes["x"]);
			currentNode->scope.ay += eval(a.attributes["y"]);
			currentNode->scope.az += eval(a.attributes["z"]);
		}
		break;
		case ACTION_VAR:
		{
			// Define a variable
			VFString constant = a.attributes["const"];
			if (constant.compare("true") != 0)
			{
				currentNode->variables[a.attributes["id"]] = a.attributes["value"];
			}
			else
			{
				float value = eval(a.attributes["value"], 1.f, currentNode);
				char str[255];
				sprintf_s(str, 255, "%f", value);
				currentNode->variables[a.attributes["id"]] = str;
			}
		}
		break;
		case ACTION_SNAP:
		{
			// Define a snap tree
			VFString id = a.attributes["id"];
			float r = eval(a.attributes["radius"]);
			VFString dir = a.attributes["dir"];
			int d;
			if (dir.compare("x") == 0)
			{
				d = 0;
			}
			else if (dir.compare("y") == 0)
			{
				d = 1;
			}
			else //if (dir.compare("z") == 0)
			{
				d = 2;
			}

			TVFVector<SnapPlane>& planes = snapPlanes[id];
			planes.push_back(SnapPlane());
			SnapPlane& plane = planes.back();
			scopeToSnapPlane(currentNode, plane, d, r);
		}
		break;
		case ACTION_SCALE:
		{
			// Scale scope
			float sx = currentNode->scope.sx;
			float sy = currentNode->scope.sy;
			float sz = currentNode->scope.sz;
			float nsx = std::max(0.0f, eval(a.attributes["x"], currentNode->scope.sx));
			float nsy = std::max(0.0f, eval(a.attributes["y"], currentNode->scope.sy));
			float nsz = std::max(0.0f, eval(a.attributes["z"], currentNode->scope.sz));
			currentNode->scope.sx = nsx;
			currentNode->scope.sy = nsy;
			currentNode->scope.sz = nsz;
			if (sx > 0.0f)
			{
				currentNode->scope.dx += currentNode->scope.px*(sx - currentNode->scope.sx)/sx;
			}
			if (sy > 0.0f)
			{
				currentNode->scope.dy += currentNode->scope.py*(sy - currentNode->scope.sy)/sy;
			}
			if (sz > 0.0f)
			{
				currentNode->scope.dz += currentNode->scope.pz*(sz - currentNode->scope.sz)/sz;
			}
		}
		break;
		case ACTION_LOCAL:
		{
			// Create an alias for the specified name
			VFString occluderId = a.attributes["volumeid"];
			char localName[512];
			sprintf_s(localName, 512, "%s %d", occluderId.c_str(), guid++);
			currentNode->localNames[occluderId] = VFString(localName);
			currentNode->hasLocalNames = true;
		}
		break;
		case ACTION_OCCLUDE:
		{
			// Registers current node as an occluder
			currentNode->occludes = true;
			VFString occluderId = a.attributes["volumeid"];

			Matrix matrix = Matrix_identity();
			updateScope(currentNode, &matrix);

			float asx = currentNode->scope.sx;
			float asy = currentNode->scope.sy;
			float asz = currentNode->scope.sz;
			Vector v[8];
			v[0] = Vector_withValues(0, 0, 0);
			v[1] = Vector_withValues(asx, 0, 0);
			v[2] = Vector_withValues(asx, 0, asz);
			v[3] = Vector_withValues(0, 0, asz);
			v[4] = Vector_withValues(0, asy, 0);
			v[5] = Vector_withValues(asx, asy, 0);
			v[6] = Vector_withValues(asx, asy, asz);
			v[7] = Vector_withValues(0, asy, asz);
			float minp[3];
			float maxp[3];
			for (int i = 0; i < 8; i++)
			{
				Vector vp = Matrix_multiplyVector(matrix, v[i]);
				if (i == 0)
				{
					minp[0] = vp.x;
					minp[1] = vp.y;
					minp[2] = vp.z;
					maxp[0] = vp.x;
					maxp[1] = vp.y;
					maxp[2] = vp.z;
				}
				else
				{
					minp[0] = std::min(minp[0], vp.x);
					minp[1] = std::min(minp[1], vp.y);
					minp[2] = std::min(minp[2], vp.z);
					maxp[0] = std::max(maxp[0], vp.x);
					maxp[1] = std::max(maxp[1], vp.y);
					maxp[2] = std::max(maxp[2], vp.z);
				}
			}

			translateName(currentNode, occluderId);

			OccluderIndex* occluderTree;
			TVFMap<VFString, OccluderIndex*>::iterator i = occluders.find(occluderId);
			if (i == occluders.end())
			{
				occluderTree = VF_NEW OccluderIndex();
				occluders[occluderId] = occluderTree;
			}
			else
			{
				occluderTree = i->second;
			}
			occluderTree->Insert(minp, maxp, currentNode);
		}
		break;
		case ACTION_INSTANCE:
		{
			// Declares an instance for the current scope
			currentNode->actions.push_back(a);

			// Get material to be used
			CNode* materialNode = currentNode;
			while (materialNode->material.empty() && materialNode->parent != NULL)
			{
				materialNode = materialNode->parent;
			}

			// Compute instance size and position matrix
			Vector size = Vector_withValues(
							  currentNode->scope.sx,
							  currentNode->scope.sy,
							  currentNode->scope.sz);
			Matrix matrix = Matrix_identity();
			updateScope(currentNode, &matrix);

			// Call callback or delegate interface to notify existence of instance
			if (instancer != NULL)
			{
				instancer(a.attributes["geometry"], matrix, size, currentNode->scope.type, materialNode->material);
			}
			else if (instanceCreator != NULL)
			{
				instanceCreator->createInstance(a.attributes["geometry"], matrix, size, currentNode->scope.type, materialNode->material);
			}
		}
		break;
		case ACTION_MATERIAL:
		{
			// Set material
			currentNode->material = a.attributes["id"];
		}
		break;
		case ACTION_CALL:
		{
			// Invokes a submodule
			VFString submodule = a.attributes["name"];
			int generation = 0;
			int maxGenerations = (int)eval(a.attributes["generations"], 1);

			// Compute for how many generations the same module has run
			CNode* parent = currentNode->parent;
			while (parent != NULL)
			{
				if (parent->module.compare(submodule) == 0)
				{
					generation++;
				}
				parent = parent->parent;
			}

			// See if generation is under the allowed limit
			if (generation < maxGenerations || (maxGenerations == -1 && generation < 1000))
			{
				// Create a new node
				CNode* newNode = VF_NEW CNode();
				currentNode->children.push_back(newNode);
				currentNode->module = submodule;
				newNode->type = currentNode->type;
				newNode->parent = currentNode;
				newNode->scope.sx = currentNode->scope.sx;
				newNode->scope.sy = currentNode->scope.sy;
				newNode->scope.sz = currentNode->scope.sz;

				// Start a new generation on the new node
				currentNode->newGeneration(newNode, submodule);
			}
		}
		break;
		case ACTION_DEFER:
		{
			// Record deferred execution for a submodule
			VFString submodule = a.attributes["name"];

			// Create new node
			CNode* newNode = VF_NEW CNode();
			currentNode->children.push_back(newNode);
			currentNode->module = submodule;
			newNode->type = currentNode->type;
			newNode->parent = currentNode;
			newNode->scope.sx = currentNode->scope.sx;
			newNode->scope.sy = currentNode->scope.sy;
			newNode->scope.sz = currentNode->scope.sz;

			// Set a deferred request for the node
			currentNode->defer(newNode, submodule);
		}
		break;
		case ACTION_ALIGN:
		{
			// Set the scope to be aligned to the specified axis
			VFString dir = a.attributes["dir"];
			if (dir.compare("x") == 0)
			{
				currentNode->scope.xaligned = true;
			}
			else if (dir.compare("y") == 0)
			{
				currentNode->scope.yaligned = true;
			}
			else if (dir.compare("z") == 0)
			{
				currentNode->scope.zaligned = true;
			}
		}
		break;
		case ACTION_CENTER_PIVOT:
		{
			// Centers scope pivot on the specified axis or plane
			VFString dir = a.attributes["dir"];
			if (dir.compare("xyz") == 0)
			{
				currentNode->scope.px = currentNode->scope.sx/2.0f;
				currentNode->scope.py = currentNode->scope.sy/2.0f;
				currentNode->scope.pz = currentNode->scope.sz/2.0f;
			}
			else if (dir.compare("x") == 0)
			{
				currentNode->scope.px = currentNode->scope.sx/2.0f;
			}
			else if (dir.compare("y") == 0)
			{
				currentNode->scope.py = currentNode->scope.sy/2.0f;
			}
			else if (dir.compare("z") == 0)
			{
				currentNode->scope.pz = currentNode->scope.sz/2.0f;
			}
			else if (dir.compare("xy") == 0)
			{
				currentNode->scope.px = currentNode->scope.sx/2.0f;
				currentNode->scope.py = currentNode->scope.sy/2.0f;
			}
			else if (dir.compare("xz") == 0)
			{
				currentNode->scope.px = currentNode->scope.sx/2.0f;
				currentNode->scope.pz = currentNode->scope.sz/2.0f;
			}
			else if (dir.compare("yz") == 0)
			{
				currentNode->scope.py = currentNode->scope.sy/2.0f;
				currentNode->scope.pz = currentNode->scope.sz/2.0f;
			}
		}
		break;
		case ACTION_RESET_PIVOT:
		{
			// Resets pivot point
			currentNode->scope.px = 0.0f;
			currentNode->scope.py = 0.0f;
			currentNode->scope.pz = 0.0f;
		}
		break;
		case ACTION_BOX:
		{
			currentNode->box.active = true;
		}
		break;
		case ACTION_LOFT_BOX:
		{
			// Sets a box as active primitive for the node
			currentNode->boxLoft.active = true;
			currentNode->boxLoft.positions = a.attributes["positions"];
			currentNode->boxLoft.widths = a.attributes["widths"];
			currentNode->boxLoft.depths = a.attributes["depths"];
		}
		break;
		case ACTION_LOFT_NGON:
		{
			// Sets a loft ngon as active primitive for the node
			currentNode->ngonLoft.active = true;
			currentNode->ngonLoft.positions = a.attributes["positions"];
			currentNode->ngonLoft.radii = a.attributes["radii"];
			currentNode->ngonLoft.sides = (int)eval(a.attributes["sides"]);

			float sx = currentNode->scope.sx;
			float sz = currentNode->scope.sz;
			currentNode->scope.sx = std::min(sx, sz);
			currentNode->scope.sz = std::min(sx, sz);
			if (sx > 0.0f)
			{
				currentNode->scope.dx += currentNode->scope.px*(sx - currentNode->scope.sx)/sx;
			}
			if (sz > 0.0f)
			{
				currentNode->scope.dz += currentNode->scope.pz*(sz - currentNode->scope.sz)/sz;
			}
		}
		break;
		case ACTION_NGON:
		{
			// Sets a ngon as active primitive for the node
			currentNode->ngon.active = true;
			currentNode->ngon.sides = (int)eval(a.attributes["sides"]);

			float sx = currentNode->scope.sx;
			float sz = currentNode->scope.sz;
			currentNode->scope.sx = std::min(sx, sz);
			currentNode->scope.sz = std::min(sx, sz);
			if (sx > 0.0f)
			{
				currentNode->scope.dx += currentNode->scope.px*(sx - currentNode->scope.sx)/sx;
			}
			if (sz > 0.0f)
			{
				currentNode->scope.dz += currentNode->scope.pz*(sz - currentNode->scope.sz)/sz;
			}
		}
		break;
		case ACTION_PRISM:
		{
			// Sets the scope type to prism
			currentNode->scope.type = SCOPE_PRISM;
		}
		break;
		case ACTION_SELECT:
		{
			// Select a primitive feature and turn it into subscopes
			VFString type = a.attributes["type"];
			VFString submodule = a.attributes["module"];

			// Get list of submodules to be applied
			TVFVector<VFString> modules;
			std::istringstream issModules(submodule);
			std::copy(std::istream_iterator<VFString>(issModules),
					  std::istream_iterator<VFString>(),
					  std::back_inserter<TVFVector<VFString> >(modules));

			if (currentNode->ngonLoft.active)
			{
				// Get positions from the primitive
				TVFVector<VFString> positions;
				std::istringstream issPositions(currentNode->ngonLoft.positions);
				std::copy(std::istream_iterator<VFString>(issPositions),
						  std::istream_iterator<VFString>(),
						  std::back_inserter<TVFVector<VFString> >(positions));

				// Get radii from the primitive
				TVFVector<VFString> radii;
				std::istringstream issRadii(currentNode->ngonLoft.radii);
				std::copy(std::istream_iterator<VFString>(issRadii),
						  std::istream_iterator<VFString>(),
						  std::back_inserter<TVFVector<VFString> >(radii));

				if (type.compare("sides") == 0)
				{
					// Select sides
					int sides = currentNode->ngonLoft.sides;
					float slice = 360.0f/sides;
					float angle = 0.0f;
					float rot = 90.f + slice/2.0f;
					float radius = std::min(currentNode->scope.sx, currentNode->scope.sz)/2.0f;
					float xc = currentNode->scope.sx/2.0f;
					float zc = currentNode->scope.sz/2.0f;

					// For each side, create subscope and apply corresponding submodule
					for (int i = 0; i < sides; i++)
					{
						float x1 = xc + radius*cos((float)(angle*piover180f));
						float z1 = zc + radius*sin((float)(angle*piover180f));
						float x2 = xc + radius*cos((float)((angle + slice)*piover180f));
						float z2 = zc + radius*sin((float)((angle + slice)*piover180f));
						float length = sqrtf((x1 - x2)*(x1 - x2) + (z1 - z2)*(z1 - z2));

						float shift = 0.0f;
						float shiftX = 0.0f;
						float shiftZ = 0.0f;
						float prevXS = currentNode->scope.sx;
						float prevZS = radius;

						for (unsigned int i = 0; i < positions.size(); i++)
						{
							float pos = eval(positions.at(i), currentNode->scope.sy);
							float depth = eval(radii.at(i), radius);
							float width = length*depth/radius;

							float dx = (prevZS - depth)*sin((float)(0.5f*slice*piover180f));
							float dz = (prevZS - depth)*cos((float)(0.5f*slice*piover180f));
							float beta = (dz != 0.0f)?atanf((pos - shift)/dz):((float)(90.f*piover180f));
							float cosBeta = cos(beta);

							// Create new node
							CNode* newNode = VF_NEW CNode();
							currentNode->children.push_back(newNode);
							newNode->type = NODE_FACE;
							newNode->parent = currentNode;
							newNode->scope.x = x1;
							newNode->scope.y = shift;
							newNode->scope.z = z1;
							float alpha = acos((x2 - x1)/length)/piover180f;
							if (x1 <= x2 && z1 <= z2)
							{
								newNode->scope.ay = alpha;
							}
							else if (x1 > x2 && z1 <= z2)
							{
								newNode->scope.ay = alpha;
							}
							else if (x1 > x2 && z1 > z2)
							{
								newNode->scope.ay = 360 - alpha;
							}
							else
							{
								newNode->scope.ay = 360 - alpha;
							}

							newNode->scope.sx = width;
							if (fabs(cosBeta) > 0.00001f)
							{
								newNode->scope.sy = fabs(dz/cosBeta);
							}
							else
							{
								newNode->scope.sy = pos;
							}
							newNode->scope.sz = 0.0f;

							// Create subnode
							CNode* subNode = VF_NEW CNode();
							newNode->children.push_back(subNode);
							subNode->parent = newNode;
							subNode->scope.x = shiftX + dx;
							subNode->scope.z = shiftZ;
							subNode->scope.sx = newNode->scope.sx;
							subNode->scope.sy = newNode->scope.sy;
							subNode->scope.sz = newNode->scope.sz;
							subNode->scope.ax = beta/piover180f - 90.0f;

							// Request generation for subnode
							newNode->newGeneration(subNode, modules.at(i % modules.size()));

							shift = pos;
							shiftX += dx;
							shiftZ += dz;
							prevXS = width;
							prevZS = depth;
						}

						angle += slice;
						rot += (x2 - x1)*(180 + slice)/fabs(x2 - x1);
					}

				}
				else if (type.compare("prism_sides") == 0)
				{
					int sides = currentNode->ngonLoft.sides;
					float slice = 360.0f/sides;
					float angle = 0.0f;
					float rot = 90.f + slice/2.0f;
					float radius = std::min(currentNode->scope.sx, currentNode->scope.sz)/2.0f;
					float xc = currentNode->scope.sx/2.0f;
					float zc = currentNode->scope.sz/2.0f;
					for (int i = 0; i < sides; i++)
					{
						float x1 = xc + radius*cos(angle*piover180f);
						float z1 = zc + radius*sin(angle*piover180f);
						float x2 = xc + radius*cos((angle + slice)*piover180f);
						float z2 = zc + radius*sin((angle + slice)*piover180f);
						float length = sqrtf((x1 - x2)*(x1 - x2) + (z1 - z2)*(z1 - z2));

						float shift = 0.0f;
						float shiftX = 0.0f;
						float shiftZ = 0.0f;
						float prevXS = currentNode->scope.sx;
						float prevZS = radius;

						for (unsigned int i = 0; i < positions.size(); i++)
						{
							float pos = eval(positions.at(i), currentNode->scope.sy);
							float depth = eval(radii.at(i), radius);
							float width = length*depth/radius;

							float dx = (prevZS - depth)*sin(0.5f*slice*piover180f);
							float dz = (prevZS - depth)*cos(0.5f*slice*piover180f);
							float alpha = acos((x2 - x1)/length)/piover180f;
							float beta = (dz != 0.0f)?atanf((pos - shift)/dz):(90.0f*piover180f);
							float cosBeta = cos(beta);

							// Left prism

							CNode* newNode1 = VF_NEW CNode();
							currentNode->children.push_back(newNode1);
							newNode1->type = NODE_FACE;
							newNode1->parent = currentNode;
							newNode1->scope.x = x1;
							newNode1->scope.y = shift;
							newNode1->scope.z = z1;

							if (x1 <= x2 && z1 <= z2)
							{
								newNode1->scope.ay = alpha;
							}
							else if (x1 > x2 && z1 <= z2)
							{
								newNode1->scope.ay = alpha;
							}
							else if (x1 > x2 && z1 > z2)
							{
								newNode1->scope.ay = 360 - alpha;
							}
							else
							{
								newNode1->scope.ay = 360 - alpha;
							}

							newNode1->scope.sx = dx;
							if (fabs(cosBeta) > 0.00001f)
							{
								newNode1->scope.sy = fabs(dz/cosBeta);
							}
							else
							{
								newNode1->scope.sy = pos;
							}
							newNode1->scope.sz = 0.0f;

							CNode* subNode1 = VF_NEW CNode();
							newNode1->children.push_back(subNode1);
							subNode1->parent = newNode1;
							subNode1->scope.type = SCOPE_PRISM_LEFT;
							subNode1->scope.x = shiftX + dx;
							subNode1->scope.z = shiftZ;
							subNode1->scope.sx = dz;
							subNode1->scope.sy = newNode1->scope.sy;
							subNode1->scope.sz = dx;
							subNode1->scope.ax = beta/piover180f - 90.0f;
							subNode1->scope.ay = 90.0f;

							newNode1->newGeneration(subNode1, modules.at(i % modules.size()));


							// Right prism
							CNode* newNode2 = VF_NEW CNode();
							currentNode->children.push_back(newNode2);
							newNode2->type = NODE_FACE;
							newNode2->parent = currentNode;
							newNode2->scope.x = x1;
							newNode2->scope.y = shift;
							newNode2->scope.z = z1;

							if (x1 <= x2 && z1 <= z2)
							{
								newNode2->scope.ay = alpha;
							}
							else if (x1 > x2 && z1 <= z2)
							{
								newNode2->scope.ay = alpha;
							}
							else if (x1 > x2 && z1 > z2)
							{
								newNode2->scope.ay = 360 - alpha;
							}
							else
							{
								newNode2->scope.ay = 360 - alpha;
							}

							newNode2->scope.sx = dx;
							if (fabs(cosBeta) > 0.00001f)
							{
								newNode2->scope.sy = fabs(dz/cosBeta);
							}
							else
							{
								newNode2->scope.sy = pos;
							}
							newNode2->scope.sz = 0.0f;

							CNode* subNode2 = VF_NEW CNode();
							newNode2->children.push_back(subNode2);
							subNode2->parent = newNode2;
							subNode2->scope.type = SCOPE_PRISM_RIGHT;
							subNode2->scope.x = shiftX + 2*dx + width;
							subNode2->scope.z = shiftZ;
							subNode2->scope.sx = dz;
							subNode2->scope.sy = newNode2->scope.sy;
							subNode2->scope.sz = dx;
							subNode2->scope.ax = beta/piover180f - 90.0f;
							subNode2->scope.ay = -270.0f;

							newNode2->newGeneration(subNode2, modules.at(i % modules.size()));

							shift = pos;
							shiftX += dx;
							shiftZ += dz;
							prevXS = width;
							prevZS = depth;
						}

						angle += slice;
						rot += (x2 - x1)*(180 + slice)/fabs(x2 - x1);
					}
				}
			}
			else if (currentNode->boxLoft.active)
			{
				// Box loft primitive is active
				TVFVector<VFString> positions;
				std::istringstream issPositions(currentNode->boxLoft.positions);
				std::copy(std::istream_iterator<VFString>(issPositions),
						  std::istream_iterator<VFString>(),
						  std::back_inserter<TVFVector<VFString> >(positions));

				TVFVector<VFString> widths;
				std::istringstream issWidths(currentNode->boxLoft.widths);
				std::copy(std::istream_iterator<VFString>(issWidths),
						  std::istream_iterator<VFString>(),
						  std::back_inserter<TVFVector<VFString> >(widths));

				TVFVector<VFString> depths;
				std::istringstream issDepths(currentNode->boxLoft.depths);
				std::copy(std::istream_iterator<VFString>(issDepths),
						  std::istream_iterator<VFString>(),
						  std::back_inserter<TVFVector<VFString> >(depths));

				if (positions.size() == widths.size() &&
						positions.size() == depths.size())
				{
					if (	type.compare("face_x") == 0 ||
							type.compare("front") == 0 ||
							type.compare("back") == 0 ||
							type.compare("sides") == 0)
					{
						float shift = 0.0f;
						float shiftZ = 0.0f;
						float prevXS = currentNode->scope.sx;
						float prevZS = currentNode->scope.sz;
						for (unsigned int i = 0; i < positions.size(); i++)
						{
							float pos = eval(positions.at(i), currentNode->scope.sy);
							float width = eval(widths.at(i), currentNode->scope.sx);
							float depth = eval(depths.at(i), currentNode->scope.sz);

							float dx = 0.5f*(prevXS - width);
							float dz = 0.5f*(prevZS - depth);
							float alpha = (dz != 0.0f)?atanf((pos - shift)/dz):(90.0f*piover180f);
							float cosAlpha = cos(alpha);

							// face 1
							if (type.compare("back") != 0)
							{
								CNode* newNodeA = VF_NEW CNode();
								currentNode->children.push_back(newNodeA);
								newNodeA->type = NODE_FACE;
								newNodeA->parent = currentNode;
								newNodeA->scope.x = dx;
								newNodeA->scope.y = shift;
								newNodeA->scope.z = shiftZ;
								newNodeA->scope.sx = currentNode->scope.sx - 2*dx;
								if (fabs(cosAlpha) > 0.00001f)
								{
									newNodeA->scope.sy = fabs(dz/cosAlpha);
								}
								else
								{
									newNodeA->scope.sy = pos - shift;
								}
								newNodeA->scope.sz = 0.0f;
								if (alpha >= 0)
								{
									newNodeA->scope.ax = alpha/piover180f - 90.0f;
								}
								else
								{
									newNodeA->scope.ax = 90.0f + alpha/piover180f;
								}
								currentNode->newGeneration(newNodeA, modules.at(i % modules.size()));
							}

							// face 2
							if (type.compare("front") != 0)
							{
								CNode* newNodeB = VF_NEW CNode();
								currentNode->children.push_back(newNodeB);
								newNodeB->type = NODE_FACE;
								newNodeB->parent = currentNode;
								newNodeB->scope.x = currentNode->scope.sx - dx;
								newNodeB->scope.y = shift;
								newNodeB->scope.z = currentNode->scope.sz - shiftZ;
								if (alpha >= 0.0f)
								{
									newNodeB->scope.ax = 90.0f - alpha/piover180f;
								}
								else
								{
									newNodeB->scope.ax = 90.0f - alpha/piover180f + 180.0f;
								}
								newNodeB->scope.ay = 180.0f;
								newNodeB->scope.sx = currentNode->scope.sx - 2*dx;
								if (fabs(cosAlpha) > 0.00001f)
								{
									newNodeB->scope.sy = fabs(dz/cosAlpha);
								}
								else
								{
									newNodeB->scope.sy = pos;
								}
								newNodeB->scope.sz = 0.0f;
								currentNode->newGeneration(newNodeB, modules.at(i % modules.size()));
							}

							shift = pos;
							shiftZ += dz;
							prevZS = depth;
						}
					}
					if (type.compare("prism_x") == 0 || type.compare("prism_sides") == 0)
					{
						float shift = 0.0f;
						float shiftX = 0.0f;
						float shiftZ = 0.0f;
						float prevXS = currentNode->scope.sx;
						float prevZS = currentNode->scope.sz;
						for (unsigned int i = 0; i < positions.size(); i++)
						{
							float pos = eval(positions.at(i), currentNode->scope.sy);
							float width = eval(widths.at(i), currentNode->scope.sx);
							float depth = eval(depths.at(i), currentNode->scope.sz);

							float dx = 0.5f*(prevXS - width);
							float dz = 0.5f*(prevZS - depth);
							float alpha = (dz != 0.0f)?atanf((pos - shift)/dz):(90.0f*piover180f);
							float sinAlpha = sin(alpha);

							// left prism 1
							CNode* newNodeA1 = VF_NEW CNode();
							currentNode->children.push_back(newNodeA1);
							newNodeA1->type = NODE_FACE;
							newNodeA1->parent = currentNode;
							newNodeA1->scope.type = SCOPE_PRISM_LEFT;

							newNodeA1->scope.x = shiftX + dx;
							newNodeA1->scope.y = shift;
							newNodeA1->scope.z = shiftZ;
							newNodeA1->scope.sz = dx;
							newNodeA1->scope.sy = pos;
							newNodeA1->scope.ay = 90.0f;
							if (fabs(sinAlpha) > 0.00001f)
							{
								newNodeA1->scope.sy = fabs((pos - shift)/sinAlpha);
							}
							else
							{
								newNodeA1->scope.sy = dz;
							}
							if (alpha >= 0)
							{
								newNodeA1->scope.ax = alpha/piover180f - 90.0f;
							}
							else
							{
								newNodeA1->scope.ax = 90.0f + alpha/piover180f;
							}

							currentNode->newGeneration(newNodeA1, modules.at(i % modules.size()));

							// left prism 2
							CNode* newNodeA2 = VF_NEW CNode();
							currentNode->children.push_back(newNodeA2);
							newNodeA2->type = NODE_FACE;
							newNodeA2->parent = currentNode;
							newNodeA2->scope.type = SCOPE_PRISM_LEFT;

							newNodeA2->scope.x = currentNode->scope.sx - (shiftX + dx);
							newNodeA2->scope.y = shift;
							newNodeA2->scope.z = currentNode->scope.sz - shiftZ;
							newNodeA2->scope.sz = dx;
							newNodeA2->scope.sy = pos;
							newNodeA2->scope.ay = -90.0f;
							if (fabs(sinAlpha) > 0.00001f)
							{
								newNodeA2->scope.sy = fabs((pos - shift)/sinAlpha);
							}
							else
							{
								newNodeA2->scope.sy = dz;
							}
							if (alpha >= 0)
							{
								newNodeA2->scope.ax = -alpha/piover180f + 90.0f;
							}
							else
							{
								newNodeA2->scope.ax = 90.0f - alpha/piover180f + 180.0f;
							}

							currentNode->newGeneration(newNodeA2, modules.at(i % modules.size()));

							// right prism 1
							CNode* newNodeB1 = VF_NEW CNode();
							currentNode->children.push_back(newNodeB1);
							newNodeB1->type = NODE_FACE;
							newNodeB1->parent = currentNode;
							newNodeB1->scope.type = SCOPE_PRISM_RIGHT;

							newNodeB1->scope.x = currentNode->scope.sx - (shiftX);
							newNodeB1->scope.y = shift;
							newNodeB1->scope.z = shiftZ;
							newNodeB1->scope.sz = dx;
							newNodeB1->scope.sy = pos;
							newNodeB1->scope.ay = 90.0f;
							if (alpha >= 0)
							{
								newNodeB1->scope.ax = alpha/piover180f - 90.0f;
							}
							else
							{
								newNodeB1->scope.ax = 90.0f + alpha/piover180f;
							}
							if (fabs(sinAlpha) > 0.00001f)
							{
								newNodeB1->scope.sy = fabs((pos - shift)/sinAlpha);
							}
							else
							{
								newNodeB1->scope.sy = dz;
							}

							currentNode->newGeneration(newNodeB1, modules.at(i % modules.size()));

							// right prism 2
							CNode* newNodeB2 = VF_NEW CNode();
							currentNode->children.push_back(newNodeB2);
							newNodeB2->type = NODE_FACE;
							newNodeB2->parent = currentNode;
							newNodeB2->scope.type = SCOPE_PRISM_RIGHT;

							newNodeB2->scope.x = shiftX;
							newNodeB2->scope.y = shift;
							newNodeB2->scope.z = currentNode->scope.sz - shiftZ;
							newNodeB2->scope.sz = dx;
							newNodeB2->scope.sy = pos;
							newNodeB2->scope.ay = -90.0f;
							if (fabs(sinAlpha) > 0.00001f)
							{
								newNodeB2->scope.sy = fabs((pos - shift)/sinAlpha);
							}
							else
							{
								newNodeB2->scope.sy = dz;
							}
							if (alpha > 0)
							{
								newNodeB2->scope.ax = -alpha/piover180f + 90.0f;
							}
							else
							{
								newNodeB2->scope.ax = 90.0f + alpha/piover180f;
							}

							currentNode->newGeneration(newNodeB2, modules.at(i % modules.size()));

							shift = pos;
							shiftX += dx;
							shiftZ += dz;
							prevXS = width;
							prevZS = depth;
						}
					}
					if (type.compare("face_z") == 0 ||
							type.compare("left") == 0 ||
							type.compare("right") == 0 ||
							type.compare("sides") == 0)
					{
						float shift = 0.0f;
						float shiftX = 0.0f;
						float prevXS = currentNode->scope.sx;
						float prevZS = currentNode->scope.sz;
						for (unsigned int i = 0; i < positions.size(); i++)
						{
							float pos = eval(positions.at(i), currentNode->scope.sy);
							float width = eval(widths.at(i), currentNode->scope.sx);
							float depth = eval(depths.at(i), currentNode->scope.sz);

							float dx = 0.5f*(prevXS - width);
							float dz = 0.5f*(prevZS - depth);
							float alpha = (dx != 0.0f)?atanf((pos - shift)/dx):(90.0f*piover180f);
							float cosAlpha = cos(alpha);

							// face 1
							if (type.compare("left") != 0)
							{
								CNode* newNodeA = VF_NEW CNode();
								currentNode->children.push_back(newNodeA);
								newNodeA->type = NODE_FACE;
								newNodeA->parent = currentNode;
								newNodeA->scope.x = shiftX;
								newNodeA->scope.y = shift;
								newNodeA->scope.z = currentNode->scope.sz - dz;
								newNodeA->scope.sx = currentNode->scope.sz - 2*dz;
								if (fabs(cosAlpha) > 0.00001f)
								{
									newNodeA->scope.sy = fabs(dx/cosAlpha);
								}
								else
								{
									newNodeA->scope.sy = pos;
								}
								newNodeA->scope.sz = 0.0f;
								if (alpha > 0)
								{
									newNodeA->scope.az = -alpha/piover180f + 90.0f;
								}
								else
								{
									newNodeA->scope.az = 90.0f + alpha/piover180f;
								}
								newNodeA->scope.ay = -90.0f;
								currentNode->newGeneration(newNodeA, modules.at(i % modules.size()));
							}

							// face2
							if (type.compare("right") != 0)
							{
								CNode* newNodeB = VF_NEW CNode();
								currentNode->children.push_back(newNodeB);
								newNodeB->type = NODE_FACE;
								newNodeB->parent = currentNode;
								newNodeB->scope.x = currentNode->scope.sx - shiftX;
								newNodeB->scope.y = shift;
								newNodeB->scope.z = dz;
								newNodeB->scope.sx = currentNode->scope.sz - 2*dz;
								if (fabs(cosAlpha) > 0.00001f)
								{
									newNodeB->scope.sy = fabs(dx/cosAlpha);
								}
								else
								{
									newNodeB->scope.sy = pos;
								}
								newNodeB->scope.sz = 0.0f;
								newNodeB->scope.az = alpha/piover180f - 90.0f;
								newNodeB->scope.ay = 90.0f;
								currentNode->newGeneration(newNodeB, modules.at(i % modules.size()));
							}

							shift = pos;
							shiftX += dx;
							prevXS = width;
						}
					}
					if (type.compare("prism_z") == 0 || type.compare("prism_sides") == 0)
					{
						float shift = 0.0f;
						float shiftX = 0.0f;
						float shiftZ = 0.0f;
						float prevXS = currentNode->scope.sx;
						float prevZS = currentNode->scope.sz;
						for (unsigned int i = 0; i < positions.size(); i++)
						{
							float pos = eval(positions.at(i), currentNode->scope.sy);
							float width = eval(widths.at(i), currentNode->scope.sx);
							float depth = eval(depths.at(i), currentNode->scope.sz);

							float dx = 0.5f*(prevXS - width);
							float dz = 0.5f*(prevZS - depth);
							float alpha = (dx != 0.0f)?atanf((pos - shift)/dx):(90.0f*piover180f);
							float cosAlpha = cos(alpha);

							// left prism 1
							CNode* newNodeA1 = VF_NEW CNode();
							currentNode->children.push_back(newNodeA1);
							newNodeA1->type = NODE_FACE;
							newNodeA1->parent = currentNode;
							newNodeA1->scope.type = SCOPE_PRISM_LEFT;

							newNodeA1->scope.x = shiftX;
							newNodeA1->scope.y = shift;
							newNodeA1->scope.z = currentNode->scope.sz - (shiftZ + dz);
							newNodeA1->scope.sz = dz;
							if (fabs(cosAlpha) > 0.00001f)
							{
								newNodeA1->scope.sy = fabs(dx/cosAlpha);
							}
							else
							{
								newNodeA1->scope.sy = (pos - shift);
							}
							newNodeA1->scope.az = 90 - alpha/piover180f;

							currentNode->newGeneration(newNodeA1, modules.at(i % modules.size()));

							// right prism 1
							CNode* newNodeB1 = VF_NEW CNode();
							currentNode->children.push_back(newNodeB1);
							newNodeB1->type = NODE_FACE;
							newNodeB1->parent = currentNode;
							newNodeB1->scope.type = SCOPE_PRISM_RIGHT;

							newNodeB1->scope.x = shiftX;
							newNodeB1->scope.y = shift;
							newNodeB1->scope.z = shiftZ;
							newNodeB1->scope.sz = dz;
							if (fabs(cosAlpha) > 0.00001f)
							{
								newNodeB1->scope.sy = fabs(dx/cosAlpha);
							}
							else
							{
								newNodeB1->scope.sy = (pos - shift);
							}
							newNodeB1->scope.az = 90 - alpha/piover180f;

							currentNode->newGeneration(newNodeB1, modules.at(i % modules.size()));

							// left prism 2
							CNode* newNodeA2 = VF_NEW CNode();
							currentNode->children.push_back(newNodeA2);
							newNodeA2->type = NODE_FACE;
							newNodeA2->parent = currentNode;
							newNodeA2->scope.type = SCOPE_PRISM_LEFT;

							newNodeA2->scope.x = currentNode->scope.sx - shiftX;
							newNodeA2->scope.y = shift;
							newNodeA2->scope.z = (shiftZ + dz);
							newNodeA2->scope.sz = dz;
							if (fabs(cosAlpha) > 0.00001f)
							{
								newNodeA2->scope.sy = fabs(dx/cosAlpha);
							}
							else
							{
								newNodeA2->scope.sy = (pos - shift);
							}
							newNodeA2->scope.az = alpha/piover180f - 90.0f;
							newNodeA2->scope.ay = 180.0f;

							currentNode->newGeneration(newNodeA2, modules.at(i % modules.size()));

							// right prism 2
							CNode* newNodeB2 = VF_NEW CNode();
							currentNode->children.push_back(newNodeB2);
							newNodeB2->type = NODE_FACE;
							newNodeB2->parent = currentNode;
							newNodeB2->scope.type = SCOPE_PRISM_RIGHT;

							newNodeB2->scope.x = currentNode->scope.sx - shiftX;
							newNodeB2->scope.y = shift;
							newNodeB2->scope.z = currentNode->scope.sz - (shiftZ);
							newNodeB2->scope.sz = dz;
							if (fabs(cosAlpha) > 0.00001f)
							{
								newNodeB2->scope.sy = fabs(dx/cosAlpha);
							}
							else
							{
								newNodeB2->scope.sy = (pos - shift);
							}
							newNodeB2->scope.az = alpha/piover180f - 90.0f;
							newNodeB2->scope.ay = 180.0f;

							currentNode->newGeneration(newNodeB2, modules.at(i % modules.size()));

							shift = pos;
							shiftX += dx;
							shiftZ += dz;
							prevXS = width;
							prevZS = depth;
						}
					}
				}
			}
			else if (currentNode->box.active)
			{
				// Primitive is a box
				if (type.compare("face_x") == 0 ||
						type.compare("front") == 0 ||
						type.compare("back") == 0 ||
						type.compare("sides") == 0)
				{
					// face 1
					if (type.compare("back") != 0)
					{
						CNode* newNodeA = VF_NEW CNode();
						currentNode->children.push_back(newNodeA);
						newNodeA->type = NODE_FACE;
						newNodeA->parent = currentNode;
						newNodeA->scope.x = 0.0f;
						newNodeA->scope.y = 0.0f;
						newNodeA->scope.z = 0.0f;
						newNodeA->scope.sx = currentNode->scope.sx;
						newNodeA->scope.sy = currentNode->scope.sy;
						newNodeA->scope.sz = 0.0f;
						currentNode->newGeneration(newNodeA, submodule);
					}

					// face 2
					if (type.compare("front") != 0)
					{
						CNode* newNodeB = VF_NEW CNode();
						currentNode->children.push_back(newNodeB);
						newNodeB->type = NODE_FACE;
						newNodeB->parent = currentNode;
						newNodeB->scope.x = currentNode->scope.sx;
						newNodeB->scope.y = 0.0f;
						newNodeB->scope.z = currentNode->scope.sz;
						newNodeB->scope.ay = 180.0f;
						newNodeB->scope.sx = currentNode->scope.sx;
						newNodeB->scope.sy = currentNode->scope.sy;
						newNodeB->scope.sz = 0.0f;
						currentNode->newGeneration(newNodeB, submodule);
					}
				}
				if (type.compare("face_z") == 0 ||
						type.compare("left") == 0 ||
						type.compare("right") == 0 ||
						type.compare("sides") == 0)
				{
					// face 1
					if (type.compare("left") != 0)
					{
						CNode* newNodeA = VF_NEW CNode();
						currentNode->children.push_back(newNodeA);
						newNodeA->type = NODE_FACE;
						newNodeA->parent = currentNode;
						newNodeA->scope.x = 0.0f;
						newNodeA->scope.y = 0.0f;
						newNodeA->scope.z = currentNode->scope.sz;
						newNodeA->scope.sx = currentNode->scope.sz;
						newNodeA->scope.sy = currentNode->scope.sy;
						newNodeA->scope.sz = 0.0f;
						newNodeA->scope.ay = -90.0f;
						currentNode->newGeneration(newNodeA, submodule);
					}

					// face 2
					if (type.compare("right") != 0)
					{
						CNode* newNodeB = VF_NEW CNode();
						currentNode->children.push_back(newNodeB);
						newNodeB->type = NODE_FACE;
						newNodeB->parent = currentNode;
						newNodeB->scope.x = currentNode->scope.sx;
						newNodeB->scope.y = 0.0f;
						newNodeB->scope.z = 0.0f;
						newNodeB->scope.ay = 90.0f;
						newNodeB->scope.sx = currentNode->scope.sz;
						newNodeB->scope.sy = currentNode->scope.sy;
						newNodeB->scope.sz = 0.0f;
						currentNode->newGeneration(newNodeB, submodule);
					}
				}
			}
			if (currentNode->ngon.active)
			{
				// Primitive is an ngon
				int sides = currentNode->ngon.sides;
				float slice = 360.0f/sides;
				float angle = 0.0f;
				float rot = 90.f + slice/2.0f;
				float radius = std::min(currentNode->scope.sx, currentNode->scope.sz)/2.0f;
				float xc = currentNode->scope.sx/2.0f;
				float zc = currentNode->scope.sz/2.0f;
				for (int i = 0; i < sides; i++)
				{
					float x1 = xc + radius*cos(angle*piover180f);
					float z1 = zc + radius*sin(angle*piover180f);
					float x2 = xc + radius*cos((angle + slice)*piover180f);
					float z2 = zc + radius*sin((angle + slice)*piover180f);
					float length = sqrtf((x1 - x2)*(x1 - x2) + (z1 - z2)*(z1 - z2));

					CNode* newNode = VF_NEW CNode();
					currentNode->children.push_back(newNode);
					newNode->type = NODE_FACE;
					newNode->parent = currentNode;
					newNode->scope.x = x1;
					newNode->scope.y = 0.0f;
					newNode->scope.z = z1;

					float alpha = acos((x2 - x1)/length)/piover180f;
					if (x1 <= x2 && z1 <= z2)
					{
						newNode->scope.ay = alpha;
					}
					else if (x1 > x2 && z1 <= z2)
					{
						newNode->scope.ay = alpha;
					}
					else if (x1 > x2 && z1 > z2)
					{
						newNode->scope.ay = 360 - alpha;
					}
					else
					{
						newNode->scope.ay = 360 - alpha;
					}

					newNode->scope.sx = length;
					newNode->scope.sy = currentNode->scope.sy;
					newNode->scope.sz = 0.0f;

					currentNode->newGeneration(newNode, modules.at(i % modules.size()));

					angle += slice;
					rot += (x2 - x1)*(180 + slice)/fabs(x2 - x1);
				}
			}
			if (currentNode->scope.type == SCOPE_PRISM)
			{
				if (type.compare("sides") == 0)
				{
					// Face 1
					CNode* newNodeA = VF_NEW CNode();
					currentNode->children.push_back(newNodeA);
					newNodeA->type = NODE_FACE;
					newNodeA->parent = currentNode;
					newNodeA->scope.x = 0.0f;
					newNodeA->scope.y = 0.0f;
					newNodeA->scope.z = 0.0f;
					newNodeA->scope.sx = currentNode->scope.sx;
					newNodeA->scope.sy = sqrt(
											 currentNode->scope.sy*currentNode->scope.sy +
											 currentNode->scope.sz*currentNode->scope.sz/4.0f);
					newNodeA->scope.sz = 0.0f;
					float alpha = acosf(currentNode->scope.sy/newNodeA->scope.sy)/piover180f;
					newNodeA->scope.ax = -alpha;
					currentNode->newGeneration(newNodeA, submodule);

					// Face 2
					CNode* newNodeB = VF_NEW CNode();
					currentNode->children.push_back(newNodeB);
					newNodeB->type = NODE_FACE;
					newNodeB->parent = currentNode;
					newNodeB->scope.x = currentNode->scope.sx;
					newNodeB->scope.y = 0.0f;
					newNodeB->scope.z = currentNode->scope.sz;
					newNodeB->scope.ax = alpha;
					newNodeB->scope.ay = 180.0f;
					newNodeB->scope.sx = currentNode->scope.sx;
					newNodeB->scope.sy = newNodeA->scope.sy;
					newNodeB->scope.sz = 0.0f;
					currentNode->newGeneration(newNodeB, submodule);

				}
				if (type.compare("front") == 0)
				{
					CNode* newNode = VF_NEW CNode();
					currentNode->children.push_back(newNode);
					newNode->type = NODE_FACE;
					newNode->parent = currentNode;
					newNode->scope = currentNode->scope;
					newNode->scope.sx = 0.0f;
					newNode->scope.dx = 0.0f;
					newNode->scope.dy = 0.0f;
					newNode->scope.dz = 0.0f;
					newNode->scope.ax = 0.0f;
					newNode->scope.ay = 0.0f;
					newNode->scope.az = 0.0f;
					currentNode->newGeneration(newNode, submodule);
				}
				else if (type.compare("back") == 0)
				{
					CNode* newNode = VF_NEW CNode();
					currentNode->children.push_back(newNode);
					newNode->type = NODE_FACE;
					newNode->parent = currentNode;
					newNode->scope = currentNode->scope;
					newNode->scope.dx = 0.0f;
					newNode->scope.dy = 0.0f;
					newNode->scope.dz = 0.0f;
					newNode->scope.ax = 0.0f;
					newNode->scope.ay = 0.0f;
					newNode->scope.az = 0.0f;
					newNode->scope.x = currentNode->scope.sx;
					newNode->scope.z = currentNode->scope.sz;
					newNode->scope.sx = 0.0f;
					newNode->scope.ay = 180.0f;
					currentNode->newGeneration(newNode, submodule);
				}
			}
		}
		break;
		case ACTION_DIVIDE:
		{
			// Divide current scope into subscopes
			VFString type = a.attributes["type"];
			if (type.compare("prism") == 0)
			{
				// It is prism division
				// Get attributes
				VFString position = a.attributes["position"];
				VFString center = a.attributes["center"];
				VFString top = a.attributes["top"];
				VFString left = a.attributes["left"];
				VFString right = a.attributes["right"];

				float p = eval(position, currentNode->scope.sy);

				if (currentNode->scope.type == SCOPE_PRISM)
				{
					float dz = currentNode->scope.sz*p/(2*currentNode->scope.sy);

					// Center box
					CNode* newNodeA = VF_NEW CNode();
					currentNode->children.push_back(newNodeA);
					newNodeA->type = NODE_FACE;
					newNodeA->parent = currentNode;
					newNodeA->scope.x = 0.0f;
					newNodeA->scope.y = 0.0f;
					newNodeA->scope.z = currentNode->scope.sz - dz;
					newNodeA->scope.sx = currentNode->scope.sz - 2*dz;
					newNodeA->scope.sy = p;
					newNodeA->scope.sz = 0.0f;
					newNodeA->scope.ay = -90.0f;
					currentNode->newGeneration(newNodeA, center);

					// Top prism
					CNode* newNodeB = VF_NEW CNode();
					currentNode->children.push_back(newNodeB);
					newNodeB->type = NODE_FACE;
					newNodeB->parent = currentNode;
					newNodeB->scope.type = SCOPE_PRISM;
					newNodeB->scope.x = 0;
					newNodeB->scope.y = p;
					newNodeB->scope.z = dz;
					newNodeB->scope.sx = 0.0f;
					newNodeB->scope.sy = currentNode->scope.sy - p;
					newNodeB->scope.sz = currentNode->scope.sz - 2*dz;
					currentNode->newGeneration(newNodeB, top);

					// Left prism
					CNode* newNodeC = VF_NEW CNode();
					currentNode->children.push_back(newNodeC);
					newNodeC->type = NODE_FACE;
					newNodeC->parent = currentNode;
					newNodeC->scope.type = SCOPE_PRISM_LEFT;
					newNodeC->scope.x = 0;
					newNodeC->scope.y = 0;
					newNodeC->scope.z = currentNode->scope.sz - dz;
					newNodeC->scope.sx = 0.0f;
					newNodeC->scope.sy = p;
					newNodeC->scope.sz = dz;
					currentNode->newGeneration(newNodeC, left);

					// Right prism
					CNode* newNodeD = VF_NEW CNode();
					currentNode->children.push_back(newNodeD);
					newNodeD->type = NODE_FACE;
					newNodeD->parent = currentNode;
					newNodeD->scope.type = SCOPE_PRISM_RIGHT;
					newNodeD->scope.x = 0;
					newNodeD->scope.y = 0;
					newNodeD->scope.z = 0;
					newNodeD->scope.sx = 0.0f;
					newNodeD->scope.sy = p;
					newNodeD->scope.sz = dz;
					currentNode->newGeneration(newNodeD, right);
				}
				else if (currentNode->scope.type == SCOPE_PRISM_LEFT)
				{
					float dz = currentNode->scope.sz*p/(currentNode->scope.sy);

					// Center box
					CNode* newNodeA = VF_NEW CNode();
					currentNode->children.push_back(newNodeA);
					newNodeA->type = NODE_FACE;
					newNodeA->parent = currentNode;
					newNodeA->scope.x = 0;
					newNodeA->scope.y = 0;
					newNodeA->scope.z = currentNode->scope.sz - dz;
					newNodeA->scope.sx = currentNode->scope.sz - dz;
					newNodeA->scope.sy = p;
					newNodeA->scope.sz = 0.0f;
					newNodeA->scope.ay = -90.0f;
					currentNode->newGeneration(newNodeA, center);

					// Top prism
					CNode* newNodeB = VF_NEW CNode();
					currentNode->children.push_back(newNodeB);
					newNodeB->type = NODE_FACE;
					newNodeB->parent = currentNode;
					newNodeB->scope.type = SCOPE_PRISM_LEFT;
					newNodeB->scope.x = 0;
					newNodeB->scope.y = p;
					newNodeB->scope.z = 0.0f;
					newNodeB->scope.sx = 0.0f;
					newNodeB->scope.sy = currentNode->scope.sy - p;
					newNodeB->scope.sz = currentNode->scope.sz - dz;
					currentNode->newGeneration(newNodeB, top);

					// Left prism
					CNode* newNodeC = VF_NEW CNode();
					currentNode->children.push_back(newNodeC);
					newNodeC->type = NODE_FACE;
					newNodeC->parent = currentNode;
					newNodeC->scope.type = SCOPE_PRISM_LEFT;
					newNodeC->scope.x = 0;
					newNodeC->scope.y = 0;
					newNodeC->scope.z = currentNode->scope.sz - dz;
					newNodeC->scope.sx = 0.0f;
					newNodeC->scope.sy = p;
					newNodeC->scope.sz = dz;
					currentNode->newGeneration(newNodeC, left);
				}
				else if (currentNode->scope.type == SCOPE_PRISM_RIGHT)
				{
					float dz = currentNode->scope.sz*p/(currentNode->scope.sy);

					// Center box
					CNode* newNodeA = VF_NEW CNode();
					currentNode->children.push_back(newNodeA);
					newNodeA->type = NODE_FACE;
					newNodeA->parent = currentNode;
					newNodeA->scope.x = 0;
					newNodeA->scope.y = 0;
					newNodeA->scope.z = currentNode->scope.sz;
					newNodeA->scope.sx = currentNode->scope.sz - dz;
					newNodeA->scope.sy = p;
					newNodeA->scope.sz = 0.0f;
					newNodeA->scope.ay = -90.0f;
					currentNode->newGeneration(newNodeA, center);

					// Top prism
					CNode* newNodeB = VF_NEW CNode();
					currentNode->children.push_back(newNodeB);
					newNodeB->type = NODE_FACE;
					newNodeB->parent = currentNode;
					newNodeB->scope.type = SCOPE_PRISM_RIGHT;
					newNodeB->scope.x = 0;
					newNodeB->scope.y = p;
					newNodeB->scope.z = dz;
					newNodeB->scope.sx = 0.0f;
					newNodeB->scope.sy = currentNode->scope.sy - p;
					newNodeB->scope.sz = currentNode->scope.sz - dz;
					currentNode->newGeneration(newNodeB, top);

					// Left prism
					CNode* newNodeC = VF_NEW CNode();
					currentNode->children.push_back(newNodeC);
					newNodeC->type = NODE_FACE;
					newNodeC->parent = currentNode;
					newNodeC->scope.type = SCOPE_PRISM_RIGHT;
					newNodeC->scope.x = 0;
					newNodeC->scope.y = 0;
					newNodeC->scope.z = 0;
					newNodeC->scope.sx = 0.0f;
					newNodeC->scope.sy = p;
					newNodeC->scope.sz = dz;
					currentNode->newGeneration(newNodeC, left);

				}
			}
			else
			{
				// It is dividing a box

				// Get division segments
				TVFVector<VFString> pattern;
				std::istringstream issPattern(a.attributes["pattern"]);
				std::copy(std::istream_iterator<VFString>(issPattern),
						  std::istream_iterator<VFString>(),
						  std::back_inserter<TVFVector<VFString> >(pattern));

				// Compute initial lengths
				float lengths[128];
				int idx = 0;
				float fixedLength = 0.f;
				for (TVFVector<VFString>::iterator i = pattern.begin();
						i != pattern.end(); ++i)
				{
					lengths[idx] = eval(*i);
					fixedLength += lengths[idx++];
				}

				// Get submodules to be applied
				TVFVector<VFString> modules;
				std::istringstream issModules(a.attributes["modules"]);
				std::copy(std::istream_iterator<VFString>(issModules),
						  std::istream_iterator<VFString>(),
						  std::back_inserter<TVFVector<VFString> >(modules));

				// Number of segments must match number of patterns
				if (pattern.size() == modules.size())
				{
					float shift = 0.f;
					int idx = 0;

					// For each segment:
					for (unsigned int i = 0; i < pattern.size(); i++)
					{
						float lengthValue = lengths[idx++];
						VFString length = pattern.at(i);
						VFString submodule = modules.at(i);

						// Create new node
						CNode* newNode = VF_NEW CNode();
						currentNode->children.push_back(newNode);
						newNode->type = currentNode->type;
						newNode->parent = currentNode;
						float nodeLength = 0.f;

						// Set new node scope according to division axis
						if (type.compare("x") == 0)
						{
							float total = std::max(0.0f, currentNode->scope.sx - fixedLength);
							nodeLength = (abs(lengthValue) < 0.0000001f)?eval(length, total):lengthValue;
							newNode->scope.x = shift;
							newNode->scope.sx = nodeLength;
							newNode->scope.sy = currentNode->scope.sy;
							newNode->scope.sz = currentNode->scope.sz;
						}
						else if (type.compare("y") == 0)
						{
							float total = std::max(0.0f, currentNode->scope.sy - fixedLength);
							nodeLength = (abs(lengthValue) < 0.0000001f)?eval(length, total):lengthValue;
							newNode->scope.y = shift;
							newNode->scope.sx = currentNode->scope.sx;
							newNode->scope.sy = nodeLength;
							newNode->scope.sz = currentNode->scope.sz;
						}
						else if (type.compare("z") == 0)
						{
							float total = std::max(0.0f, currentNode->scope.sz - fixedLength);
							nodeLength = (abs(lengthValue) < 0.0000001f)?eval(length, total):lengthValue;
							newNode->scope.z = shift;
							newNode->scope.sx = currentNode->scope.sx;
							newNode->scope.sy = currentNode->scope.sy;
							newNode->scope.sz = nodeLength;
						}
						shift += nodeLength;

						// Schedule new generation on the new node
						if (nodeLength > 0.0f)
						{
							currentNode->newGeneration(newNode, submodule);
						}
					}
				}
			}
		}
		break;
		case ACTION_REPEAT:
		{
			// Repeat a sequence of modules at regular intervals

			// Get modules
			TVFVector<VFString> modules;
			std::istringstream issModules(a.attributes["module"]);
			std::copy(std::istream_iterator<VFString>(issModules),
					  std::istream_iterator<VFString>(),
					  std::back_inserter<TVFVector<VFString> >(modules));

			// Get snap plane ids
			TVFVector<VFString> snapIds;
			std::istringstream issSnap(a.attributes["snap"]);
			std::copy(std::istream_iterator<VFString>(issSnap),
					  std::istream_iterator<VFString>(),
					  std::back_inserter<TVFVector<VFString> >(snapIds));

			// Compute scope length based on repeat direction
			VFString dir = a.attributes["dir"];
			int d, scopeD;
			float scopeLen = 0.f;
			if (dir.compare("x") == 0)
			{
				d = 0;
				scopeD = 2;
				scopeLen = currentNode->scope.sx;
			}
			else if (dir.compare("y") == 0)
			{
				d = 1;
				scopeD = 2;
				scopeLen = currentNode->scope.sy;
			}
			else if (dir.compare("z") == 0)
			{
				d = 2;
				scopeD = 0;
				scopeLen = currentNode->scope.sz;
			}

			TVFVector<float> subscopes;

			// If no snap is required, insert lengh
			if (snapIds.empty() || snapIds.at(0).compare("*") == 0)
			{
				subscopes.push_back(scopeLen);
			}
			else
			{
				// If snap is required, compute lengh as distance to closest
				// snap plane intersectopm
				SnapPlane plane;
				scopeToSnapPlane(currentNode, plane, d, 0);
				Vector segments[4][2];
				scopeToSnapSegments(currentNode, segments, d, 0);
				TVFVector<float> points;
				findSnapPoints(snapIds.at(0), segments, points);
				if (!points.empty())
				{
					float prev = 0.0f;
					float remainingLength = scopeLen;

					for (TVFVector<float>::iterator i = points.begin();
							i != points.end(); ++i)
					{
						float p = *i;
						float l = p - prev;
						if (l > 0.000001f)
						{
							subscopes.push_back(l);
						}
						remainingLength -= l;
						prev = p;
					}
					if (remainingLength > 0.0f)
					{
						subscopes.push_back(remainingLength);
					}
				}
				else
				{
					subscopes.push_back(scopeLen);
				}
			}

			// Iterate over the found segments
			float scopeShift = 0.0f;
			for (TVFVector<float>::iterator i = subscopes.begin(); i != subscopes.end(); ++i)
			{
				// Get segment's length
				scopeLen = *i;

				// Estimate how many times will be repeated
				float length = eval(a.attributes["length"], scopeLen);
				int times = (length > 0)?(int)(scopeLen/length):0;
				times = std::max(times, 1);
				VFString flex = a.attributes["flex"];
				if (flex.compare("fit") == 0)
				{
					if (times*length < scopeLen)
					{
						times++;
					}
				}

				// Repeat
				float totLen = 0.0f;
				for (int i = 0; i < times; i++)
				{
					// Create node
					CNode* newNode = VF_NEW CNode();
					currentNode->children.push_back(newNode);
					newNode->type = currentNode->type;
					newNode->parent = currentNode;
					newNode->scope.sx = currentNode->scope.sx;
					newNode->scope.sy = currentNode->scope.sy;
					newNode->scope.sz = currentNode->scope.sz;

					// Assign scope dimensions based on repeat axis
					if (flex.compare("flex") == 0)
					{
						if (dir.compare("x") == 0)
						{
							newNode->scope.sx = scopeLen/times;
							newNode->scope.x = i*newNode->scope.sx + scopeShift;
						}
						else if (dir.compare("y") == 0)
						{
							newNode->scope.sy = scopeLen/times;
							newNode->scope.y = i*newNode->scope.sy + scopeShift;
						}
						else if (dir.compare("z") == 0)
						{
							newNode->scope.sz = scopeLen/times;
							newNode->scope.z = i*newNode->scope.sz + scopeShift;
						}
					}
					else if (flex.compare("fit") == 0)
					{
						float segmentLength = std::min(length, scopeLen - totLen);
						if (dir.compare("x") == 0)
						{
							newNode->scope.sx = segmentLength;
							newNode->scope.x = totLen + scopeShift;
						}
						else if (dir.compare("y") == 0)
						{
							newNode->scope.sy = segmentLength;
							newNode->scope.y = totLen + scopeShift;
						}
						else if (dir.compare("z") == 0)
						{
							newNode->scope.sz = segmentLength;
							newNode->scope.z = totLen + scopeShift;
						}
						totLen += length;
					}

					// Schedule the node for execution
					VFString moduleName = modules.at(i % modules.size());
					currentNode->newGeneration(newNode, moduleName);
				}
				scopeShift += scopeLen;
			}
		}
		break;
		default:
		{
			currentNode->actions.push_back(a);
		}
		break;
		}
	}
}

