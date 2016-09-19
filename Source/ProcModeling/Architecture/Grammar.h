/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

#include <map>
#include <vector>
#include <set>
#include <string>
#include <iterator>

#include "Common/Vector.h"
#include "Common/MatrixAlg.h"
#include "Common/RTree.h"
#include "Common/WhiteNoise.h"

namespace HyVoxel
{
	namespace Architecture
	{
		/// Enumerates instructions for the grammar virtual machine
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

		/// Enumerates node types for a grammar program
		enum NodeType
		{
			NODE_ROOT,
			NODE_NULL,
			NODE_VOLUME,
			NODE_FACE,
			NODE_EDGE,
			NODE_POINT
		};

		/// Debug information
		class CDebugInfo
		{
		public:
			VFString filename;
			int line;
			int column;
			CDebugInfo() : filename(""), line(-1), column(-1) {}
		};

		enum DebugAction
		{
			DEBUG_RUN,
			DEBUG_STOP,
			DEBUG_STEP,
			DEBUG_INSTANCE
		};

		class CEvaluator;
		class IDebugger
		{
		public:
			virtual DebugAction debug(const char* filename, int line, int column, CEvaluator* evaluator) = 0;
		};

		/// Defines a single action and its parameters
		class CAction
		{
		public:
			/// Action type
			ActionType type;
			/// A list of parameters for the action
			TVFMap<VFString, VFString> attributes;
			/// Debug information for the action
			CDebugInfo debugInfo;
		};

		/// Defines a scope object
		class CScope
		{
		public:
			CScope() :
				type(SCOPE_BOX),
				x(0), y(0), z(0),
				sx(0), sy(0), sz(0),
				ax(0), ay(0), az(0),
				dx(0), dy(0), dz(0),
				px(0), py(0), pz(0),
				xaligned(false),
				yaligned(false),
				zaligned(false) {}
		public:
			/// Scope origin along X axis
			float x;
			/// Scope origin along Y axis
			float y;
			/// Scope origin along Z axis
			float z;
			/// Scope translation along X axis
			float dx;
			/// Scope translation along Y axis
			float dy;
			/// Scope translation along Z axis
			float dz;
			/// Scope pivot along X axis
			float px;
			/// Scope pivot along Y axis
			float py;
			/// Scope pivot along Z axis
			float pz;
			/// Scope size along X axis
			float sx;
			/// Scope size along Y axis
			float sy;
			/// Scope size along Z axis
			float sz;
			/// Scope rotation angle over X axis
			float ax;
			/// Scope rotation angle over Y axis
			float ay;
			/// Scope rotation angle over Z axis
			float az;
			/// Flags whether the scope is aligned along the X axis
			bool xaligned;
			/// Flags whether the scope is aligned along the Y axis
			bool yaligned;
			/// Flags whether the scope is aligned along the Z axis
			bool zaligned;
			/// Scope type (whether it is a box or a prism, etc)
			ScopeType type;
		};

		/// Defines a box primitive
		class CBoxVolume
		{
		public:
			CBoxVolume() : active(false) {}
			/// Flags whether the primitive is active or not
			bool active;
		};

		/// Defines a ngon primitive
		class CNGonVolume
		{
		public:
			CNGonVolume() : sides(0), radius(0), active(false) {}
			/// Sides in the ngon
			int sides;
			/// Ngon Radius
			float radius;
			/// Flags whether the primitive is active or not
			bool active;
		};

		/// Defines a box loft primitive
		class CBoxLoftVolume
		{
		public:
			CBoxLoftVolume() : active(false) {}
			/// Flags whether the primitive is active or not
			bool active;
			/// Contains the list of positions for the loft
			VFString positions;
			/// Contains the list of widths for the loft
			VFString widths;
			/// Contains the list of depths for the loft
			VFString depths;
		};

		/// Defines a ngon loft primitive
		class CNGonLoftVolume
		{
		public:
			CNGonLoftVolume() : active(false), sides(0) {}
			/// Flags whether the primitive is active or not
			bool active;
			/// Number of sides in the ngon loft
			int sides;
			/// A list of positions for the loft
			VFString positions;
			/// A list of radii for the loft
			VFString radii;
		};

		/// Defines a node in the program tree
		class CNode
		{
		public:
			CNode() : parent(NULL), generation(0), occludes(false), matrixCacheValid(false), hasLocalNames(false), type(0), cachedMatrix() {}
			~CNode();
			/// Node type (see ActionType)
			int type;
			/// Pointer to parent node
			CNode* parent;
			/// List of children nodes
			TVFVector<CNode*> children;
			/// List of child nodes that still require processing
			TVFVector<std::pair<CNode*, VFString> > recursePending;
			/// List of child nodes that requested deferred processing
			TVFVector<std::pair<CNode*, VFString> > deferred;
			/// Variables defined for the node and children
			TVFMap<VFString, VFString> variables;
			/// List of actions to be performed by the node
			TVFVector<CAction> actions;
			/// Scope for the node
			CScope scope;
			/// Box primitive for the node
			CBoxVolume box;
			/// Ngon primitive for the node
			CNGonVolume ngon;
			/// Box loft primitive for the node
			CBoxLoftVolume boxLoft;
			/// Ngon loft primitive for the node
			CNGonLoftVolume ngonLoft;
			/// Counts the tree generation for the node
			int generation;
			/// Stores the ID of the grammar module that originated the node
			VFString module;
			/// Stores the material ID set for the node
			VFString material;
			/// Indicate whether this module may occlude other modules
			bool occludes;
			/// Contain a list of aliases for occlusion and snap identifiers so they can be local to a portion of the execution tree
			TVFMap<VFString, VFString> localNames;
			/// Indicates whether the node has local aliases
			bool hasLocalNames;
			/// To speed up processing, the transformation matrix can be cached at some point
			Matrix cachedMatrix;
			/// Indicates whether the cache transformation maxtrix is still valid
			bool matrixCacheValid;
		public:
			/// Starts a new execution subtree generation from this child, applying the specified grammar module
			void newGeneration(CNode* child, VFString module);
			/// Defers a new execution subtree generation from this child, applying the specified grammar module
			void defer(CNode* child, VFString module);
		};

		class CModule;

		/// Encapsulates a series of actions to be performed over a node. This is equivalent to a "module" statement in the language.
		class CRule
		{
		public:
			/// Defines a variable
			CRule& var(VFString identifier, VFString value, bool constant);
			/// Defines a snap plane
			CRule& snap(VFString dir, VFString id, VFString radius);
			/// Defines a conditional
			CRule& condition(VFString exp);
			/// Defines a scale statement
			CRule& scale(VFString x, VFString y, VFString z);
			/// Defines a translation statement
			CRule& move(VFString x, VFString y, VFString z);
			/// Defines a rotation over the X axis
			CRule& rotate_x(VFString angle);
			/// Defines a rotation over the Y axis
			CRule& rotate_y(VFString angle);
			/// Defines a rotation over the Z axis
			CRule& rotate_z(VFString angle);
			/// Aligns the scope along the specified axis
			CRule& align(VFString dir);
			/// Moves the scope's pivot to the scope's geometric center
			CRule& center();
			/// Moves the scope's pivot to the scope's geometric center along the X axis
			CRule& center_x();
			/// Moves the scope's pivot to the scope's geometric center along the Y axis
			CRule& center_y();
			/// Moves the scope's pivot to the scope's geometric center along the Y axis
			CRule& center_z();
			/// Moves the scope's pivot to the scope's geometric center along the XY plane
			CRule& center_xy();
			/// Moves the scope's pivot to the scope's geometric center along the XZ plane
			CRule& center_xz();
			/// Moves the scope's pivot to the scope's geometric center along the YZ plane
			CRule& center_yz();

			/// Defines a box primitive
			CRule& box();
			/// Defines a ngon primitive
			CRule& ngon(VFString sides);
			/// Defines a prism primitive
			CRule& prism(VFString type);
			/// Defines a box loft primitive
			CRule& loft_box(VFString positions, VFString widths, VFString depths);
			/// Defines a ngon loft primitive
			CRule& loft_ngon(VFString sides, VFString positions, VFString radii);

			/// Invokes a module
			CRule& module(VFString name, VFString generations = "-1");
			/// Deferred invocation for a module
			CRule& defer(VFString name);
			/// Instances a geometry mesh
			CRule& instance(VFString geometry);
			/// Set the node's material
			CRule& material(VFString id);

			/// Select a primitive feature and invoke module within it
			CRule& select(VFString stype, VFString module);

			/// Adds occlusion test
			CRule& if_occluded_rule(VFString scope, VFString volumeId, VFString margin);

			/// Declares the current scope may occlude other scopes
			CRule& occludes(VFString volumeid);

			/// Creates an alias for the specified Id so queries can be localized from this point on
			CRule& local(VFString volumeid);

			/// Divide the current scope along the X axis
			CRule& subdivide_x(VFString pattern, VFString modules);
			/// Divide the current scope along the Y axis
			CRule& subdivide_y(VFString pattern, VFString modules);
			/// Divide the current scope along the Z axis
			CRule& subdivide_z(VFString pattern, VFString modules);
			/// Divide a prism scope
			CRule& subdivide_prism(VFString position, VFString center, VFString top, VFString left, VFString right);
			/// Repeat module along X axis
			CRule& repeat_x(VFString length, VFString module, VFString flex = "flex", VFString snap = "");
			/// Repeat module along Y axis
			CRule& repeat_y(VFString length, VFString module, VFString flex = "flex", VFString snap = "");
			/// Repeat module along Z axis
			CRule& repeat_z(VFString length, VFString module, VFString flex = "flex", VFString snap = "");

			/// Defines a parameter
			CRule& parameter(VFString name, VFString value, VFString module, VFString flex = "flex");

			/// Pushes the current state into the stack
			CRule& push();
			/// Pops the current state into the stack
			CRule& pop();

			/// Used for fluent inteface decalaration. Ends a series of rule declarations and returns a reference to the module.
			CModule& end();
		public:
			/// Pointer to the module that owns the rule
			CModule* parent;
			/// Pointer to the node where the rule is being exectuted
			CNode* context;
			/// List of actions for the rule
			TVFVector<CAction> actions;
			/// List of conditional statements for the rule
			TVFVector<VFString> conditionals;
			/// List of occluding IDs for occlusion tests
			TVFVector<VFString> occlusionVolumeIds;
			/// List of occluding margins for occlusion tests
			TVFVector<VFString> occlusionMargins;
			/// List of occluding modes for occlusion tests (partial, full, front, back)
			TVFVector<VFString> occlusionScopes;
		};

		/// Defines a program
		class CModule
		{
		public:
			/// Declares the starting rule object
			CModule& axiom(VFString name);
			/// Adds a rule
			CModule& add(CRule& rule);
			/// Starts the definition of a new rule object
			CRule& begin(VFString name);
		public:
			/// Maps rules to identifiers. A single identifier can be shared by many rules
			TVFMap<VFString, TVFVector<CRule> > rules;
			/// Stores the base axiom
			VFString baseAxiom;
		public:
			/// Saves the program to a file name
			void saveToFile(const char* filename);
			/// Loads a program from a file name
			void loadFromFile(const char* filename, char* debugfile = NULL);
			void loadFromString(std::string& strBuf);
		};

		/// A grammar is a collection of programs
		class CGrammar
		{
		public:
			/// Starts the definition of a program.
			CModule& define(VFString name);
		public:
			/// A grammar contains several module programs. Each program is identified by a VFString.
			TVFMap<VFString, CModule> modules;
		};

		/// A callback to request the creation of an instance. The IInstanceCreator interface provides an alternative to this callback
		typedef void CreateInstance(
			/// Instance ID
			VFString id,
			/// Translation and rotation matrices for the instance
			Matrix& matrix,
			/// Scale for the instance
			Vector& size,
			/// Scope type
			ScopeType scopeType,
			/// Material ID for the instance
			VFString material);

		/// This interfaces allows the grammar evaluator to delegate creation of instances
		class IInstanceCreator
		{
		public:
			/// Requests the creation of an instance. The IInstanceCreator interface provides an alternative to this callback
			virtual void createInstance(
				/// Instance ID
				VFString id,
				/// Translation and rotation matrices for the instance
				Matrix& matrix,
				/// Scale for the instance
				Vector& size,
				/// Scope type
				ScopeType scopeType,
				/// Material ID for the instance
				VFString material) = 0;
		};

		/// Contains four co-planar 3D points that define a snap plane
		class SnapPlane
		{
		public:
			Vector s1p1;
			Vector s1p2;
			Vector s2p1;
			Vector s2p2;
		};

		/// RTree to accelerate occlussion test.
		typedef RTree<CNode*, float, 3, float, 20> OccluderIndex;

		/// This object can be used to evaluate a grammar
		class CEvaluator
		{
		public:
			CEvaluator();
			/// Root node for the evaluation tree
			CNode root;
		public:
			/// Runs a grammar program
			void run(
				/// A grammar object containing all the programs potentially involved in the evaluation
				CGrammar& grammar,
				/// The ID of the starting module.
				VFString module,
				/// A callback function to create instances.
				CreateInstance* instancer,
				/// A cap on the number of generations the evaluator will go through. Use -1 to let the program terminate freely.
				int maxgeneration = -1);
		public:
			/// Interface to a debugger
			IDebugger* debugger;
			/// Current evaluation node
			CNode* currentNode;
			/// Current grammar being evaluated
			CGrammar* currentGrammar;
			/// Current module being evaluated
			CModule* currentModule;
			/// Evaluates an expresion
			float eval(
				/// Expresion definition
				VFString exp,
				/// A total size to be used as reference for any percentage magnitudes in the expression
				float total = -1.0f,
				/// A node to provide context to the evaluation. Variables declared at this node level will be used first.
				CNode* startNode = NULL);
			/// Evaluates the specified rule and then any other rules it may spawn
			void recurse(VFString rule);
			/// Evaluate any pending rules.
			bool evaluatePending(CNode* node, int maxgeneration = -1);
			/// Evaluate any rules that may have been deferred.
			bool evaluateDeferred(CNode* node, int maxgeneration = -1);

			/// Given four segments in a scope, this function return the intersection points between the segments and all snap planes identifed with the provided ID.
			void findSnapPoints(
				/// Snap plane identifier
				VFString& id,
				/// Scope segments
				Vector segments[4][2],
				/// Resulting intersection points
				TVFVector<float>& points);
			/// Tests whether the provided node is occluded
			bool isOccluded(
				/// Node that is to be tested for occlusion
				CNode* node,
				/// Occluder ID
				VFString volumeId,
				/// Indicates whether partial occlusion should be considered
				bool partial);
		public:
			/// A set of parameters supplied to the evaluator
			TVFMap<VFString, float> parameters;
			/// Stores snap planes according to their IDs
			TVFMap<VFString, TVFVector<SnapPlane> > snapPlanes;
			/// Pointer to the create instance callback
			CreateInstance* instancer;
			/// Pointer to the create instance interface
			IInstanceCreator* instanceCreator;
			/// Store occludes according to their IDs
			TVFMap<VFString, OccluderIndex*> occluders;
			/// Keeps an ever-incrementing counter used to generate unique aliases for local features
			int guid;
			/// Position of the root node in world coordinates. It is used for the rnd command.
			double rootOffsetX;
			double rootOffsetY;
			double rootOffsetZ;
		public:
			DebugAction debugAction;
		};

		void updateScope(CNode* node, HyVoxel::Matrix* matrix, bool useCache = false);
	};

}