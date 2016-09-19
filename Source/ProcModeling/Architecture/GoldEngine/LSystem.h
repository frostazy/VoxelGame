#pragma once

namespace HyVoxel
{
	enum SymbolIndex
	{
		SYM_EOF = 0, // (EOF)
		SYM_ERROR = 1, // (Error)
		SYM_COMMENT = 2, // Comment
		SYM_NEWLINE = 3, // NewLine
		SYM_WHITESPACE = 4, // Whitespace
		SYM_TIMESDIV = 5, // '*/'
		SYM_DIVTIMES = 6, // '/*'
		SYM_DIVDIV = 7, // '//'
		SYM_MINUS = 8, // '-'
		SYM_PERCENT = 9, // '%'
		SYM_LPAREN = 10, // '('
		SYM_RPAREN = 11, // ')'
		SYM_TIMES = 12, // '*'
		SYM_COMMA = 13, // ','
		SYM_DIV = 14, // '/'
		SYM_COLON = 15, // ':'
		SYM_LBRACKET = 16, // '['
		SYM_RBRACKET = 17, // ']'
		SYM_LBRACE = 18, // '{'
		SYM_PIPE = 19, // '|'
		SYM_RBRACE = 20, // '}'
		SYM_PLUS = 21, // '+'
		SYM_LT = 22, // '<'
		SYM_LTEQ = 23, // '<='
		SYM_LTGT = 24, // '<>'
		SYM_EQEQ = 25, // '=='
		SYM_GT = 26, // '>'
		SYM_GTEQ = 27, // '>='
		SYM_ALIGN = 28, // align
		SYM_BACK = 29, // back
		SYM_BOTTOM = 30, // bottom
		SYM_BOX = 31, // box
		SYM_CENTER = 32, // center
		SYM_CONST = 33, // const
		SYM_DEFER = 34, // defer
		SYM_DEFINE = 35, // define
		SYM_DIVIDE = 36, // divide
		SYM_FACE = 37, // face
		SYM_FLOAT = 38, // Float
		SYM_FRONT = 39, // front
		SYM_FULL = 40, // full
		SYM_IDENTIFIER = 41, // Identifier
		SYM_INCLUDE = 42, // include
		SYM_INSTANCE = 43, // instance
		SYM_INTEGER = 44, // Integer
		SYM_LEFT = 45, // left
		SYM_LOCAL = 46, // local
		SYM_LOFT = 47, // loft
		SYM_MATERIAL = 48, // material
		SYM_MODULE = 49, // module
		SYM_MOVE = 50, // move
		SYM_NGON = 51, // ngon
		SYM_OCCLUDED = 52, // occluded
		SYM_OCCLUDES = 53, // occludes
		SYM_PARTIAL = 54, // partial
		SYM_PRISM = 55, // prism
		SYM_PRISMSIDES = 56, // prismsides
		SYM_REPEAT = 57, // repeat
		SYM_RIGHT = 58, // right
		SYM_ROTATE = 59, // rotate
		SYM_SCALE = 60, // scale
		SYM_SELECT = 61, // select
		SYM_SIDES = 62, // sides
		SYM_SNAP = 63, // snap
		SYM_STRINGLITERAL = 64, // StringLiteral
		SYM_TOP = 65, // top
		SYM_X = 66, // x
		SYM_XY = 67, // xy
		SYM_XZ = 68, // xz
		SYM_Y = 69, // y
		SYM_YZ = 70, // yz
		SYM_Z = 71, // z
		SYM_ADDEXP = 72, // <Add Exp>
		SYM_ASSIGMENT = 73, // <Assigment>
		SYM_ASSIGMENTCONST = 74, // <AssigmentConst>
		SYM_DIRECTION = 75, // <Direction>
		SYM_DIRECTION2D = 76, // <Direction2D>
		SYM_EXPRESSION = 77, // <Expression>
		SYM_FUNCTIONCALL = 78, // <FunctionCall>
		SYM_INCLUDE2 = 79, // <Include>
		SYM_INCLUDELIST = 80, // <IncludeList>
		SYM_INLINEMODULE = 81, // <InlineModule>
		SYM_MODULE2 = 82, // <Module>
		SYM_MODULELIST = 83, // <ModuleList>
		SYM_MODULEREF = 84, // <ModuleRef>
		SYM_MODULEREF1D = 85, // <ModuleRef1D>
		SYM_MODULESELECTOR = 86, // <ModuleSelector>
		SYM_MODULESEQUENCE1D = 87, // <ModuleSequence1D>
		SYM_MULTEXP = 88, // <Mult Exp>
		SYM_NEGATEEXP = 89, // <Negate Exp>
		SYM_OCCLUSIONEXPRESSION = 90, // <OcclusionExpression>
		SYM_OCCLUSIONSUBSCOPE = 91, // <OcclusionSubscope>
		SYM_PARAMETERLIST = 92, // <ParameterList>
		SYM_PERCENTAGE = 93, // <Percentage>
		SYM_PROGRAM = 94, // <Program>
		SYM_RULE = 95, // <Rule>
		SYM_RULEALIGN = 96, // <RuleAlign>
		SYM_RULEBOX = 97, // <RuleBox>
		SYM_RULECALL = 98, // <RuleCall>
		SYM_RULECENTER = 99, // <RuleCenter>
		SYM_RULEDEFER = 100, // <RuleDefer>
		SYM_RULEDIVIDE = 101, // <RuleDivide>
		SYM_RULEFACE = 102, // <RuleFace>
		SYM_RULEINSTANCE = 103, // <RuleInstance>
		SYM_RULELIST = 104, // <RuleList>
		SYM_RULELOCAL = 105, // <RuleLocal>
		SYM_RULELOFT = 106, // <RuleLoft>
		SYM_RULELOFTNGON = 107, // <RuleLoftNGon>
		SYM_RULEMATERIAL = 108, // <RuleMaterial>
		SYM_RULEMOVE = 109, // <RuleMove>
		SYM_RULENGON = 110, // <RuleNGon>
		SYM_RULEOCCLUDE = 111, // <RuleOcclude>
		SYM_RULEPRISM = 112, // <RulePrism>
		SYM_RULEREPEAT = 113, // <RuleRepeat>
		SYM_RULEROTATE = 114, // <RuleRotate>
		SYM_RULESCALE = 115, // <RuleScale>
		SYM_RULESELECT = 116, // <RuleSelect>
		SYM_RULESNAP = 117, // <RuleSnap>
		SYM_SELECTOR = 118, // <Selector>
		SYM_SUBSCOPE = 119, // <Subscope>
		SYM_VALUE = 120, // <Value>
		SYM_VECTOR1D = 121, // <Vector1D>
		SYM_VECTOR2D = 122, // <Vector2D>
		SYM_VECTOR2DLIST = 123, // <Vector2DList>
		SYM_VECTOR3D = 124, // <Vector3D>
		SYM_VECTOR3DLIST = 125, // <Vector3DList>
	};

	enum ProductionConstants
	{
		PROD_EXPRESSION_GT = 0, // <Expression> ::= <Expression> '>' <Add Exp>
		PROD_EXPRESSION_LT = 1, // <Expression> ::= <Expression> '<' <Add Exp>
		PROD_EXPRESSION_LTEQ = 2, // <Expression> ::= <Expression> '<=' <Add Exp>
		PROD_EXPRESSION_GTEQ = 3, // <Expression> ::= <Expression> '>=' <Add Exp>
		PROD_EXPRESSION_EQEQ = 4, // <Expression> ::= <Expression> '==' <Add Exp>
		PROD_EXPRESSION_LTGT = 5, // <Expression> ::= <Expression> '<>' <Add Exp>
		PROD_EXPRESSION = 6, // <Expression> ::= <Add Exp>
		PROD_ADDEXP_PLUS = 7, // <Add Exp> ::= <Add Exp> '+' <Mult Exp>
		PROD_ADDEXP_MINUS = 8, // <Add Exp> ::= <Add Exp> '-' <Mult Exp>
		PROD_ADDEXP = 9, // <Add Exp> ::= <Mult Exp>
		PROD_MULTEXP_TIMES = 10, // <Mult Exp> ::= <Mult Exp> '*' <Negate Exp>
		PROD_MULTEXP_DIV = 11, // <Mult Exp> ::= <Mult Exp> '/' <Negate Exp>
		PROD_MULTEXP = 12, // <Mult Exp> ::= <Negate Exp>
		PROD_NEGATEEXP_MINUS = 13, // <Negate Exp> ::= '-' <Value>
		PROD_NEGATEEXP = 14, // <Negate Exp> ::= <Value>
		PROD_PERCENTAGE_PERCENT = 15, // <Percentage> ::= <Value> '%'
		PROD_FUNCTIONCALL_IDENTIFIER_LPAREN_RPAREN = 16, // <FunctionCall> ::= Identifier '(' <ParameterList> ')'
		PROD_PARAMETERLIST = 17, // <ParameterList> ::= <Expression>
		PROD_PARAMETERLIST_COMMA = 18, // <ParameterList> ::= <Expression> ',' <ParameterList>
		PROD_PARAMETERLIST2 = 19, // <ParameterList> ::= 
		PROD_VALUE_IDENTIFIER = 20, // <Value> ::= Identifier
		PROD_VALUE_LPAREN_RPAREN = 21, // <Value> ::= '(' <Expression> ')'
		PROD_VALUE = 22, // <Value> ::= <Percentage>
		PROD_VALUE2 = 23, // <Value> ::= <FunctionCall>
		PROD_VALUE_INTEGER = 24, // <Value> ::= Integer
		PROD_VALUE_FLOAT = 25, // <Value> ::= Float
		PROD_DIRECTION_X = 26, // <Direction> ::= x
		PROD_DIRECTION_Y = 27, // <Direction> ::= y
		PROD_DIRECTION_Z = 28, // <Direction> ::= z
		PROD_DIRECTION2D_XY = 29, // <Direction2D> ::= xy
		PROD_DIRECTION2D_YZ = 30, // <Direction2D> ::= yz
		PROD_DIRECTION2D_XZ = 31, // <Direction2D> ::= xz
		PROD_VECTOR3D_LBRACKET_COMMA_COMMA_RBRACKET = 32, // <Vector3D> ::= '[' <Expression> ',' <Expression> ',' <Expression> ']'
		PROD_VECTOR2D_LBRACKET_COMMA_RBRACKET = 33, // <Vector2D> ::= '[' <Expression> ',' <Expression> ']'
		PROD_VECTOR1D_LBRACKET_RBRACKET = 34, // <Vector1D> ::= '[' <Expression> ']'
		PROD_VECTOR1D_LBRACKET_PIPE_IDENTIFIER_RBRACKET = 35, // <Vector1D> ::= '[' <Expression> '|' Identifier ']'
		PROD_VECTOR3DLIST = 36, // <Vector3DList> ::= <Vector3D> <Vector3DList>
		PROD_VECTOR3DLIST2 = 37, // <Vector3DList> ::= <Vector3D>
		PROD_VECTOR2DLIST = 38, // <Vector2DList> ::= <Vector2D> <Vector2DList>
		PROD_VECTOR2DLIST2 = 39, // <Vector2DList> ::= <Vector2D>
		PROD_PROGRAM = 40, // <Program> ::= <IncludeList> <ModuleList>
		PROD_INCLUDELIST = 41, // <IncludeList> ::= <Include> <IncludeList>
		PROD_INCLUDELIST2 = 42, // <IncludeList> ::= 
		PROD_INCLUDE_INCLUDE_STRINGLITERAL = 43, // <Include> ::= include StringLiteral
		PROD_MODULELIST = 44, // <ModuleList> ::= <Module> <ModuleList>
		PROD_MODULELIST2 = 45, // <ModuleList> ::= 
		PROD_ASSIGMENT_DEFINE_IDENTIFIER = 46, // <Assigment> ::= define Identifier <Expression>
		PROD_ASSIGMENTCONST_CONST_IDENTIFIER = 47, // <AssigmentConst> ::= const Identifier <Expression>
		PROD_MODULE_IDENTIFIER = 48, // <Module> ::= Identifier <Selector> <InlineModule>
		PROD_SELECTOR_COLON = 49, // <Selector> ::= ':' <ModuleSelector>
		PROD_SELECTOR = 50, // <Selector> ::= 
		PROD_MODULESELECTOR_COMMA = 51, // <ModuleSelector> ::= <Expression> ',' <ModuleSelector>
		PROD_MODULESELECTOR = 52, // <ModuleSelector> ::= <Expression>
		PROD_MODULESELECTOR_COMMA2 = 53, // <ModuleSelector> ::= <OcclusionExpression> ',' <ModuleSelector>
		PROD_MODULESELECTOR2 = 54, // <ModuleSelector> ::= <OcclusionExpression>
		PROD_OCCLUSIONEXPRESSION_OCCLUDED_IDENTIFIER = 55, // <OcclusionExpression> ::= occluded <OcclusionSubscope> Identifier <Expression>
		PROD_OCCLUSIONSUBSCOPE_FULL = 56, // <OcclusionSubscope> ::= full
		PROD_OCCLUSIONSUBSCOPE_PARTIAL = 57, // <OcclusionSubscope> ::= partial
		PROD_OCCLUSIONSUBSCOPE_FRONT = 58, // <OcclusionSubscope> ::= front
		PROD_OCCLUSIONSUBSCOPE_BACK = 59, // <OcclusionSubscope> ::= back
		PROD_OCCLUSIONSUBSCOPE_TOP = 60, // <OcclusionSubscope> ::= top
		PROD_OCCLUSIONSUBSCOPE_BOTTOM = 61, // <OcclusionSubscope> ::= bottom
		PROD_INLINEMODULE_LBRACE_RBRACE = 62, // <InlineModule> ::= '{' <RuleList> '}'
		PROD_MODULEREF_IDENTIFIER = 63, // <ModuleRef> ::= Identifier
		PROD_MODULEREF = 64, // <ModuleRef> ::= <InlineModule>
		PROD_SUBSCOPE = 65, // <Subscope> ::= <Direction2D>
		PROD_SUBSCOPE_SIDES = 66, // <Subscope> ::= sides
		PROD_SUBSCOPE_PRISM = 67, // <Subscope> ::= prism
		PROD_SUBSCOPE_TOP = 68, // <Subscope> ::= top
		PROD_SUBSCOPE_BOTTOM = 69, // <Subscope> ::= bottom
		PROD_SUBSCOPE_FRONT = 70, // <Subscope> ::= front
		PROD_SUBSCOPE_BACK = 71, // <Subscope> ::= back
		PROD_SUBSCOPE_LEFT = 72, // <Subscope> ::= left
		PROD_SUBSCOPE_RIGHT = 73, // <Subscope> ::= right
		PROD_SUBSCOPE_PRISMSIDES = 74, // <Subscope> ::= prismsides
		PROD_RULELIST = 75, // <RuleList> ::= <Rule> <RuleList>
		PROD_RULELIST2 = 76, // <RuleList> ::= 
		PROD_RULE = 77, // <Rule> ::= <RuleCall>
		PROD_RULE2 = 78, // <Rule> ::= <RuleDefer>
		PROD_RULE3 = 79, // <Rule> ::= <Assigment>
		PROD_RULE4 = 80, // <Rule> ::= <AssigmentConst>
		PROD_RULE5 = 81, // <Rule> ::= <RuleInstance>
		PROD_RULE6 = 82, // <Rule> ::= <RuleRepeat>
		PROD_RULE7 = 83, // <Rule> ::= <RuleDivide>
		PROD_RULE8 = 84, // <Rule> ::= <RuleCenter>
		PROD_RULE9 = 85, // <Rule> ::= <RuleMove>
		PROD_RULE10 = 86, // <Rule> ::= <RuleFace>
		PROD_RULE11 = 87, // <Rule> ::= <RuleRotate>
		PROD_RULE12 = 88, // <Rule> ::= <RuleScale>
		PROD_RULE13 = 89, // <Rule> ::= <RuleBox>
		PROD_RULE14 = 90, // <Rule> ::= <RulePrism>
		PROD_RULE15 = 91, // <Rule> ::= <RuleNGon>
		PROD_RULE16 = 92, // <Rule> ::= <RuleLoft>
		PROD_RULE17 = 93, // <Rule> ::= <RuleLoftNGon>
		PROD_RULE18 = 94, // <Rule> ::= <RuleSelect>
		PROD_RULE19 = 95, // <Rule> ::= <RuleOcclude>
		PROD_RULE20 = 96, // <Rule> ::= <RuleLocal>
		PROD_RULE21 = 97, // <Rule> ::= <RuleSnap>
		PROD_RULE22 = 98, // <Rule> ::= <RuleAlign>
		PROD_RULE23 = 99, // <Rule> ::= <RuleMaterial>
		PROD_RULECALL_MODULE_IDENTIFIER = 100, // <RuleCall> ::= module Identifier
		PROD_RULECALL_MODULE = 101, // <RuleCall> ::= module <InlineModule>
		PROD_RULEDEFER_DEFER_IDENTIFIER = 102, // <RuleDefer> ::= defer Identifier
		PROD_RULEINSTANCE_INSTANCE_STRINGLITERAL = 103, // <RuleInstance> ::= instance StringLiteral
		PROD_RULEREPEAT_REPEAT = 104, // <RuleRepeat> ::= repeat <Direction> <ModuleRef1D>
		PROD_RULEDIVIDE_DIVIDE = 105, // <RuleDivide> ::= divide <Direction> <ModuleSequence1D>
		PROD_RULEDIVIDE_DIVIDE_PRISM = 106, // <RuleDivide> ::= divide prism <Expression> <ModuleRef> <ModuleRef> <ModuleRef> <ModuleRef>
		PROD_RULELOFT_LOFT = 107, // <RuleLoft> ::= loft <Vector3DList>
		PROD_RULELOFTNGON_LOFT_NGON = 108, // <RuleLoftNGon> ::= loft ngon <Expression> <Vector2DList>
		PROD_MODULEREF1D = 109, // <ModuleRef1D> ::= <Vector1D> <ModuleRef>
		PROD_MODULESEQUENCE1D = 110, // <ModuleSequence1D> ::= <ModuleRef1D> <ModuleSequence1D>
		PROD_MODULESEQUENCE1D2 = 111, // <ModuleSequence1D> ::= <ModuleRef1D>
		PROD_RULEMOVE_MOVE = 112, // <RuleMove> ::= move <Vector3D>
		PROD_RULEMOVE_MOVE2 = 113, // <RuleMove> ::= move <Direction> <Expression>
		PROD_RULEMOVE_MOVE3 = 114, // <RuleMove> ::= move <Direction2D> <Vector2D>
		PROD_RULEROTATE_ROTATE = 115, // <RuleRotate> ::= rotate <Direction> <Expression>
		PROD_RULEFACE_FACE_IDENTIFIER = 116, // <RuleFace> ::= face Identifier
		PROD_RULESELECT_SELECT = 117, // <RuleSelect> ::= select <Subscope> <ModuleRef>
		PROD_RULEBOX_BOX = 118, // <RuleBox> ::= box
		PROD_RULEPRISM_PRISM = 119, // <RulePrism> ::= prism
		PROD_RULENGON_NGON = 120, // <RuleNGon> ::= ngon <Expression>
		PROD_RULECENTER_CENTER = 121, // <RuleCenter> ::= center
		PROD_RULECENTER_CENTER2 = 122, // <RuleCenter> ::= center <Direction>
		PROD_RULECENTER_CENTER3 = 123, // <RuleCenter> ::= center <Direction2D>
		PROD_RULESCALE_SCALE = 124, // <RuleScale> ::= scale <Direction> <Expression>
		PROD_RULESCALE_SCALE2 = 125, // <RuleScale> ::= scale <Direction2D> <Vector2D>
		PROD_RULESCALE_SCALE3 = 126, // <RuleScale> ::= scale <Vector3D>
		PROD_RULEOCCLUDE_OCCLUDES_IDENTIFIER = 127, // <RuleOcclude> ::= occludes Identifier
		PROD_RULELOCAL_LOCAL_IDENTIFIER = 128, // <RuleLocal> ::= local Identifier
		PROD_RULESNAP_SNAP_IDENTIFIER = 129, // <RuleSnap> ::= snap Identifier <Direction> <Expression>
		PROD_RULEALIGN_ALIGN = 130, // <RuleAlign> ::= align <Direction>
		PROD_RULEMATERIAL_MATERIAL_IDENTIFIER = 131, // <RuleMaterial> ::= material Identifier
	};
}