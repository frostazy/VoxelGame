#pragma once

#include "Parser.h"

#include <string>
#include <vector>

namespace GoldCPP {

	struct SourceFile
	{
		std::string filename;
		std::string modSrc;
	};

	struct GrammarError
	{
		std::string filename;
		int line;
		int column;
		std::string message;
	};

	class ModCompiler
	{
	public:
		ModCompiler(const char* inEgtPath = NULL)
			: isEgtLoaded(false)
		{
			if (inEgtPath)
				egtPath = inEgtPath;
		}

		bool Compile(std::vector<SourceFile>& files, std::string& outCode, std::string& outDbg, std::vector<GrammarError>& errors);

	private:
		void LoadEgtIfNeeded();

		Parser parser;
		bool isEgtLoaded;
		std::string egtPath;
	};

}