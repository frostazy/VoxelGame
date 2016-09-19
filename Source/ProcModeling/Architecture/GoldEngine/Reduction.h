#ifndef GOLDCPP_REDUCTION_H
#define GOLDCPP_REDUCTION_H

#include "Token.h"

namespace GoldCPP
{
  struct Production;

  struct Reduction
  {
    TokenList Branches;
    Production *Parent;
    void *User;

    Reduction(size_t n) :
      Branches(n, NULL),
      Parent(NULL),
      User(NULL)
    {}

	int Count()
	{
		return (int)Branches.Count();
	}

	Token& get_Data(int idx)
	{
		return *(Branches[idx]);
	}
  };
}

#endif // GOLDCPP_REDUCTION_H





