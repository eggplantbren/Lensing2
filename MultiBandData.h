#ifndef _MultiBandData_
#define _MultiBandData_

#include "Data.h"
#include <vector>

namespace Lensing2
{

class MultiBandData
{
	private:
		std::vector<Data> bands;

	public:
		MultiBandData();

		static MultiBandData instance;
};

} // namespace Lensing2

#endif

