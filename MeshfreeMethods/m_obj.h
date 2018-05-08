/*! \file m_obj.h */

//@HEADER
// ************************************************************************
//
//                         M4_Technologies
//                 Copyright (2017) by Wentao Xu
//
// M4 stands for Multiscale Mesh-based Meshfree Method
// Any questions, please contact:
// Wentao Xu   wx2151@columbia.edu
//
// ************************************************************************
//@HEADER

#ifndef _M4_OBJ_
#define _M4_OBJ_

#include <string>
#pragma warning (disable:4244)

namespace msl {

	class MsObj {
	protected:

		std::string scope_name_;

		std::string class_name_;

		void throwException(const std::string& method, const std::string& message,
			const std::string& file = "") const;

		std::string getException(const std::string& method, const std::string& message,
			const std::string& file = "") const;

	public:
		std::string getClassName() const {
			return class_name_;
		}

		static std::string exceptionStr(const std::string& scope, const std::string& class_name,
			const std::string& method, const std::string& message,
			const std::string& file = "");

		MsObj() : scope_name_("msl"), class_name_("MsObj") {}
		virtual ~MsObj() {}
	};
}

#endif // !_M4_Obj_