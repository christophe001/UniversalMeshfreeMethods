/*! \file u_exception.h */

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

#ifndef _M4_EXCEPTION_
#define _M4_EXCEPTION_

#include <string>
#pragma warning (disable:4244)
#pragma warning (disable:4267)

namespace msl {

	//#####################################################################
	//# Class m_exception
	//#####################################################################
	//! extend std::exception for exception handling 

	class m_exception : public std::exception {

	protected:

		std::string scope_name_;       //! Name of the scope of the class that generates the exception.

		std::string class_name_;       //! Name of the class that generates the exception.

		std::string method_;           //! Name of the method that generates the exception.

		std::string message_;          //! Message of the exception.

		std::string file_;             //! File related to the exception.

	public:

		m_exception& operator= (const m_exception&) = delete;

		m_exception(const std::string& scope_name, const std::string& class_name,
			const std::string& method, const std::string& message,
			const std::string& file);

		virtual ~m_exception() throw() {}

		std::string Info() const;

		virtual const char* what() const throw() {
			static char buffer[1024];
			std::string str = Info();
			for (size_t i = 0; i < str.size(); i++)
				buffer[i] = str[i];
			buffer[str.size()] = '\0';
			return buffer;
		}
	};
}

#endif // !_M4_EXCEPTION_