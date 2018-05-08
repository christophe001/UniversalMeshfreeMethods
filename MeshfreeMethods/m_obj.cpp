/*! \file m_obj.cpp */

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

#include "m_obj.h"
#include "u_exception.h"

namespace msl {

	void MsObj::throwException(const std::string& method, const std::string& message,
		const std::string& file) const {
		throw m_exception(scope_name_, class_name_, method, message, file);
	}

	// ===================================================================================

	std::string MsObj::getException(const std::string& method, const std::string& message,
		const std::string& file) const {
		return m_exception(scope_name_, class_name_, method, message, file).Info();
	}

	// ====================================================================================

	std::string MsObj::exceptionStr(const std::string& scope, const std::string& class_name,
		const std::string& method, const std::string& message,
		const std::string& file) {
		std::string msg;
		msg += "\n******************** Exception occurred ********************\n";
		msg = msg + "# Scope: " + scope + "   Class: " + class_name;
		msg = msg + "   Method: " + method + "\n";
		if (!message.empty())
			msg = msg + "# Message:\n\t" + message + "\n";
		if (!file.empty())
			msg = msg + "# File:\n\t" + file + "\n";
		msg += "************************************************************\n";
		return msg;
	}
}