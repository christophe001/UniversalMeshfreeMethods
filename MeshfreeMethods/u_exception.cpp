/*! \file u_exception.cpp */

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


#include "u_exception.h"

namespace msl {
	m_exception::m_exception(const std::string& scope_name, const std::string& class_name,
		const std::string& method, const std::string& message,
		const std::string& file) :
		scope_name_(scope_name), class_name_(class_name), method_(method),
		message_(message), file_(file) {}

	//================================================================================

	std::string m_exception::Info() const {
		std::string msg;
		msg += "\n******************** Exception occurred ********************\n";
		msg = msg + "# Scope: " + scope_name_ + "   Class: " + class_name_;
		msg = msg + "   Method: " + method_ + "\n";
		if (!message_.empty())
			msg = msg + "# Message:\n\t" + message_ + "\n";
		if (!file_.empty())
			msg = msg + "# File:\n\t" + file_ + "\n";
		msg += "************************************************************\n";
		return msg;
	}

}