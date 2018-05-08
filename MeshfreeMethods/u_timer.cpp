/*! \file u_timer.cpp */

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

#include "u_timer.h"
#include <sstream>
#include <iomanip> 
#include <ctime>

using namespace std;
using namespace std::chrono;


namespace msl {
	//! Utilities get_cpu_time --------------------------------------------
	void getCpuTime(BasicTimer::TimeDuration& curr) {
		curr = high_resolution_clock::now().time_since_epoch();
	}
	//! BasicTimer --------------------------------------------------------
	void BasicTimer::start() {
		is_stopped_ = false;
		getCpuTime(time_);
	}

	void BasicTimer::stop() {
		if (isStopped())
			return;
		is_stopped_ = true;
		TimeDuration curr;
		getCpuTime(curr);
		time_ = curr - time_;
	}

	void BasicTimer::resume() {
		if (isStopped()) {
			TimeDuration curr(time_);
			start();
			time_ -= curr;
		}
	}

	BasicTimer::TimeDuration BasicTimer::elapsed() const {
		if (isStopped())
			return time_;
		TimeDuration current;
		getCpuTime(current);
		current -= time_;
		return current;
	}


	//! Class Timer used within application -------------------------------
	Timer& Timer::instance() {
		static Timer timer;
		return timer;
	}

	//! get_current_time, return current time in a string format.
	string Timer::getCurrentTime() {
		time_point<system_clock> tp = system_clock::now();
		time_t curr = system_clock::to_time_t(tp);
		string current_time(ctime(&curr));
		return current_time;
	}

	//! get_day_time, return day time in a string format.
	//! copy day_time part from get_current_time().
	string Timer::getDayTime() {
		string curr_t = getCurrentTime();
		string day_t(curr_t.begin() + 11, curr_t.begin() + 19);
		return day_t;
	}

	std::string Timer::getDate() {
		string curr_t = getCurrentTime();
		string date(curr_t.begin(), curr_t.begin() + 11);
		string year(curr_t.begin() + 19, curr_t.end());
		return date + year;
	}

	//! start_timer
	void Timer::startTimer(const string& name) {
		if (timers_.find(name) == timers_.end())
			timers_[name].start();
		else timers_[name].resume();
	}

	//! get_elapsed
	Timer::TimeTable Timer::getElapsed(const string& name) {
		BasicTimer::TimeDuration td;
		td = timers_[name].elapsed();
		Timer::TimeTable tt;
		tt.hr_ = duration_cast<hours>(td);
		tt.min_ = duration_cast<minutes>(td % hours(1));
		tt.sec_ = duration_cast<seconds>(td % minutes(1));
		tt.msec_ = duration_cast<milliseconds>(td % seconds(1));
		return tt;
	}

	string Timer::format(const string& name) {
		Timer::TimeTable tt = Timer::getElapsed(name);
		stringstream ss;
		ss << tt.hr_.count() << setfill('0')
			<< ":" << setw(2) << tt.min_.count()
			<< ":" << setw(2) << tt.sec_.count()
			<< "." << setw(3) << tt.msec_.count();
		return ss.str();
	}
}