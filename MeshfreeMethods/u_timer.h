/*! \file u_timer.h */

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

#ifndef _M4_TIMER_
#define _M4_TIMER_

#include <ctime>
#include <chrono>
#include <string>
#include <map>
#pragma warning (disable:4244)
#pragma warning (disable:4996)


namespace msl {

	//#####################################################################
	//# Class BasicTimer
	//#####################################################################
	//! Basic timer type used in the global Timer class 

	class BasicTimer {
	public:
		typedef std::chrono::high_resolution_clock::duration TimeDuration;

		BasicTimer() { start(); };
		bool isStopped() const { return is_stopped_; }
		void start();
		void stop();
		void resume();
		TimeDuration elapsed() const;

	private:
		bool is_stopped_;
		TimeDuration time_;
	};

	//! Function to retrieve current cpu time
	void getCpuTime(BasicTimer::TimeDuration& current);

	//#####################################################################
	//# Class Timer
	//#####################################################################
	//! Timer class for global time management
	// Singleton design pattern implemented for class Timer,
	// manages a set of TimeKeeper objects.

	class Timer {

	public:

		//! Destructor. 
		~Timer() {}

		//! Singeton.
		static Timer& instance();

		//! Print out current date time
		static std::string getCurrentTime();

		//! Print out local day time, in the format "hh:mm:ss"
		static std::string getDayTime();

		//! Print out local date
		static std::string getDate();

		//! Time table to represent the time duration.
		struct TimeTable {
			std::chrono::hours hr_;
			std::chrono::minutes min_;
			std::chrono::seconds sec_;
			std::chrono::milliseconds msec_;
		};

		//! Starts specified timer, creates the timer if it does not exist.
		void startTimer(const std::string& name);

		//! Reset a stopped timer.
		void resetTimer(const std::string& name) { timers_[name].start(); }

		//! Stops specified timer.
		void stopTimer(const std::string& name) { timers_[name].stop(); }

		//! Get elapsed time data of specified timer 
		// TimeTable (hour, min, sec, millisec).
		TimeTable getElapsed(const std::string& name);

		std::string format(const std::string& name);

		std::string stopAndFormat(const std::string& name) { timers_[name].stop(); return format(name); }

		//! @name Private and unimplemented to prevent use
		Timer(const Timer&) = delete;

		Timer& operator=(const Timer&) = delete;

	protected:

		//! Map that associates a name with a TimeKeeper.
		std::map<std::string, BasicTimer> timers_;

	private:

		//! Private constructor, prevent from use.
		Timer() {}
	};
}

#endif // ! _M4_TIMER_