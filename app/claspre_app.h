// 
// Copyright (c) 2012, Benjamin Kaufmann
// 
// This file is part of Claspre. See http://potassco.sourceforge.net/
// 
// Claspre is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Claspre is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Claspre; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
#ifndef CLASPRE_APP_H_INCLUDED
#define CLASPRE_APP_H_INCLUDED

#ifdef _MSC_VER
#pragma warning (disable : 4200) // nonstandard extension used : zero-sized array
#pragma once
#endif
#include <program_opts/app_options.h>
#include <clasp/clasp_facade.h>
#include <clasp/util/timer.h>
#include <clasp/util/platform.h>
#include <string>
#include <vector>
#include <iosfwd>
#include <memory>
#include <stdio.h>
#include <signal.h>
/////////////////////////////////////////////////////////////////////////////////////////
// Output macros and app exit codes
/////////////////////////////////////////////////////////////////////////////////////////
#define WRITE_STDERR(type,sys,msg) ( fflush(stdout), fprintf(stderr, "*** %-5s: (%s): %s\n", (type),(sys),(msg)), fflush(stderr) )
#define ERROR_OUT(sys,msg)   WRITE_STDERR("ERROR", (sys), (msg))
#define INFO_OUT(sys,msg)    WRITE_STDERR("Info", (sys), (msg))
#define WARNING_OUT(sys,msg) WRITE_STDERR("Warn", (sys), (msg))

// exit codes
#define S_SATISFIABLE   10  // problem is satisfiable
#define S_UNSATISFIABLE 20  // problem was proved to be unsatisfiable
#define S_UNKNOWN       0   // satisfiablity of problem not knwon; search was interrupted or did not start
#define S_ERROR EXIT_FAILURE// internal error, except out of memory
#define S_MEMORY        127 // out of memory
namespace Claspre {
using Clasp::ClaspFacade;
using Clasp::Solver;
using Clasp::Timer;
using Clasp::ProcessTime;
using Clasp::RealTime;
/////////////////////////////////////////////////////////////////////////////////////////
// Claspre::Options
/////////////////////////////////////////////////////////////////////////////////////////
struct Options {
	Options();
	typedef ProgramOptions::StringSeq StringSeq;
	typedef std::pair<int, int>       SolveLimit;
	void initOptions(ProgramOptions::OptionContext& root);
	bool validateOptions(const ProgramOptions::ParsedOptions& parsed, ProgramOptions::Messages&);
	void applyDefaults(Clasp::Input::Format f);
	static bool parsePositional(const std::string& s, std::string& out);
	Clasp::ClaspConfig clasp;
	StringSeq          input;     // list of input files - only first used!
	SolveLimit         limit;     // current solve limit
	uint32             timeout;   // timeout in seconds (default: none=-1)
	bool               fastExit;  // force fast exit (no dtors)
};

/////////////////////////////////////////////////////////////////////////////////////////
// Claspre::Output
/////////////////////////////////////////////////////////////////////////////////////////
class Output {
public:
	struct Result {
		double totalTime;
		double solveTime;
		double modelTime;
		double unsatTime;
		double cpuTime;
		int    signal;
		bool   complete;
	};
	explicit Output(const Options& o, const ClaspFacade& f)
		: opts_(o)
		, facade_(f) {}
	// called once solving has finished
	void onExit(const Result& r);
	// called on each state change (read->preprocess->solve)
	void onState(bool enter);
	// called once the solver is ready to solve
	void onProgramPrepared(const Solver& s);
	// called on *before* each (re-)start
	void onRestart(const Solver& s, uint64 conflictLimit, uint32 learntLimit);
	// called *after* each deletion operation
	void onDeletion(const Solver& s, uint64 conflictLimit, uint32 learntLimit);
	// called on each model
	void onModel(const Solver& s, const Clasp::Enumerator& en);
private:
	Output& operator=(const Output&);
	typedef Clasp::Input::Format InputFormat;
	const Options&     opts_;
	const ClaspFacade& facade_;
	// <TBD>
	// add more data here
};
/////////////////////////////////////////////////////////////////////////////////////////
// Claspre::Application
/////////////////////////////////////////////////////////////////////////////////////////
#define CLASPRE_VERSION "2.0"
class Application : public ProgramOptions::AppOptions, public ClaspFacade::Callback, public Clasp::ProgressReport {
public:
	static Application& instance();    // returns the singleton instance
	static void sigHandler(int sig);   // signal/timeout handler
	void   installSigHandlers();       // adds handlers for SIGINT, SIGTERM, SIGALRM
	int    run(int argc, char** argv); // "entry-point"
	void   printWarnings()const;
private:
	Application();
	Application(const Application&);
	const Application& operator=(const Application&);
	// -------------------------------------------------------------------------------------------
	// AppOptions interface
	void    printHelp(const ProgramOptions::OptionContext& root)    ;
	void    printVersion(const ProgramOptions::OptionContext& root) ;
	void    initOptions(ProgramOptions::OptionContext& root) {
		app_.initOptions(root);
	}
	bool    validateOptions(const ProgramOptions::OptionContext&, const ProgramOptions::ParsedOptions& vm, ProgramOptions::Messages& m) {
		return app_.validateOptions(vm, m);
	}
	// -------------------------------------------------------------------------------------------
	// ClaspFacade::Callback interface
	void state(ClaspFacade::Event e, ClaspFacade& f);
	void event(const Solver& s, ClaspFacade::Event e, ClaspFacade& f);
	void warning(const char* msg) { messages.warning.push_back(msg); }
	// ProgressReport interface
	void reportProgress(const Clasp::SolveEvent& ev);
	void reportProgress(const Clasp::PreprocessEvent& ev);
	// -------------------------------------------------------------------------------------------
	std::istream& getStream();
	void killAlarm();
	void kill(int sig);
	int  blockSignals();
	void unblockSignals(bool deliverPending);
	int  exception(int status, const char* what);
	void appTerminate(int exitCode) const;
	// -------------------------------------------------------------------------------------------  
	// Status information & output
	void printResult(int sig);
	// -------------------------------------------------------------------------------------------  
	Options                app_;
	Timer<ProcessTime>     cpuTotalTime_;
	Timer<RealTime>        timer_[ClaspFacade::num_states]; // one for each state
	double                 timeToFirst_;                    // time to first model
	double                 timeToLast_;                     // time to last model
	ClaspFacade*           facade_;
	std::auto_ptr<Output>  out_;
	volatile sig_atomic_t  blocked_;
	volatile sig_atomic_t  pending_;
};

const char* const app_banner = "claspre version " CLASPRE_VERSION " (based on clasp " CLASP_VERSION ")";

}
#endif
