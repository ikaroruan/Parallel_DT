#ifndef LOCK_CPP_H
#define LOCK_CPP_H

#include <mutex>
#include <thread>

extern std::mutex critical_mutex;

class Lock
{
	public:

	Lock();
	
	bool lock();
	bool try_lock();
	bool unlock();
	bool locked();
	
	//std::thread::id owner();
	//void owner(std::thread::id value);
	//bool is_locked();
	//void is_locked(bool value);

	private:
	std::mutex _mtx;
	std::thread::id _owner;
	bool _locked;
};

#endif
