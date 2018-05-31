#include "Lock.h"

Lock::Lock()
{
	_locked = false;
}

bool Lock::lock()
{
	// Check if the thread already owns the lock.
	if(_locked && _owner == std::this_thread::get_id()){
		return true;
	}

	_mtx.lock();
	_locked = true;
	_owner = std::this_thread::get_id();
	return true;
}

bool Lock::try_lock()
{
	// Check if the thread already owns the lock.
	if(_locked && _owner == std::this_thread::get_id()){
		return true;
	}

	if(_mtx.try_lock()){
		_owner = std::this_thread::get_id();
		_locked = true;
		return true;
	}

	return false;
}

bool Lock::unlock()
{
	// Only the owner thread can unlock the mutex.
	if(_owner == std::this_thread::get_id()){
		_mtx.unlock();
		_locked = false;
		return true;
	}

	return false;
}

bool Lock::locked()
{
	return _locked;
}
